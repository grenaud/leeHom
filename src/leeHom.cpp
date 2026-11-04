#include <iostream>
#include <fstream>
#include <thread>
#include <mutex>
#include <vector>
#include <string>
#include <cstdio>
#include <cstring>
#include <time.h>

#include <zlib.h>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include "MergeTrimReads.h"
#include "PutProgramInHeaderHTS.h"
#include "FastQParser.h"
#include "BamRecord.h"
#include "Pipeline.h"

#include "libgab.h"
#include "version.h"

using namespace std;

/*
  Overview:

  The main thread reads the input (BAM or fastq) and cuts it into chunks. A chunk is a few thousand reads/clusters.
  Several worker threads process the chunks, each one with its own MergeTrimReads object, and a writer
  thread writes the results. The chunks are written in the order in which they were read (see OrderedPipeline),
  hence the output is in the same order as the input whatever the number of threads.

  For BAM the compression/decompression is done by htslib's thread pool; for fastq the workers compress
  the text of their chunk as an independent gzip member (a concatenation of gzip members is a valid gzip file).
*/

// https://samtools.github.io/hts-specs/SAMtags.pdf
// RX:Z:sequence+ Sequence bases from the unique molecular identifier. These could be either corrected or
// uncorrected. Unlike MI, the value may be non-unique in the file. Should be comprised of a sequence of
// bases. In the case of multiple unique molecular identifiers (e.g., one on each end of the template) the
// recommended implementation concatenates all the barcodes with a hyphen (-) between the different
// barcodes.
// If the bases represent corrected bases, the original sequence can be stored in OX (similar to OQ storing
// the original qualities of bases.)
//
// QX:Z:qualities+ Phred quality of the unique molecular identifier sequence in the RX tag. Same encoding
// as QUAL, i.e., Phred score + 33. The qualities here may have been corrected (Raw bases and qualities
// can be stored in OX and BZ respectively.) The lengths of the QX and the RX tags must match. In the
// case of multiple unique molecular identifiers (e.g., one on each end of the template) the recommended
// implementation concatenates all the quality strings with a space (' ') between the different strings.

MergeTrimReads * mtr;
string options_adapter_F_BAM="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACIIIIIIIATCTCGTATGCCGTCTTCTGCTTG";
string options_adapter_S_BAM="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT";
string options_adapter_chimera_BAM="ACACTCTTTCCCTACACGTCTGAACTCCAG,ACACTCTTTCCCACACGTCTGAACTCCAGT,ACACTCTTTCCCTACACACGTCTGAACTCC,CTCTTTCCCTACACGTCTGAACTCCAGTCA,GAAGAGCACACGTCTGAACTCCAGTCACII,GAGCACACGTCTGAACTCCAGTCACIIIII,GATCGGAAGAGCACACGTCTGAACTCCAGT,AGATCGGAAGAGCACACGTCTGAACTCCAG,AGAGCACACGTCTGAACTCCAGTCACIIII,ACACGTCTGAACTCCAGTCACIIIIIIIAT,GTGCACACGTCTGAACTCCAGTCACIIIII,AGCACACGTCTGAACTCCAGTCACIIIIII,CGTATGCCGTCTTCTGCTTGAAAAAAAAAA";

unsigned int sizeChunkDefault = 5000;
unsigned int sizeChunk        = 5000;
int          numberOfThreads  = 1;


//! Everything the run functions need to know, filled by main()
struct Settings {
    int          gzLevel;     //compression level of the output, -1: default
    bool         keepOrig;
    unsigned int umif;
    unsigned int umir;
    unsigned int umiTf;
    unsigned int umiTr;
    bool         detectAdapt;
    bool         trimKey;
    string       key1;
    string       key2;
    bool         printLog;
    string       logFileName;

    //BAM
    bool         allowAligned;
    bool         produceUnCompressedBAM;
    string       bamFile;
    string       bamFileOUT;
    bool         fastqToBam;        //fastq input, BAM output
    int          fqQualOffset;      //offset of the qualities of the fastq input

    //fastq
    bool         singleEndModeFQ;
    string       fastqfile1;
    string       fastqfile2;
    string       fastqoutfile;
    bool         fastqoutfileb;
    string       fastqoutfilese;
    string       fastqoutfilesef;
    string       fastqoutfilepe1;
    string       fastqoutfilepe1f;
    string       fastqoutfilepe2;
    string       fastqoutfilepe2f;
    string       intfastqoutfile;
    //the output stream (single, pair 1, pair 2, failed single, failed pair 1, failed pair 2) -> file
    int          sidx[6];
};

enum { S_SINGLE=0, S_PAIR1=1, S_PAIR2=2, S_SINGLEF=3, S_PAIR1F=4, S_PAIR2F=5 };


//! A chunk of BAM records, the output chunks have the same type
struct BamChunk {
    uint64_t         firstRecord;     //index in the input of the first record
    vector<bam1_t *> recs;

    BamChunk():firstRecord(0){}
    ~BamChunk(){
	for(size_t i=0;i<recs.size();i++)
	    bam_destroy1(recs[i]);
    }
private:
    BamChunk(const BamChunk &);
    BamChunk & operator=(const BamChunk &);
};

//! A chunk of fastq records
struct FqChunk {
    uint64_t         firstRecord;
    vector<fqrecord> recs;

    FqChunk():firstRecord(0){}
};

//! The text to write to each of the output streams for a chunk, gzip compressed
struct FqOutChunk {
    string data[6];
};


//! Compresses in as a complete gzip member (header and trailer)
static void gzipMember(const string & in, string & out, int level){
    z_stream zs;
    memset(&zs, 0, sizeof(zs));
    if(deflateInit2(&zs, level, Z_DEFLATED, 15+16, 8, Z_DEFAULT_STRATEGY) != Z_OK){
	cerr<<"Unable to initialize zlib"<<endl;
	exit(1);
    }
    out.resize(deflateBound(&zs, in.size()));
    zs.next_in   = (Bytef *) in.data();
    zs.avail_in  = in.size();
    zs.next_out  = (Bytef *) &out[0];
    zs.avail_out = out.size();
    if(deflate(&zs, Z_FINISH) != Z_STREAM_END){
	cerr<<"Unable to compress the output"<<endl;
	exit(1);
    }
    out.resize(zs.total_out);
    deflateEnd(&zs);
}


//! CPU time used so far by the calling thread, in seconds
static double threadCpuSeconds(){
    struct timespec ts;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts);
    return double(ts.tv_sec) + double(ts.tv_nsec)*1e-9;
}

//! With LEEHOM_STATS set, the CPU time of each stage is reported at the end
static mutex  statsMutex;
static double statsReader = 0.0, statsWorkers = 0.0, statsWriter = 0.0;
static bool   statsEnabled(){ static const bool b = (getenv("LEEHOM_STATS") != NULL); return b; }
static void   statsAdd(double & where, double seconds){ lock_guard<mutex> lock(statsMutex); where += seconds; }
static void   statsReport(){
    if(statsEnabled())
	cerr<<"CPU seconds: reading "<<statsReader<<", processing (all workers) "<<statsWorkers<<", writing "<<statsWriter<<endl;
}


static inline void fqPut(string & b, const string & head, const string & cmt, const string & seq, const string & qual){
    b.append(head);
    b.append(cmt);
    b.push_back('\n');
    b.append(seq);
    b.append("\n+\n");
    b.append(qual);
    b.push_back('\n');
}


bool setPairedDetectAdapters(vector<string>         & adapterSeqs_fwd,
			     vector<string> & qualadaptSeq_fwd,
			     vector<string>         & adapterSeqs_rev,
			     vector<string> & qualadaptSeq_rev){

    //collecting the adapters
    string fwd_adapt=mtr->mlConsensus(adapterSeqs_fwd,qualadaptSeq_fwd);
    
    float fwdIfrac;
    int   fwdI=0;
    for(unsigned int i=0;i<fwd_adapt.size();i++){	
	if(fwd_adapt[i] == 'I'){
	    fwdI++;
	}
    }
    fwdIfrac = (float( fwdI) /float( fwd_adapt.size() ) );
    float thresholdUndefined = 0.2;

    bool errorFound=false;
    if(fwdIfrac < thresholdUndefined){
	cerr<<"Inferred fwd adapter: "<<fwd_adapt<<endl;
	mtr->setFadapter(fwd_adapt);
    }else{
	cerr<<"Cannot infer the forward adapter, only found "<<fwd_adapt<<" fraction undefined: "<<(fwdIfrac*100)<<"%"<<" threshold: "<<(thresholdUndefined*100)<<"%"<<endl;
	errorFound=true;
	//return false;
    }
    
    string rev_adapt=mtr->mlConsensus(adapterSeqs_rev,qualadaptSeq_rev);

    float revIfrac;
    int   revI=0;
    for(unsigned int i=0;i<rev_adapt.size();i++){	
	if(rev_adapt[i] == 'I'){
	    revI++;
	}
    }
    revIfrac = (float( revI) /float( rev_adapt.size() ) );
    if(revIfrac < thresholdUndefined){
	cerr<<"Inferred rev adapter: "<<rev_adapt<<endl;
	mtr->setRadapter(rev_adapt);		
    }else{
	cerr<<"Cannot infer the reverse adapter, only found "<<rev_adapt<<" fraction undefined: "<<(revIfrac*100)<<"%"<<" threshold: "<<(thresholdUndefined*100)<<"%"<<endl;
	errorFound=true;
	//return false;
    }

    if(errorFound)
	return false;
    
    return true;

}


bool detectAdaptBAM(const BamChunk * currentChunk){

    BamRecord al;
    BamRecord al2;
    bool al2Null=true;
    vector<string> adapterSeqs_fwd;
    vector<string> qualadaptSeq_fwd;
    vector<string> adapterSeqs_rev;
    vector<string> qualadaptSeq_rev;

    //iterating over each record in the chunk
    for(unsigned i=0;i<currentChunk->recs.size();i++){
	
	al.fromBam(currentChunk->recs[i]);
	
	if(al.IsPaired() && 
	   al2Null ){
	    al2=al;
	    al2Null=false;
	    continue;
	}else{
	    if(al.IsPaired() && 
	       !al2Null){
				
		//bool  result =  mtr->processPair(al,al2);
		mtr->inferadaptPairBAM(al ,al2,adapterSeqs_fwd,qualadaptSeq_fwd,adapterSeqs_rev,qualadaptSeq_rev);
		
	    }else{ 

		//mtr->processSingle(al);
		cerr<<"The automatic detection of the adapters is not yet coded for single-end reads"<<endl;
		exit(1);
	    } //end single end
	    al2Null=true;
	}//second pair		    

    }

    return setPairedDetectAdapters(adapterSeqs_fwd,qualadaptSeq_fwd,adapterSeqs_rev,qualadaptSeq_rev);
}


bool detectAdaptFQ(const FqChunk * currentChunk){
    fqrecord fr;

    vector<string> adapterSeqs_fwd;
    vector<string> qualadaptSeq_fwd;
    vector<string> adapterSeqs_rev;
    vector<string> qualadaptSeq_rev;

    //iterating over each record in the chunk
    for(unsigned i=0;i<currentChunk->recs.size();i++){
	
	fr = currentChunk->recs[i];
	
	if( fr.paired ){
	    mtr->inferadaptPairFQ(fr,adapterSeqs_fwd,qualadaptSeq_fwd,adapterSeqs_rev,qualadaptSeq_rev);
	}else{
	    //mtr->processSingle(al);
	    cerr<<"The automatic detection of the adapters is not yet coded for single-end reads"<<endl;
	    exit(1);
	}//second pair
		    

    }

    return setPairedDetectAdapters(adapterSeqs_fwd,qualadaptSeq_fwd,adapterSeqs_rev,qualadaptSeq_rev);
}


//! Processes a chunk of BAM records
static BamChunk * processBamChunk(MergeTrimReads * m, const Settings & S, const BamChunk & in){

    const bool         keepOrig = S.keepOrig;
    const unsigned int umif     = S.umif;
    const unsigned int umir     = S.umir;
    const unsigned int umtf     = S.umiTf;
    const unsigned int umtr     = S.umiTr;
    const bool         umi      = ((umif+umir)>0);

    BamChunk * chunkToWrite = new BamChunk();
    chunkToWrite->firstRecord = in.firstRecord;
    chunkToWrite->recs.reserve(in.recs.size() + in.recs.size()/4);

    BamRecord al;
    BamRecord al2;
    BamRecord orig;
    BamRecord orig2;
    bool al2Null=true;
    string umifS;
    string umirS;
    string umifQ;
    string umirQ;

    //adds a record to the output
    #define emit(rec) { bam1_t * bout = bam_init1(); (rec).toBam(bout); chunkToWrite->recs.push_back(bout); }

    for(unsigned i=0;i<in.recs.size();i++){
	//nonsense work
	// int countIdentical=0;
	// for(unsigned int j=0;j<currentChunk->dataToProcess[i].QueryBases.size();j++){
	//     for(unsigned int k=0;k<currentChunk->dataToProcess[i].QueryBases.size();k++){
	// 	for(unsigned int l=0;l<currentChunk->dataToProcess[i].QueryBases.size();l++){
	// 	    if(currentChunk->dataToProcess[i].QueryBases[j] == currentChunk->dataToProcess[i].QueryBases[k]){
	// 		countIdentical++;
	// 	    }
	// 	}
	//     }
	// }
	
	al.fromBam(in.recs[i]);
	//cout<<rankThread<<" al  "<<al.Name<<" "<<al.IsFirstMate()<<endl;
	//if(!al2Null)
	//cout<<rankThread<<" al2 "<<al2.Name<<" "<<al2.IsFirstMate()<<endl;
	
	if(al.IsPaired() && 
	   al2Null ){
	    al2=al;
	    al2Null=false;
	    continue;
	}else{
	    if(al.IsPaired() && 
	       !al2Null){
				
		if(keepOrig){
		    orig =al;
		}

		//fastq input: the same random choices as when writing fastq, where the index is that of the pair
		m->setRandomContext( S.fastqToBam ? ((in.firstRecord + i)/2) : (in.firstRecord + i) );
		bool  result =  m->processPair(al,al2);

		if( result ){//was merged

		    if(keepOrig){
			orig2 = al2;
			al.Name = al.Name+"_M";
		    }

		    if( umi ){

			if( al.QueryBases.length() > (umif+umir+umtf+umtr)){
			    if(umif != 0){	    
				umifS         = al.QueryBases.substr(0,umif);
				umifQ         = al.Qualities.substr( 0,umif);

				al.QueryBases = al.QueryBases.substr(umif+umtf);
				al.Qualities  = al.Qualities.substr( umif+umtf);

				if(umir == 0){	//no reverse UMI, add the tags now
				    if(!al.AddTag("RX","Z",umifS)){ cerr << "Error while editing RX tags for UMI:"<<umifS<<"#"<< endl; exit(1);   }				    
				    if(!al.AddTag("QX","Z",umifQ)){ cerr << "Error while editing QX tags for UMI:"<<umifQ<<"#"<< endl; exit(1);   }
				}

			    }

			    if(umir != 0){
				umirS         = al.QueryBases.substr(al.QueryBases.length()-umir,  umir);
				umirQ         = al.Qualities.substr( al.Qualities.length() -umir,  umir);
			    
				al.QueryBases = al.QueryBases.substr(0,al.QueryBases.length()-umir-umtr);
				al.Qualities  = al.Qualities.substr(0, al.Qualities.length() -umir-umtr);

				if(umif == 0){	//no forward UMI, just add the reverse
				    				    
				    if(!al.AddTag("RX","Z",umirS)){ cerr << "Error while editing RX tags for UMI:"<<umirS<<"#"<< endl;	exit(1);   }				    
				    if(!al.AddTag("QX","Z",umirQ)){ cerr << "Error while editing QX tags for UMI:"<<umirQ<<"#"<< endl;	exit(1);   }

				}else{
				    //add both umif and umir
				    				    				    
				    if(!al.AddTag("RX","Z",umifS+"-"+umirS)){ cerr << "Error while editing RX tags for UMI:"<<umifS<<"-"<<umirS<<"#"<< endl; exit(1);    }				    
				    if(!al.AddTag("QX","Z",umifQ+" "+umirQ)){ cerr << "Error while editing QX tags for UMI:"<<umirQ<<" "<<umirQ<<"#"<< endl; exit(1);    }
				    
				}
			    }
			}else{//the new fragment is shorter than both umis
			    //leave but fail QC
			    al.SetIsDuplicate(true);
			}
		    }
		    
		    emit(al);
		    //writer.SaveAlignment(al);
		    
		    if(keepOrig){
			orig.SetIsDuplicate(true);
			orig2.SetIsDuplicate(true);
			// writer.SaveAlignment(orig2);
			// writer.SaveAlignment(orig);
			emit(orig2);
			emit(orig);
		    }

		    //the second record is empty
		}else{
		    //keep the sequences as pairs
		    
		    if( al.QueryBases.length()  < (umif+umtf)){ cerr << "Error fwd read length is inferior to the sum of the UMIs"<< endl; exit(1); }
		    if( al2.QueryBases.length() < (umir+umtr)){ cerr << "Error rev read length is inferior to the sum of the UMIs"<< endl; exit(1); }

		    if(umif != 0){	    
			//al2 is the fwd read
			umifS          = al2.QueryBases.substr(0,umif);
			umifQ          = al2.Qualities.substr( 0,umif);
			
			al2.QueryBases = al2.QueryBases.substr(umif+umtf);
			al2.Qualities  = al2.Qualities.substr( umif+umtf);
			
			if(umir == 0){	//no reverse UMI, add the tags now
			    if(!al.AddTag("RX","Z",umifS)){  cerr << "Error while editing RX tags for UMI:"<<umifS<<"#"<< endl; exit(1);   }				    
			    if(!al.AddTag("QX","Z",umifQ)){  cerr << "Error while editing QX tags for UMI:"<<umifQ<<"#"<< endl; exit(1);   }
			    if(!al2.AddTag("RX","Z",umifS)){ cerr << "Error while editing RX tags for UMI:"<<umifS<<"#"<< endl; exit(1);   }				    
			    if(!al2.AddTag("QX","Z",umifQ)){ cerr << "Error while editing QX tags for UMI:"<<umifQ<<"#"<< endl; exit(1);   }
			}			
		    }

		    if(umir != 0){
			//al is the reverse read
			umirS          = reverseComplement(al.QueryBases.substr(0,umir));
			umirQ          =                   al.Qualities.substr( 0,umir);
			reverse(umirQ.begin(),umirQ.end());
			
			al.QueryBases  = al.QueryBases.substr(umir+umtr);
			al.Qualities   = al.Qualities.substr( umir+umtr);

			
			if(umif == 0){	//no forward UMI, just add the reverse
			    if(!al.AddTag("RX","Z",umirS)){  cerr << "Error while editing RX tags for UMI:"<<umirS<<"#"<< endl;	exit(1);   }				    
			    if(!al.AddTag("QX","Z",umirQ)){  cerr << "Error while editing QX tags for UMI:"<<umirQ<<"#"<< endl;	exit(1);   }
			    if(!al2.AddTag("RX","Z",umirS)){ cerr << "Error while editing RX tags for UMI:"<<umirS<<"#"<< endl;	exit(1);   }				    
			    if(!al2.AddTag("QX","Z",umirQ)){ cerr << "Error while editing QX tags for UMI:"<<umirQ<<"#"<< endl;	exit(1);   }
			}else{
			    //add both umif and umir
			    if(!al.AddTag("RX","Z",umifS+"-"+umirS)){  cerr << "Error while editing RX tags for UMI:"<<umifS<<"-"<<umirS<<"#"<< endl; exit(1);    }				    
			    if(!al.AddTag("QX","Z",umifQ+" "+umirQ)){  cerr << "Error while editing QX tags for UMI:"<<umirQ<<" "<<umirQ<<"#"<< endl; exit(1);    }
			    if(!al2.AddTag("RX","Z",umifS+"-"+umirS)){ cerr << "Error while editing RX tags for UMI:"<<umifS<<"-"<<umirS<<"#"<< endl; exit(1);    }				    
			    if(!al2.AddTag("QX","Z",umifQ+" "+umirQ)){ cerr << "Error while editing QX tags for UMI:"<<umirQ<<" "<<umirQ<<"#"<< endl; exit(1);    }
			}

		    }
		    
		    // writer.SaveAlignment(al2);		    
		    // writer.SaveAlignment(al);
		    emit(al2);
		    emit(al);

		}

		//
		//  SINGLE END
		//
	    }else{ 

		if(keepOrig){
		    orig =al;
		    al.Name = al.Name+"_M";
		}
		unsigned int orgL=al.QueryBases.length();		
		m->setRandomContext(in.firstRecord + i);
		m->processSingle(al);
		unsigned int newL=al.QueryBases.length();
				
		if(keepOrig){
		    //write duplicate
		    if(orgL  != newL){
			orig.SetIsDuplicate(true);
			emit(orig);
		    }
		}
		//if(!al.AddTag("RX","Z",umifS)){ cerr << "Error while editing tags for UMI:"<<umifS<<"#"<< endl;  exit(1); }

		
		if(umif != 0){	    
		    if(al.QueryBases.length() < (umif+umir+umtf+umtr) ){//the new fragment is shorter than both umis, cannot trim 
			al.SetIsDuplicate(true);  //leave but fail QC
		    }else{				    
			umifS         = al.QueryBases.substr(0,umif);
			umifQ         = al.Qualities.substr(0, umif);

			al.QueryBases = al.QueryBases.substr(umif+umtf);
			al.Qualities  = al.Qualities.substr( umif+umtf);

			if(umir == 0){	//no reverse UMI, add the tags now
			    if(!al.AddTag("RX","Z",umifS)){  cerr << "Error while editing RX tags for UMI:"<<umifS<<"#"<< endl; exit(1);   }				    
			    if(!al.AddTag("QX","Z",umifQ)){  cerr << "Error while editing QX tags for UMI:"<<umifQ<<"#"<< endl; exit(1);   }
			}	
		    }
		}

		if(umir != 0){	    
		    if(al.QueryBases.length() < (umif+umir+umtf+umtr) ){//the new fragment is shorter than both umis, cannot trim 
			al.SetIsDuplicate(true);  //leave but fail QC
		    }else{				    
			if(orgL  != newL){//was trimmed
			    umirS         = al.QueryBases.substr(al.QueryBases.length()-umir,  umir);
			    umirQ         = al.Qualities.substr( al.Qualities.length() -umir,  umir);
			
			    al.QueryBases = al.QueryBases.substr(0,al.QueryBases.length()-umir-umtr);
			    al.Qualities  = al.Qualities.substr( 0,al.Qualities.length() -umir-umtr);
			}else{//if wasn't trimmed we do not know the rev UMI, add Ns
			    umirS         = string(umir, 'N');
			    umirQ         = string(umir, '!');
			}
		
			if(umif == 0){	//no forward UMI, just add the reverse
			    if(!al.AddTag("RX","Z",umirS)){  cerr << "Error while editing RX tags for UMI:"<<umirS<<"#"<< endl;	exit(1);   }				    
			    if(!al.AddTag("QX","Z",umirQ)){  cerr << "Error while editing QX tags for UMI:"<<umirQ<<"#"<< endl;	exit(1);   }
			}else{
			    //add both umif and umir
			    if(!al.AddTag("RX","Z",umifS+"-"+umirS)){  cerr << "Error while editing RX tags for UMI:"<<umifS<<"-"<<umirS<<"#"<< endl; exit(1);    }				    
			    if(!al.AddTag("QX","Z",umifQ+" "+umirQ)){  cerr << "Error while editing QX tags for UMI:"<<umirQ<<" "<<umirQ<<"#"<< endl; exit(1);    }
			}
	
		    }
		}

		//writer.SaveAlignment(al);
		emit(al);



	    } //end single end
	    al2Null=true;
	}//second pair
		    

    }

    #undef emit
    return chunkToWrite;
}


//! Processes a chunk of fastq records, returns the compressed text for each output file
static FqOutChunk * processFqChunk(MergeTrimReads * m, const Settings & S, FqChunk & in){

    const bool         singleEndModeFQ = S.singleEndModeFQ;
    const unsigned int umif   = S.umif;
    const unsigned int umir   = S.umir;
    const unsigned int umtf   = S.umiTf;
    const unsigned int umtr   = S.umiTr;
    const bool         trimKey = S.trimKey;
    const bool         umi    = ((umif+umir+umtf+umtr)>0);

    string text[6];

    string umifS;
    string umirS;
    string umifQ;
    string umirQ;
    string cmt;

    //formats a processed record, and tallies what happened to it
    #define emitFq(r) {                                                                                           \
	const fqrecord & rr = (r);                                                                                \
	if(rr.paired){                                                                                            \
	    m->incrementCountall();                                                                               \
	    if(rr.code != ' '){ /*keys or chimeras*/                                                              \
		if(rr.code == 'K'){                                                                               \
		    m->incrementCountfkey();                                                                      \
		}else{                                                                                            \
		    if(rr.code == 'D'){                                                                           \
			m->incrementCountchimera();                                                               \
		    }else{                                                                                        \
			if(rr.code == 'U'){ /*UMI problem, sequence is too short*/                                \
			    m->incrementCountUMIp();                                                              \
			}else{                                                                                    \
			    cerr << "leehom: Wrong return code =\""<<rr.code<<"\""<<endl;                         \
			    exit(1);                                                                              \
			}                                                                                         \
		    }                                                                                             \
		}                                                                                                 \
		fqPut(text[S.sidx[S_PAIR1F]], "@"+rr.d1+"/1", rr.cmt, rr.s1, rr.q1);                              \
		fqPut(text[S.sidx[S_PAIR2F]], "@"+rr.d2+"/2", rr.cmt, rr.s2, rr.q2);                              \
	    }else{                                                                                                \
		if(rr.sequence != ""){ /*new sequence*/                                                           \
		    fqPut(text[S.sidx[S_SINGLE]], "@"+rr.d1, rr.cmt, rr.sequence, rr.quality);                    \
		    if( (rr.sequence.length()+umif+umir+umtf+umtr) > max(rr.s1.length(),rr.s2.length()) ){        \
			m->incrementCountmergedoverlap();                                                         \
		    }else{                                                                                        \
			m->incrementCountmerged();                                                                \
		    }                                                                                             \
		}else{ /*keep as is*/                                                                             \
		    m->incrementCountnothing();                                                                   \
		    if(trimKey){                                                                                  \
			fqPut(text[S.sidx[S_PAIR1]], "@"+rr.d1+"/1", rr.cmt, rr.s1.substr(S.key1.length()), rr.q1.substr(S.key1.length())); \
			fqPut(text[S.sidx[S_PAIR2]], "@"+rr.d2+"/2", rr.cmt, rr.s2.substr(S.key2.length()), rr.q2.substr(S.key2.length())); \
		    }else{                                                                                        \
			fqPut(text[S.sidx[S_PAIR1]], "@"+rr.d1+"/1", rr.cmt, rr.s1, rr.q1);                       \
			fqPut(text[S.sidx[S_PAIR2]], "@"+rr.d2+"/2", rr.cmt, rr.s2, rr.q2);                       \
		    }                                                                                             \
		}                                                                                                 \
	    }                                                                                                     \
	}else{ /*single-end*/                                                                                     \
	    if(rr.code != ' '){ /*either chimera or missing key*/                                                 \
		if(rr.code == 'K'){                                                                               \
		    m->incrementCountfkey();                                                                      \
		}else{                                                                                            \
		    if(rr.code == 'D'){                                                                           \
			m->incrementCountchimera();                                                               \
		    }else{                                                                                        \
			if(rr.code == 'U'){ /*UMI problem, sequence is too short*/                                \
			    m->incrementCountUMIp();                                                              \
			}else{                                                                                    \
			    cerr << "leehom: Wrong return code =\""<<rr.code<<"\""<<endl;                         \
			    exit(1);                                                                              \
			}                                                                                         \
		    }                                                                                             \
		}                                                                                                 \
		fqPut(text[S.sidx[S_SINGLEF]], rr.d1, rr.cmt, rr.s1, rr.q1);                                      \
	    }else{                                                                                                \
		if(rr.sequence != ""){ /*new sequence*/                                                           \
		    m->incrementCounttrimmed();                                                                   \
		    fqPut(text[S.sidx[S_SINGLE]], rr.d1, rr.cmt, rr.sequence, rr.quality);                        \
		}else{                                                                                            \
		    m->incrementCountnothing();                                                                   \
		    if(trimKey){                                                                                  \
			fqPut(text[S.sidx[S_SINGLE]], rr.d1, rr.cmt, rr.s1.substr(S.key1.length()), rr.q1.substr(S.key1.length())); \
		    }else{                                                                                        \
			fqPut(text[S.sidx[S_SINGLE]], rr.d1, rr.cmt, rr.s1, rr.q1);                               \
		    }                                                                                             \
		}                                                                                                 \
	    }                                                                                                     \
	}                                                                                                         \
    }

    for(unsigned i=0;i<in.recs.size();i++){
	m->setRandomContext(in.firstRecord + i);
	//cerr<<"process "<<in.recs[i].d1<<endl;
	fqrecord toadd;


	
	if(!singleEndModeFQ){
	    

	    merged result=	m->process_PE(in.recs[i].s1,in.recs[i].q1,
						in.recs[i].s2,in.recs[i].q2);
	    
	    

	    if(umi){
		cmt.clear();
	        if( result.sequence != ""){//has new sequence

		    if( result.sequence.length() < (umif+umir+umtf+umtr) ){//the resulting read is too short
			if( result.code == ' '){
			    result.code = 'U';
			}

			if(umif != 0 && umir != 0){
			    cmt = cmt+"\tRX:Z:"+string(umif, 'N')+"-"+string(umir, 'N');
			    cmt = cmt+"\tQX:Z:"+string(umif, '!')+" "+string(umir, '!');
			}

			if(umif != 0 && umir == 0){
			    cmt = cmt+"\tRX:Z:"+string(umif, 'N');
			    cmt = cmt+"\tQX:Z:"+string(umif, '!');
			}

			if(umif == 0 && umir != 0){
			    cmt = cmt+"\tRX:Z:"+string(umir, 'N');
			    cmt = cmt+"\tQX:Z:"+string(umir, '!');
			}
			
		    }else{//resulting read is at least the size of both UMIs

			//if UMI forward
			if(umif != 0){	    
			    umifS           = result.sequence.substr(0,umif);
			    umifQ           = result.quality.substr(0, umif);

			    result.sequence = result.sequence.substr(umif+umtf);
			    result.quality  = result.quality.substr( umif+umtf);
			    
			    if(umir == 0){//add them now
				cmt = cmt+"\tRX:Z:"+umifS;
				cmt = cmt+"\tQX:Z:"+umifQ;
			    }			    
			}
			
			//if UMI reverse
			if(umir != 0){
			    umirS           = result.sequence.substr(  result.sequence.length()-umir,  umir);
			    umirQ           = result.quality.substr(   result.quality.length() -umir,  umir);
			    
			    result.sequence = result.sequence.substr(0,result.sequence.length()-umir-umtr);
			    result.quality  = result.quality.substr(0, result.quality.length() -umir-umtr);
			    
			    if(umif != 0){
				cmt = cmt+"\tRX:Z:"+umifS+"-"+umirS;
				cmt = cmt+"\tQX:Z:"+umifQ+" "+umirQ;
			    }else{//no fwd UMI				
				cmt = cmt+"\tRX:Z:"+umirS;
				cmt = cmt+"\tQX:Z:"+umirQ;
			    }
			}
		    }

		}else{
		    //kept as pairs
		    
		    if(umif != 0){
			//cut
			umifS                                 = in.recs[i].s1.substr(0,  umif);
			umifQ                                 = in.recs[i].q1.substr(0,  umif);
			
			in.recs[i].s1 = in.recs[i].s1.substr(umif+umtf);
			in.recs[i].q1 = in.recs[i].q1.substr(umif+umtf);	 	

			if(umir == 0){//add them now
			    cmt = cmt+"\tRX:Z:"+umifS;
			    cmt = cmt+"\tQX:Z:"+umifQ;
			}	
		    }
		    

		    if(umir != 0){
			umirS                                 = reverseComplement(in.recs[i].s2.substr(0,  umir));
			umirQ                                 =                   in.recs[i].q2.substr(0,  umir);
			reverse(umirQ.begin(),umirQ.end());

			in.recs[i].s2 = in.recs[i].s2.substr(0,in.recs[i].s2.length()-umir-umtr);
			in.recs[i].q2 = in.recs[i].q2.substr(0,in.recs[i].q2.length()-umir-umtr);
			
			if(umif != 0){
			    cmt = cmt+"\tRX:Z:"+umifS+"-"+umirS;
			    cmt = cmt+"\tQX:Z:"+umifQ+" "+umirQ;
			}else{//no fwd UMI	
			    cmt = cmt+"\tRX:Z:"+umirS;
			    cmt = cmt+"\tQX:Z:"+umirQ;
			}
		    }

		}

	    }//end UMI
	    
	    //m->incrementCountall();
	    toadd.code         = result.code;
	    
	    toadd.sequence     = result.sequence;
	    toadd.quality      = result.quality;
	    toadd.cmt          = cmt;

	    toadd.d1           = in.recs[i].d1;
	    toadd.s1           = in.recs[i].s1;
	    toadd.q1           = in.recs[i].q1;
	    
	    toadd.d2           = in.recs[i].d2;
	    toadd.s2           = in.recs[i].s2;
	    toadd.q2           = in.recs[i].q2;

	    toadd.paired       = true;

	    // continue;



	}else{//V single mode below V
		
		
	    merged result=m->process_SR(in.recs[i].s1,
					  in.recs[i].q1); //*(fo1->getSeq()),*(fo1->getQual()));


	    if(umi){
		cmt.clear();
	        if( result.sequence != ""){//has new sequence

		    if( result.sequence.length() < (umif+umir+umtf+umtr) ){//the resulting read is too short
			if( result.code == ' '){
			    result.code = 'U';
			}

			if(umif != 0 && umir != 0){
			    cmt = cmt+"\tRX:Z:"+string(umif, 'N')+"-"+string(umir, 'N');
			    cmt = cmt+"\tQX:Z:"+string(umif, '!')+" "+string(umir, '!');
			}

			if(umif != 0 && umir == 0){
			    cmt = cmt+"\tRX:Z:"+string(umif, 'N');
			    cmt = cmt+"\tQX:Z:"+string(umif, '!');
			}

			if(umif == 0 && umir != 0){
			    cmt = cmt+"\tRX:Z:"+string(umir, 'N');
			    cmt = cmt+"\tQX:Z:"+string(umir, '!');
			}
			
		    }else{//resulting read is at least the size of both UMIs

			//if UMI forward
			if(umif != 0){	    
			    umifS           = result.sequence.substr(0,umif);
			    umifQ           = result.quality.substr(0, umif);

			    result.sequence = result.sequence.substr(umif+umtf);
			    result.quality  = result.quality.substr( umif+umtf);
			    
			    if(umir == 0){//add them now
				cmt = cmt+"\tRX:Z:"+umifS;
				cmt = cmt+"\tQX:Z:"+umifQ;
			    }			    
			}
			
			//if UMI reverse
			if(umir != 0){
			    umirS           = result.sequence.substr(  result.sequence.length()-umir,  umir);
			    umirQ           = result.quality.substr(   result.quality.length() -umir,  umir);
			    
			    result.sequence = result.sequence.substr(0,result.sequence.length()-umir-umtr);
			    result.quality  = result.quality.substr(0, result.quality.length() -umir-umtr);
			    
			    if(umif != 0){
				cmt = cmt+"\tRX:Z:"+umifS+"-"+umirS;
				cmt = cmt+"\tQX:Z:"+umifQ+" "+umirQ;
			    }else{//no fwd UMI				
				cmt = cmt+"\tRX:Z:"+umirS;
				cmt = cmt+"\tQX:Z:"+umirQ;
			    }
			}
		    }

		}else{
		    //kept as single-end read
		    
		    if(umif != 0){
			//cut
			umifS                                 = in.recs[i].s1.substr(0,  umif);
			umifQ                                 = in.recs[i].q1.substr(0,  umif);
			
			in.recs[i].s1 = in.recs[i].s1.substr(umif+umtf);
			in.recs[i].q1 = in.recs[i].q1.substr(umif+umtf);	 	

			if(umir == 0){//add them now
			    cmt = cmt+"\tRX:Z:"+umifS;
			    cmt = cmt+"\tQX:Z:"+umifQ;
			}	
		    }
		    

		    if(umir != 0){//we cannot infer the UMI for the rev
			umirS                                 = string(umir, 'N');
			umirQ                                 = string(umir, '!');
						
			if(umif != 0){
			    cmt = cmt+"\tRX:Z:"+umifS+"-"+umirS;
			    cmt = cmt+"\tQX:Z:"+umifQ+" "+umirQ;
			}else{//no fwd UMI	
			    cmt = cmt+"\tRX:Z:"+umirS;
			    cmt = cmt+"\tQX:Z:"+umirQ;
			}
		    }

		}

	    }//end UMI

	   

	    m->incrementCountall();
	    toadd.code = result.code;
	    toadd.paired  =  false;


	    
	    toadd.sequence     = result.sequence;
	    toadd.quality      = result.quality;
	    toadd.cmt          = cmt;
	    toadd.d1           = in.recs[i].d1;
	    toadd.s1           = in.recs[i].s1;
	    toadd.q1           = in.recs[i].q1;
	    


	}
	


	//totalSeqs++;
    
	emitFq(toadd);


    }

    #undef emitFq

    FqOutChunk * out = new FqOutChunk();
    for(int k=0;k<6;k++)
	if(!text[k].empty())
	    gzipMember(text[k], out->data[k], S.gzLevel);
    return out;
}


//! Worker: processes the chunks until there are none left
template<class In, class Out, class Fn>
static void workerLoop(OrderedPipeline<In,Out> * pipe,
		       Fn process,
		       vector<MergeTrimReads *> * instances,
		       mutex * mtx){
    MergeTrimReads * m = NULL;
    uint64_t rank;
    In *     in;
    struct StatsAtExit { ~StatsAtExit(){ if(statsEnabled()) statsAdd(statsWorkers, threadCpuSeconds()); } } statsAtExit;
    while(pipe->take(rank, in)){
	if(m == NULL){
	    //created now and not when the thread starts, the adapters might have been inferred from the first chunk
	    m = new MergeTrimReads(*mtr);
	    lock_guard<mutex> lock(*mtx);
	    instances->push_back(m);
	}
	Out * out = process(m, *in);
	delete in;
	pipe->put(rank, out);
    }
}


//! Adds the tally of the workers to the main object and prints the report
static void reportTallies(const Settings & S, vector<MergeTrimReads *> & instances){
    for(size_t i=0;i<instances.size();i++){
	mtr->addCounts(*instances[i]);
	delete instances[i];
    }
    instances.clear();

    cerr <<mtr->reportSingleLine()<<endl;
    statsReport();

    if(S.printLog){
	ofstream fileLog;
	fileLog.open(S.logFileName.c_str());

	if (fileLog.is_open()){
	    fileLog <<mtr->reportMultipleLines() <<endl;
	}else{
	    cerr << "Unable to print to file "<<S.logFileName<<endl;
	}
	fileLog.close();
    }
}


/////////////////////////////////////////////////////////////////////
//                              BAM                                //
/////////////////////////////////////////////////////////////////////

typedef OrderedPipeline<BamChunk,BamChunk> BamPipeline;

static void bamWriterLoop(BamPipeline * pipe, samFile * fpOut, sam_hdr_t * hdrOut){
    struct StatsAtExit { ~StatsAtExit(){ if(statsEnabled()) statsAdd(statsWriter, threadCpuSeconds()); } } statsAtExit;
    BamChunk * c;
    while(pipe->nextOutput(c)){
	for(size_t i=0;i<c->recs.size();i++){
	    if(sam_write1(fpOut, hdrOut, c->recs[i]) < 0){
		cerr << "Error while writing the output BAM file"<<endl;
		exit(1);
	    }
	}
	delete c;
	pipe->written();
    }
}


//! Reads the next record (a pair, or a single read) from the fastq files. Returns 0, or 1 if the files do not match.
static int nextFqRecord(FastQParser * fqp1, FastQParser * fqp2, bool singleEndModeFQ, fqrecord & fqtoadd){
	    FastQObj * fo1=fqp1->getData();
	    vector<string> def1=allTokens( *(fo1->getID()), ' '  );
	    string def1s=def1[0];
	

	    FastQObj * fo2;
	    string def2s;
	    string ext2s;
	    
	    if(!singleEndModeFQ){
		if(!fqp2->hasData()){
		    cerr << "ERROR: Discrepency between fastq files at record " <<  *(fo1->getID()) <<endl;
		    return 1;
		}

		fo2=fqp2->getData();
		vector<string> def2=allTokens( *(fo2->getID()), ' ' );
		def2s=def2[0];


		if(strEndsWith(def1s,"/1")){
		    def1s=def1s.substr(0,def1s.size()-2);
		}
		if(strEndsWith(def2s,"/2")){
		    def2s=def2s.substr(0,def2s.size()-2);
		}

		if(strBeginsWith(def1s,"@")){
		    def1s=def1s.substr(1,def1s.size()-1);
		}
		if(strBeginsWith(def2s,"@")){
		    def2s=def2s.substr(1,def2s.size()-1);
		}


		if(def1s != def2s){
		    cerr << "ERROR: Discrepency between fastq files, different names " << *(fo1->getID()) <<" and "<< *(fo2->getID()) <<endl;
		    return 1;
		}

		//PAIRED COMPUTATION
		fqtoadd.paired = true;

		fqtoadd.d1     = def1s;//*(fo1->getID()  );
		fqtoadd.s1     = *(fo1->getSeq() );
		fqtoadd.q1     = *(fo1->getQual());

		fqtoadd.d2     = def2s;//*(fo2->getID()  );
		fqtoadd.s2     = *(fo2->getSeq() );
		fqtoadd.q2     = *(fo2->getQual());
		
	    }else{		
		//SINGLE COMPUTATION
		fqtoadd.paired = false;

		fqtoadd.d1     = def1s;//*(fo1->getID());
		fqtoadd.s1     = *(fo1->getSeq());
		fqtoadd.q1     = *(fo1->getQual());

	    }

    return 0;
}

//! Fills b with an unmapped read taken from a fastq record, the qualities are converted from ASCII to phred
static void fastqToBamRecord(bam1_t * b, const string & name, const string & seq, const string & qual, uint16_t flag, int offset){
    if(seq.size() != qual.size()){
	cerr << "ERROR: the sequence and the qualities of read "<<name<<" do not have the same length"<<endl;
	exit(1);
    }
    const string qname = strBeginsWith(name,"@") ? name.substr(1) : name;
    vector<uint8_t> q(qual.size());
    for(size_t i=0;i<qual.size();i++){
	int v = int((unsigned char) qual[i]) - offset;
	q[i]  = uint8_t( v < 0 ? 0 : (v > 93 ? 93 : v) );
    }
    if(bam_set1(b, qname.size(), qname.c_str(), flag, -1, -1, 0, 0, NULL, -1, -1, 0,
		seq.size(), seq.c_str(), (const char *) (q.empty() ? NULL : &q[0]), 0) < 0){
	cerr << "ERROR: cannot build the BAM record of read "<<name<<endl;
	exit(1);
    }
}

static int runBam(Settings & S, int argc, char * argv[]){

    if( S.bamFileOUT == ""  ){
	S.bamFileOUT = "/dev/stdout" ;
    }

    samFile *     fpIn  = NULL;
    sam_hdr_t *   hdrIn = NULL;
    FastQParser * fqp1  = NULL;
    FastQParser * fqp2  = NULL;
    bam1_t *      pendingMate = NULL;     //fastq input: second read of the pair just read

    if(S.fastqToBam){
	if(S.fastqfile1 == ""){
	    cerr<<"ERROR : Must specify as least the first file for fastq"<<endl;
	    return 1;
	}
	if(S.singleEndModeFQ && S.detectAdapt){
	    cerr<<"The automatic detection of the adapters is not yet coded for single-end reads"<<endl;
	    return 1;
	}
	fqp1 = new FastQParser (S.fastqfile1);
	if(!S.singleEndModeFQ)
	    fqp2 = new FastQParser (S.fastqfile2);
	const string hdrText = "@HD\tVN:1.6\tSO:unsorted";
	hdrIn = sam_hdr_parse(hdrText.size(), hdrText.c_str());
    }else{
	fpIn = sam_open(S.bamFile.c_str(), "r");
	if(fpIn == NULL){
	    cerr << "Could not open input BAM file" << S.bamFile<<endl;
	    return 1;
	}
	hdrIn = sam_hdr_read(fpIn);
    }
    if(hdrIn == NULL){
	cerr << "Could not read the header of the input BAM file" << S.bamFile<<endl;
	return 1;
    }

    //adding this program to the header
    string pID          = "leeHom";
    string pName        = "leeHom";
    string pCommandLine = "";
    for(int i=0;i<(argc);i++){
	pCommandLine += (string(argv[i])+" ");
    }
    string headerText = (sam_hdr_length(hdrIn) > 0) ? string(sam_hdr_str(hdrIn), sam_hdr_length(hdrIn)) : string("");
    while(!headerText.empty() && headerText[headerText.size()-1] == '\n')   //putProgramInHeaderHTS adds its own line breaks
	headerText.erase(headerText.size()-1);
    putProgramInHeaderHTS(&headerText,pID,pName,pCommandLine,returnGitHubVersion(string(argv[0]),".."));
    headerText += "\n";

    sam_hdr_t * hdrOut = sam_hdr_parse(headerText.size(), headerText.c_str());
    if(hdrOut == NULL || sam_hdr_nref(hdrOut) != sam_hdr_nref(hdrIn)){
	cerr << "Could not build the header of the output BAM file" <<endl;
	return 1;
    }

    //we will not call bgzip with full compression, good for piping into another program to
    //lessen the load on the CPU
    string bamMode = "wb";
    if(S.produceUnCompressedBAM)
	bamMode = "wb0";
    else if(S.gzLevel >= 0)
	bamMode += char('0'+S.gzLevel);
    samFile * fpOut = sam_open(S.bamFileOUT.c_str(), bamMode.c_str());
    if(fpOut == NULL){
	cerr << "Could not open output BAM file "<<S.bamFileOUT << endl;
	return 1;
    }
    if(sam_hdr_write(fpOut, hdrOut) < 0){
	cerr << "Could not write the header of the output BAM file "<<S.bamFileOUT << endl;
	return 1;
    }

    if(numberOfThreads > 1){
	if(fpIn != NULL)
	    hts_set_threads(fpIn,  numberOfThreads);
	hts_set_threads(fpOut, numberOfThreads);
    }

    //next record of the input: 1 got one, 0 end of the input
    auto readRecord = [&](bam1_t * b) -> int {
	if(fpIn != NULL)
	    return (sam_read1(fpIn, hdrIn, b) >= 0) ? 1 : 0;
	if(pendingMate != NULL){
	    bam_copy1(b, pendingMate);
	    bam_destroy1(pendingMate);
	    pendingMate = NULL;
	    return 1;
	}
	if(!fqp1->hasData())
	    return 0;
	fqrecord r;
	if(nextFqRecord(fqp1, fqp2, S.singleEndModeFQ, r) != 0)
	    exit(1);
	if(r.paired){
	    fastqToBamRecord(b, r.d1, r.s1, r.q1, BAM_FPAIRED|BAM_FUNMAP|BAM_FMUNMAP|BAM_FREAD1, S.fqQualOffset);
	    pendingMate = bam_init1();
	    fastqToBamRecord(pendingMate, r.d2, r.s2, r.q2, BAM_FPAIRED|BAM_FUNMAP|BAM_FMUNMAP|BAM_FREAD2, S.fqQualOffset);
	}else{
	    fastqToBamRecord(b, r.d1, r.s1, r.q1, BAM_FUNMAP, S.fqQualOffset);
	}
	return 1;
    };

    BamPipeline pipe(3*numberOfThreads+3);
    vector<MergeTrimReads *> instances;
    mutex                    instancesMutex;

    vector<thread> workers;
    for(int i=0;i<numberOfThreads;i++){
	workers.push_back( thread( workerLoop<BamChunk,BamChunk,function<BamChunk*(MergeTrimReads*,const BamChunk&)> >,
				   &pipe,
				   function<BamChunk*(MergeTrimReads*,const BamChunk&)>( [&S](MergeTrimReads * m, const BamChunk & c){ return processBamChunk(m, S, c); } ),
				   &instances,
				   &instancesMutex) );
    }
    thread writerThread(bamWriterLoop, &pipe, fpOut, hdrOut);

    unsigned int rank    = 0;
    unsigned int counter = 0;
    uint64_t     nread   = 0;

    BamChunk * currentChunk = new BamChunk();
    bam1_t *   b            = bam_init1();

    while( readRecord(b) ){

	if( !(b->core.flag & BAM_FUNMAP) || bam_aux_get(b,"NM") != NULL || bam_aux_get(b,"MD") != NULL ){
	    if(!S.allowAligned){
		cerr << "Reads should not be aligned" << endl;
		exit(1);
	    }else{
		//should we remove tags ?
	    }
	}

	currentChunk->recs.push_back(b);
	b = bam_init1();
	nread++;

	if( counter== ( (S.detectAdapt&&(rank==0))?(sizeChunkDefault*10-1):((sizeChunk-1))) ){

	    if(S.detectAdapt){
		bool adpd=detectAdaptBAM( currentChunk );
		if(adpd){//worked
		    S.detectAdapt=false;
		}else{
		    cerr << "Failed to auto detect adapters, please contact the vendor" << endl;
		    exit(1);
		}
	    }

	    pipe.submit(currentChunk);
	    rank++;

	    currentChunk = new BamChunk();
	    currentChunk->firstRecord = nread;
	    counter=0;
	}else{
	    counter++;
	}
    }
    bam_destroy1(b);
    if(fpIn != NULL)
	sam_close(fpIn);
    delete fqp1;
    delete fqp2;

    //last chunk
    pipe.submit(currentChunk);
    pipe.closeInput();

    for(size_t i=0;i<workers.size();i++)
	workers[i].join();
    writerThread.join();
    if(statsEnabled()) statsAdd(statsReader, threadCpuSeconds());

    sam_close(fpOut);
    sam_hdr_destroy(hdrIn);
    sam_hdr_destroy(hdrOut);

    reportTallies(S, instances);
    return 0;
}



/////////////////////////////////////////////////////////////////////
//                             FASTQ                               //
/////////////////////////////////////////////////////////////////////

typedef OrderedPipeline<FqChunk,FqOutChunk> FqPipeline;

static void fqWriterLoop(FqPipeline * pipe, FILE * const * sinks, bool * wrote){
    struct StatsAtExit { ~StatsAtExit(){ if(statsEnabled()) statsAdd(statsWriter, threadCpuSeconds()); } } statsAtExit;
    FqOutChunk * c;
    while(pipe->nextOutput(c)){
	for(int k=0;k<6;k++){
	    if(!c->data[k].empty()){
		if(fwrite(c->data[k].data(), 1, c->data[k].size(), sinks[k]) != c->data[k].size()){
		    cerr << "Error while writing the output fastq file"<<endl;
		    exit(1);
		}
		wrote[k]=true;
	    }
	}
	delete c;
	pipe->written();
    }
}


static int runFastq(Settings & S){

    const bool singleEndModeFQ = S.singleEndModeFQ;

    if( S.bamFileOUT != ""  || S.produceUnCompressedBAM || S.allowAligned){
	cerr<<"ERROR : Cannot specify options like -o, -u or --allowAligned for fastq"<<endl;
	return 1;
    }

    if(S.fastqfile1 == ""){
	cerr<<"ERROR : Must specify as least the first file for fastq"<<endl;
	return 1;
    }

    if(singleEndModeFQ && S.detectAdapt){
	cerr<<"The automatic detection of the adapters is not yet coded for single-end reads"<<endl;
	return 1;
    }

    //names of the output files for each stream
    string names[6];
    const bool interleaved = (S.fastqoutfile.empty() && !S.fastqoutfileb);

    if(singleEndModeFQ){
	if(!S.fastqoutfile.empty()){
	    names[S_SINGLE]  = S.fastqoutfile+".fq.gz";
	    names[S_SINGLEF] = S.fastqoutfile+".fail.fq.gz";
	}else{
	    if(S.fastqoutfileb){
		names[S_SINGLE]  = S.fastqoutfilese;
		names[S_SINGLEF] = S.fastqoutfilesef;
	    }else{
		names[S_SINGLE]  = S.intfastqoutfile;
	    }
	}
    }else{
	if(!S.fastqoutfile.empty()){
	    names[S_SINGLE]  = S.fastqoutfile+".fq.gz";
	    names[S_PAIR1]   = S.fastqoutfile+"_r1.fq.gz";
	    names[S_PAIR2]   = S.fastqoutfile+"_r2.fq.gz";
	    names[S_SINGLEF] = S.fastqoutfile+".fail.fq.gz";
	    names[S_PAIR1F]  = S.fastqoutfile+"_r1.fail.fq.gz";
	    names[S_PAIR2F]  = S.fastqoutfile+"_r2.fail.fq.gz";
	}else{
	    if(S.fastqoutfileb){
		names[S_SINGLE]  = S.fastqoutfilese;
		names[S_PAIR1]   = S.fastqoutfilepe1;
		names[S_PAIR2]   = S.fastqoutfilepe2;
		names[S_SINGLEF] = S.fastqoutfilesef;
		names[S_PAIR1F]  = S.fastqoutfilepe1f;
		names[S_PAIR2F]  = S.fastqoutfilepe2f;
	    }else{
		names[S_SINGLE]  = S.intfastqoutfile;
	    }
	}
    }

    //with a single output file all the streams go to it, in the order of the reads
    for(int k=0;k<6;k++)
	S.sidx[k] = interleaved ? 0 : k;

    FILE * sinks[6]  = {NULL,NULL,NULL,NULL,NULL,NULL};
    bool   wrote[6]  = {false,false,false,false,false,false};
    for(int k=0;k<6;k++){
	if(!names[k].empty()){
	    sinks[k] = fopen(names[k].c_str(), "wb");
	    if(sinks[k] == NULL){
		cerr<<"Cannot write to file "<<names[k]<<endl;
		return 1;
	    }
	}
    }

    FastQParser * fqp1;
    FastQParser * fqp2 = NULL;

    fqp1 = new FastQParser (S.fastqfile1);
    if(!singleEndModeFQ)
	fqp2 = new FastQParser (S.fastqfile2);

    FqPipeline pipe(3*numberOfThreads+3);
    vector<MergeTrimReads *> instances;
    mutex                    instancesMutex;

    vector<thread> workers;
    for(int i=0;i<numberOfThreads;i++){
	workers.push_back( thread( workerLoop<FqChunk,FqOutChunk,function<FqOutChunk*(MergeTrimReads*,FqChunk&)> >,
				   &pipe,
				   function<FqOutChunk*(MergeTrimReads*,FqChunk&)>( [&S](MergeTrimReads * m, FqChunk & c){ return processFqChunk(m, S, c); } ),
				   &instances,
				   &instancesMutex) );
    }
    thread writerThread(fqWriterLoop, &pipe, (FILE * const *) sinks, wrote);

    unsigned int rank    = 0;
    unsigned int counter = 0;
    uint64_t     nread   = 0;
    bool         detectAdapt = S.detectAdapt;

    FqChunk * currentChunkfq = new FqChunk();

	while(fqp1->hasData()){

	    fqrecord fqtoadd;
	    if(nextFqRecord(fqp1, fqp2, singleEndModeFQ, fqtoadd) != 0)
		return 1;

	nread++;

	if(counter== ( (detectAdapt && (rank==0))?(sizeChunkDefault*10-1):((sizeChunk-1))) ){

	    currentChunkfq->recs.push_back(fqtoadd);
	    if(detectAdapt){
		bool adpd=detectAdaptFQ( currentChunkfq );
		if(adpd){//worked
		    detectAdapt=false;
		}else{
		    cerr << "Failed to auto detect adapters, please contact the vendor" << endl;
		    return 1;
		}
	    }

	    pipe.submit(currentChunkfq);
	    rank++;

	    currentChunkfq = new FqChunk();
	    currentChunkfq->firstRecord = nread;
	    counter=0;
	}else{
	    currentChunkfq->recs.push_back(fqtoadd);
	    counter++;
	}

    }//for each data record

    delete fqp1;
    if(!singleEndModeFQ){
	delete fqp2;
    }

    //last chunk
    pipe.submit(currentChunkfq);
    pipe.closeInput();

    for(size_t i=0;i<workers.size();i++)
	workers[i].join();
    writerThread.join();

    if(statsEnabled()) statsAdd(statsReader, threadCpuSeconds());

    for(int k=0;k<6;k++){
	if(sinks[k] != NULL){
	    if(!wrote[k]){ //a gzip file with no data still has to contain one (empty) member
		string empty;
		gzipMember("", empty, S.gzLevel);
		fwrite(empty.data(), 1, empty.size(), sinks[k]);
	    }
	    if(fclose(sinks[k]) != 0){
		cerr<<"Error while closing the file "<<names[k]<<endl;
		return 1;
	    }
	}
    }

    reportTallies(S, instances);
    return 0;
}



int main (int argc, char *argv[]) {

    if(sizeChunk%2 != 0){
	cerr<<"The size of the chunks has to be an even number"<<endl;
	return 1;
    }

    bool produceUnCompressedBAM=false;
    bool verbose=false;
    bool ancientDNA=false;
    bool keepOrig=false;

    string adapter_F      =options_adapter_F_BAM;
    string adapter_S      =options_adapter_S_BAM;
    string adapter_chimera=options_adapter_chimera_BAM;
    string key="";
    string ikey="";
    bool   trimKey=false;
    
    bool   allowMissing=false;
    int    trimCutoff=1;

    bool allowAligned=false;
    bool printLog=false;
    string logFileName;


    string bamFile;
    string bamFileOUT="";

    string key1;
    string key2;
    
    bool useDist=false;
    double location=-1.0;
    double scale   =-1.0;
    bool detectAdapt=false;

    bool specifiedAdapt=false;
    bool fastqFormat=false;
    
    
    string fastqfile1   = "";
    string fastqfile2   = "";
    string fastqoutfile = "";

    bool   fastqoutfileb=false;
    string fastqoutfilese   = "";
    string fastqoutfilesef  = "";
    string fastqoutfilepe1  = "";
    string fastqoutfilepe1f = "";
    string fastqoutfilepe2  = "";
    string fastqoutfilepe2f = "";

    string intfastqoutfile = "";
    
    bool singleEndModeFQ=true;
    bool fastqToBam=false;      //fastq input, BAM output
    int    qualOffset     = 33;

    bool seedSpecified=false;
    unsigned int seed=0;
    
    int          gzLevel=-1;       //-1: default compression level
    int          quickBases=0;     //0: quick mode off
    double       quickFraction=0.4;

    unsigned int umif=0;
    unsigned int umir=0;
    unsigned int umiTf=0;//trim after
    unsigned int umiTr=0;
    

    const string usage1=string(string(argv[0])+
			      
			      " [options] BAMfile"+"\n"+
			      "\nThis program takes an unaligned BAM where mates are consecutive\nor fastq files and trims and merges reads\n"+

			      "\n\tYou can specify a unaligned bam file or one or two fastq :\n"+			      
			      "\t\t"+"-fq1 [file]" +"\t\t"+"First fastq file"+"\n"+
			      "\t\t"+"-fq2 [file]" +"\t\t"+"Second  fastq file (for paired-end)"+"\n"+
			      "\t\t"+"-fqo [string]" +"\t\t"+"Output fastq prefix to file (type --help for more options)"+"\n"+
			      "\t\t"+"                \t\t"+"With fastq input, -o or -u write an unaligned BAM (mates consecutive) instead"+"\n");

    const string usage2=string(string("or you can use specify manually those:\n")+
			       "\t\t"+"-fqose   [file]" +"\t\t"+"Output file for fastq single-end (or merged)"+"\n"+
			       "\t\t"+"-fqosef  [file]" +"\t\t"+"Output file for fastq single-end (or merged) that failed"+"\n"+
			       "\t\t"+"-fqope1  [file]" +"\t\t"+"Output file for fastq first paired-end (or merged)"+"\n"+
			       "\t\t"+"-fqope1f [file]" +"\t\t"+"Output file for fastq first paired-end (or merged) that failed"+"\n"+
			       "\t\t"+"-fqope2  [file]" +"\t\t"+"Output file for fastq first paired-end (or merged)"+"\n"+
			       "\t\t"+"-fqope2f [file]" +"\t\t"+"Output file for fastq first paired-end (or merged) that failed"+"\n");
			       

    const string usage3=string(string("")+
			       "\t\t"+"-fqt [file]" +"\t\t"+"Output interweaved fastq"+"\n\n"+

			      //"\t"+"-p , --PIPE"+"\n\t\t"+"Read BAM from and write it to PIPE"+"\n"+
			      "\t"+"-o , --outfile" +"\t\t"+"Output (BAM format)."+"\n"+


			      "\t"+"-u            " +"\t\t"+"Produce uncompressed bam (good for pipe)"+"\n"+

			      //	"\t"+" , --outprefix" +"\n\t\t"+"Prefix for output files (default '"+outprefix+"')."+"\n"+
			      //"\t"+" , --SAM" +"\n\t\t"+"Output SAM not BAM."+"\n"+
			      "\t"+"--aligned" +"\t\t"+"Allow reads to be aligned (default "+boolStringify(allowAligned)+")"+"\n"+
			      "\t"+"-v , --verbose" +"\t\t"+"Turn all messages on (default "+boolStringify(verbose)+")"+"\n"+
			      "\t"+"--version" +"\t\t"+"Just print version and quit"+"\n"+
			      "\t"+"--log [log file]" +"\t"+"Print a tally of merged reads to this log file (default only to stderr)"+"\n"+
			      "\t"+"--phred64" +"\t\t"+"Use PHRED 64 as the offset for QC scores (default : PHRED33)"+"\n"+
			      "\t"+"-t [threads]" +"\t\t"+"Use multiple cores (default : "+stringify(numberOfThreads)+")"+"\n"+
			      
			      "\n\t"+"Paired End merging/Single Read trimming  options"+"\n"+
			      "\t\t"+"You can specify either:"+"\n"+
			      "\t\t\t"+"--ancientdna"+"\t\t\t"+"ancient DNA (default "+boolStringify(ancientDNA)+")"+"\n"+
			      "\t\t"+"            "+"\t\t\t\t"+"this allows for partial overlap"+"\n"+
			      "\n\t\t"+"or if you know your size length distribution:"+"\n"+
			      "\t\t\t"+"--loc"+"\t\t\t\t"+"Location for lognormal dist. (default none)"+"\n"+
			      "\t\t\t"+"--scale"+"\t\t\t\t"+"Scale for lognormal dist. (default none)"+"\n"+
			      //			      "\t\t\t\t\t\t\tGood for merging ancient DNA reads into a single sequence\n\n"
			      "\n\t\t"+"--keepOrig"+"\t\t\t\t"+"Write original reads if they are trimmed or merged  (default "+boolStringify(keepOrig)+")"+"\n"+
			      "\t\t\t\t\t\t\tSuch reads will be marked as PCR duplicates\n\n"
			      "\n\t\t"+"--seed [seed]"+"\t\t\t\t"+"Set a seed for the random bases (default randomly generated)"+"\n"+
			      "\n"+
			      "\t\t"+"-f , --adapterFirstRead" +"\t\t\t"+"Adapter that is observed after the forward read (def. Multiplex: "+options_adapter_F_BAM.substr(0,30)+")"+"\n"+
			      "\t\t"+"-s , --adapterSecondRead" +"\t\t"+"Adapter that is observed after the reverse read (def. Multiplex: "+options_adapter_S_BAM.substr(0,30)+")"+"\n"+
			      "\t\t"+"     --auto" +"\t\t\t\t"+"Auto detect adapters, requires at least 1M reads, does NOT currently support UMIs (def. : "+boolStringify(detectAdapt)+")"+"\n"+

			      "\t\t"+"-c , --FirstReadChimeraFilter" +"\t\t"+"If the forward read looks like this sequence, the cluster is filtered out.\n\t\t\t\t\t\t\tProvide several sequences separated by comma (def. Multiplex: "+options_adapter_chimera_BAM.substr(0,30)+")"+"\n"+
			      "\t\t"+"--ikey"+"\t\t\t\t"+"Key sequence at the beginning of the insert. (default '"+ikey+"')"+"\n"+
			      "\t\t"+"-k , --key"+"\t\t\t\t"+"Key sequence with which each sequence starts. Comma separate for forward and reverse reads. (default '"+key+"')"+"\n"+
			      "\t\t"+"--trimkey"     +"\t\t\t\t"+"Trim the key sequence even for untrimmed. (default '"+boolStringify(trimKey)+"')"+"\n"+

			      "\t\t"+"-i , --allowMissing"+"\t\t\t"+"Allow one base in one key to be missing or wrong. (default "+boolStringify(allowMissing)+")"+"\n"+

			      "\t\t"+"--umif [bp]"+"\t\t\t\t"+"Extract bp for the UMI for the forward read. (default "+stringify(umif)+")"+"\n"+
			      "\t\t"+"--umir [bp]"+"\t\t\t\t"+"Extract bp for the UMI for the reverse read. (default "+stringify(umir)+")"+"\n"+
			      
			      "\t\t"+"--umtf [bp]"+"\t\t\t\t"+"Remove overhangs of bp after the UMI for the forward read. (default "+stringify(umiTf)+")"+"\n"+
			      "\t\t"+"--umtr [bp]"+"\t\t\t\t"+"Remove overhangs of bp after the UMI for the reverse read. (default "+stringify(umiTr)+")"+"\n"+

			      "\t\t"+"--gzlevel [0-9]"+"\t\t\t\t"+"Compression level of the fastq or BAM output (default zlib's own, 6). 1 is much faster for files about 20% bigger, 0 is uncompressed (ignored with -u)"+"\n"+
			      "\t\t"+"--quickfraction [0-1]"+"\t\t\t"+"How far between the likelihood of random sequences (0) and of a perfect overlap (1) a candidate must be after --quickbases bases to be computed in full (default 0.4)"+"\n"+
				      "\t\t"+"--trimCutoff"+"\t\t\t\t"+"Lowest number of adapter bases to be observed for single Read trimming (default "+stringify(trimCutoff)+")");

    if( (argc== 1) ||
	(argc== 2 && string(argv[1]) == "-h") ){
    	cout<<"Usage:"<<endl;
    	cout<<""<<endl;
    	cout<<usage1<<usage3<<endl;
    	return 1;
    }

    if( 
    	(argc== 2 && string(argv[1]) == "-help") ||
    	(argc== 2 && string(argv[1]) == "--help") ){
    	cout<<"Usage:"<<endl;
    	cout<<""<<endl;
    	cout<<usage1<<usage2<<usage3<<endl;
    	return 1;
    }


    if( (argc== 2 && string(argv[1]) == "--version") ){
	cout<<"leeHom version: "<<VERSION<<endl;
	return 0; 
    }

    
    int lastarg=(argc-1);
    for(int i=1;i<lastarg;i++){ //all but the last arg
	//int    numberOfThreads   = 1;

	if(strcmp(argv[i],"-t") == 0 ){
	    numberOfThreads =destringify<int>(argv[i+1]);
	    //the size of the chunks does not depend on the number of threads
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"-z") == 0 ){
	    sizeChunk=destringify<int>(argv[i+1]);
	    i++;
	    continue;
	}

	// if(strcmp(argv[i],"-t") == 0 ){
	//     numberOfThreads=destringify<int>(argv[i+1]);
	//     i++;
	//     continue;
	// }


	if(strcmp(argv[i],"-fq1") == 0 ){
	    fastqfile1=string(argv[i+1]);
	    fastqFormat=true;
	    lastarg=argc;
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"-fq2") == 0 ){
	    fastqfile2=string(argv[i+1]);
	    fastqFormat=true;
	    lastarg=argc;
	    singleEndModeFQ=false;
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"-fqo") == 0 ){
	    fastqoutfile=string(argv[i+1]);
	    fastqFormat=true;
	    lastarg=argc;
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"-fqose") == 0 ){
	    fastqoutfilese=string(argv[i+1]);
	    fastqFormat=true;
	    fastqoutfileb=true;
	    lastarg=argc;
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"-fqosef") == 0 ){
	    fastqoutfilesef=string(argv[i+1]);
	    fastqFormat=true;
	    fastqoutfileb=true;
	    lastarg=argc;
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"-fqope1") == 0 ){
	    fastqoutfilepe1=string(argv[i+1]);
	    fastqFormat=true;
	    fastqoutfileb=true;
	    lastarg=argc;
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"-fqope1f") == 0 ){
	    fastqoutfilepe1f=string(argv[i+1]);
	    fastqFormat=true;
	    fastqoutfileb=true;
	    lastarg=argc;
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"-fqope2") == 0 ){
	    fastqoutfilepe2=string(argv[i+1]);
	    fastqFormat=true;
	    fastqoutfileb=true;
	    lastarg=argc;
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"-fqope2f") == 0 ){
	    fastqoutfilepe2f=string(argv[i+1]);
	    fastqFormat=true;
	    fastqoutfileb=true;
	    lastarg=argc;
	    i++;
	    continue;
	}


	
	if(strcmp(argv[i],"-fqt") == 0 ){
	    intfastqoutfile=string(argv[i+1]);
	    fastqFormat=true;
	    lastarg=argc;
	    i++;
	    continue;
	}




	if(strcmp(argv[i],"--log") == 0 ){
	    logFileName =string(argv[i+1]);
	    printLog=true;
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"-p") == 0 || strcmp(argv[i],"--PIPE") == 0 ){
	    cerr<<"This version no longer works with pipe, exiting"<<endl;
	    return 1;	    
	}

	if(strcmp(argv[i],"-u") == 0  ){
	    produceUnCompressedBAM=true;
	    continue;
	}

	if(strcmp(argv[i],"--aligned") == 0  ){
	    allowAligned=true;
	    continue;
	}

	if(strcmp(argv[i],"--auto") == 0  ){
	    detectAdapt=true;
	    continue;
	}

	if(string(argv[i]) == "--phred64"  ){
	    qualOffset=64;
	    continue;
	}

	if(string(argv[i]) == "--version"  ){
	    cout<<"leeHom version: "<<VERSION<<endl;
	    return 0; 
	}

	if(strcmp(argv[i],"-o") == 0 || strcmp(argv[i],"--outfile") == 0 ){
	    bamFileOUT =string(argv[i+1]);
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"-v") == 0 || strcmp(argv[i],"--verbose") == 0 ){
	    verbose=true;
	    continue;
	}

	if(strcmp(argv[i],"--ancientdna") == 0 ){
	    ancientDNA=true;
	    continue;
	}

	if(strcmp(argv[i],"--keepOrig") == 0 ){
	    keepOrig=true;
	    continue;
	}

	if(strcmp(argv[i],"--loc") == 0 ){
	    location =destringify<double>(argv[i+1]);
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"--scale") == 0 ){
	    scale =destringify<double>(argv[i+1]);
	    i++;
	    continue;
	}



	if(strcmp(argv[i],"-f") == 0 || strcmp(argv[i],"--adapterFirstRead") == 0 ){
	    adapter_F =string(argv[i+1]);
	    specifiedAdapt=true;
	    i++;
	    continue;
	}


	if(strcmp(argv[i],"-s") == 0 || strcmp(argv[i],"--adapterSecondRead") == 0 ){
	    adapter_S =string(argv[i+1]);
	    specifiedAdapt=true;
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"-c") == 0 || strcmp(argv[i],"--FirstReadChimeraFilter") == 0 ){
	    adapter_chimera =string(argv[i+1]);
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"--ikey") == 0 ){
	    ikey =string(argv[i+1]);
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"-k") == 0 || strcmp(argv[i],"--keys") == 0 ){
	    key =string(argv[i+1]);
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"--trimkey") ==0 ){
	    trimKey=true;
	    continue;
	}

	if(strcmp(argv[i],"-i") == 0 || strcmp(argv[i],"--allowMissing") == 0 ){
	    allowMissing=true;
	    continue;
	}

	if(strcmp(argv[i],"--seed") == 0 ){
	    seedSpecified=true;
	    seed = destringify<unsigned int>(argv[i+1]);
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"--umif") == 0 ){
	    umif = destringify<unsigned int>(argv[i+1]);
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"--umir") == 0 ){
	    umir = destringify<unsigned int>(argv[i+1]);
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"--umtf") == 0 ){
	    umiTf = destringify<unsigned int>(argv[i+1]);
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"--umtr") == 0 ){
	    umiTr = destringify<unsigned int>(argv[i+1]);
	    i++;
	    continue;
	}

	
	if(strcmp(argv[i],"--quick") == 0 ){
	    if(quickBases == 0)
		quickBases=12;
	    continue;
	}

	if(strcmp(argv[i],"--gzlevel") == 0 ){
	    gzLevel = atoi(argv[i+1]);
	    if(gzLevel < 0 || gzLevel > 9){
		cerr<<"--gzlevel must be between 0 and 9"<<endl;
		return 1;
	    }
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"--quickfraction") == 0 ){
	    quickFraction = atof(argv[i+1]);
	    if(quickFraction < 0.0 || quickFraction > 1.0){
		cerr<<"--quickfraction must be between 0 and 1"<<endl;
		return 1;
	    }
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"--quickbases") == 0 ){
	    quickBases = atoi(argv[i+1]);
	    if(quickBases < 1){
		cerr<<"--quickbases must be at least 1"<<endl;
		return 1;
	    }
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"--trimCutoff") == 0 ){
	    trimCutoff=atoi(argv[i+1]);
	    i++;
	    continue;
	}
	
	cerr<<"Unknown option "<<argv[i] <<" exiting"<<endl;
	return 1;	    
    }


    bamFile=argv[argc-1];

    fastqToBam = fastqFormat && (!bamFileOUT.empty() || produceUnCompressedBAM);

    if(!fastqFormat){
	if( (bamFile != "-") &&
	    strBeginsWith(bamFile,"-") ){
	    cerr<<"If running in bam input/output mode, the last argument should be the bam file."<<endl;
	    return 1;
	}

	if( (sizeChunk%2) != 0){
	    cerr<<"If running in bam input/output mode, the size of the chunks (option -z) must be even."<<endl;
	    return 1;	    
	}
    }else if(fastqToBam){
	if(!fastqoutfile.empty() || !intfastqoutfile.empty() || fastqoutfileb){
	    cerr<<"Cannot write both BAM (-o or -u) and fastq output."<<endl;
	    return 1;
	}
	if(allowAligned){
	    cerr<<"ERROR : Cannot specify --allowAligned for fastq"<<endl;
	    return 1;
	}
	if( (sizeChunk%2) != 0){
	    cerr<<"If writing BAM output, the size of the chunks (option -z) must be even."<<endl;
	    return 1;
	}
    }else{
	int boolCfq= (!fastqoutfile.empty()) + (!intfastqoutfile.empty()) + fastqoutfileb;

	if(boolCfq==0){
	    cerr<<"If the input is fastq, the output must be specified: -fqo or -fqt for fastq, -o or -u for BAM."<<endl;
	    return 1;	    	    
	}

	if(boolCfq>1){	
	    cerr<<"If running in fastq input/output mode, cannot specify both -fqo and -fqt (or the name of the 6 fastq files)."<<endl;
	    return 1;	    	    
	}
	// if(!fastqoutfile.empty() && !intfastqoutfile.empty()){
	//     cerr<<"If running in fastq input/output mode, cannot specify both -fqo and -fqt."<<endl;
	//     return 1;	    	    
	// }
	
	// if(fastqoutfile.empty() && intfastqoutfile.empty()){
	// }

	if(fastqoutfileb){
	    if(singleEndModeFQ){
		boolCfq= (!fastqoutfilese.empty()) + (!fastqoutfilesef.empty());
		if(boolCfq!=2){	
		    cerr<<"If running in fastq and single-end and specifying the names of the files, 2 files have to be specified, found "<<boolCfq<<endl;
		    return 1;	    	    
		}
	    }else{
		
		boolCfq= ((!fastqoutfilese.empty()) + (!fastqoutfilesef.empty()) + (!fastqoutfilepe1.empty()) + (!fastqoutfilepe1f.empty())  + (!fastqoutfilepe2.empty()) + (!fastqoutfilepe2f.empty()));

		if(boolCfq!=6){	
		    cerr<<"If running in fastq and paired-end and specifying the names of the files, 6 files have to be specified, found "<<boolCfq<<endl;
		    return 1;
		}

	    }
	}

    }


    
    if(specifiedAdapt && detectAdapt){
	cerr<<"Cannot specify adapters and automatically detect them"<<endl;
	return 1;	    
    }
    
    if( (location != -1.0 && scale == -1.0) ||
	(location == -1.0 && scale != -1.0) ){
	cerr<<"Cannot specify --location without specifying --scale"<<endl;
	return 1;	    
    }
	
    if( (location != -1.0 && scale != -1.0) ){
	useDist=true;
	    
	if(ancientDNA){
	    cerr<<"Cannot specify --location/--scale and --ancientDNA"<<endl;
	    return 1;	    
	}
    }
    
    if(key != ""){
	size_t found=key.find(",");
	if (found == string::npos){ //single end reads
	    key1=key;
	    key2="";
	} else{                     //paired-end
	    key1=key.substr(0,found);
	    key2=key.substr(found+1,key.length()-found+1);
	}
    }

    if(umiTf>0 && umif==0){
	cerr<<"Cannot specify --umtf without --umif"<<endl;
	return 1;	    
    }

    if(umiTr>0 && umir==0){
	cerr<<"Cannot specify --umtr without --umir"<<endl;
	return 1;	    
    }


    const int fqQualOffset = qualOffset;
    if(fastqToBam)
	qualOffset = 33;    //the qualities are converted to those of a BAM file

    mtr = new MergeTrimReads(adapter_F,adapter_S,adapter_chimera,
			     key1,key2,trimKey,ikey,
			     trimCutoff,allowMissing,ancientDNA,location,scale,useDist,qualOffset,seed,seedSpecified);

    mtr->setQuick(quickBases, quickFraction);

    Settings S;
    S.gzLevel                = gzLevel;
    S.keepOrig               = keepOrig;
    S.umif                   = umif;
    S.umir                   = umir;
    S.umiTf                  = umiTf;
    S.umiTr                  = umiTr;
    S.detectAdapt            = detectAdapt;
    S.trimKey                = trimKey;
    S.key1                   = key1;
    S.key2                   = key2;
    S.printLog               = printLog;
    S.logFileName            = logFileName;
    S.allowAligned           = allowAligned;
    S.produceUnCompressedBAM = produceUnCompressedBAM;
    S.bamFile                = bamFile;
    S.bamFileOUT             = bamFileOUT;
    S.fastqToBam             = fastqToBam;
    S.fqQualOffset           = fqQualOffset;
    S.singleEndModeFQ        = singleEndModeFQ;
    S.fastqfile1             = fastqfile1;
    S.fastqfile2             = fastqfile2;
    S.fastqoutfile           = fastqoutfile;
    S.fastqoutfileb          = fastqoutfileb;
    S.fastqoutfilese         = fastqoutfilese;
    S.fastqoutfilesef        = fastqoutfilesef;
    S.fastqoutfilepe1        = fastqoutfilepe1;
    S.fastqoutfilepe1f       = fastqoutfilepe1f;
    S.fastqoutfilepe2        = fastqoutfilepe2;
    S.fastqoutfilepe2f       = fastqoutfilepe2f;
    S.intfastqoutfile        = intfastqoutfile;
    for(int k=0;k<6;k++)
	S.sidx[k] = k;

    if(fastqFormat && !fastqToBam)
	return runFastq(S);

    return runBam(S, argc, argv);
}
