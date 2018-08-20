#include <iostream>
#include <string>
#include <cstring>
#include <sys/stat.h>
#include <gzstream.h>

#include <api/SamHeader.h>
#include <api/BamMultiReader.h>
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <api/BamAux.h>

#include "MergeTrimReads.h"
#include "PutProgramInHeader.h"
#include "FastQParser.h"

// #include "JSON.h"

#include "utils.h"

////////////////////////////////
// TO DO
//
////////////////////////////////
using namespace std;

typedef struct { 
    ogzstream single;
    ogzstream pairr1;
    ogzstream pairr2;

    ogzstream singlef;
    ogzstream pairr1f;
    ogzstream pairr2f;

} fqwriters;


//using namespace BamTools;
// using namespace __MergeTrimReads__;



// void initializeDefaultSequences(string configFile){
//     string line;
//     ifstream myFile;
//     string content="";
//     myFile.open(configFile.c_str(), ios::in);

//     if (myFile.is_open()){
// 	while ( getline (myFile,line)){
// 	    content+=line;
// 	}
// 	myFile.close();
//     }else{
// 	cerr << "Unable to open config file "<<configFile<<endl;
// 	exit(1);
//     }


//     JSONValue *value = JSON::Parse(content.c_str());
//     if (value == NULL){
// 	cerr<<"Failed to parse JSON file"<<endl;
// 	exit(1);
//     }

    
// }
string options_adapter_F_BAM="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACIIIIIIIATCTCGTATGCCGTCTTCTGCTTG";
string options_adapter_S_BAM="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT";
string options_adapter_chimera_BAM="ACACTCTTTCCCTACACGTCTGAACTCCAG,ACACTCTTTCCCACACGTCTGAACTCCAGT,ACACTCTTTCCCTACACACGTCTGAACTCC,CTCTTTCCCTACACGTCTGAACTCCAGTCA,GAAGAGCACACGTCTGAACTCCAGTCACII,GAGCACACGTCTGAACTCCAGTCACIIIII,GATCGGAAGAGCACACGTCTGAACTCCAGT,AGATCGGAAGAGCACACGTCTGAACTCCAG,AGAGCACACGTCTGAACTCCAGTCACIIII,ACACGTCTGAACTCCAGTCACIIIIIIIAT,GTGCACACGTCTGAACTCCAGTCACIIIII,AGCACACGTCTGAACTCCAGTCACIIIIII,CGTATGCCGTCTTCTGCTTGAAAAAAAAAA";

inline bool isBamAlignEmpty(const BamAlignment & toTest){
    return ( toTest.Name.empty() &&
	     toTest.Length == 0 );    
}


int main (int argc, char *argv[]) {

    bool produceUnCompressedBAM=false;
    bool verbose=false;
    bool ancientDNA=false;
    bool keepOrig=false;

    string adapter_F=options_adapter_F_BAM;
    string adapter_S=options_adapter_S_BAM;
    string adapter_chimera=options_adapter_chimera_BAM;
    string key="";
    bool   trimKey=false;

    bool allowMissing=false;
    int trimCutoff=1;

    bool allowAligned=false;
    bool printLog=false;
    string logFileName;

    BamReader reader;
    BamWriter writer;

    string bamFile;
    string bamFileOUT="";

    string key1;
    string key2;
    
    bool useDist=false;
    double location=-1.0;
    double scale   =-1.0;

    bool fastqFormat=false;
    string fastqfile1   = "";
    string fastqfile2   = "";
    string fastqoutfile = "";
    bool singleEndModeFQ=true;
    int    qualOffset     = 33;

    const string usage=string(string(argv[0])+
			      
    			      " [options] BAMfile"+"\n"+
    			      "\nThis program takes an unaligned BAM where mates are consecutive\nor fastq files and trims and merges reads\n"+

    			      "\n\tYou can specify a unaligned bam file or one or two fastq :\n"+			      
    			      "\t\t"+"-fq1" +"\t\t"+"First fastq"+"\n"+
    			      "\t\t"+"-fq2" +"\t\t"+"Second  fastq file (for paired-end)"+"\n"+
    			      "\t\t"+"-fqo" +"\t\t"+"Output fastq prefix"+"\n\n"+
    			      //"\t"+"-p , --PIPE"+"\n\t\t"+"Read BAM from and write it to PIPE"+"\n"+
    			      "\t"+"-o , --outfile" +"\t\t"+"Output (BAM format)."+"\n"+


    			      "\t"+"-u            " +"\t\t"+"Produce uncompressed bam (good for pipe)"+"\n"+

    			      //	"\t"+" , --outprefix" +"\n\t\t"+"Prefix for output files (default '"+outprefix+"')."+"\n"+
    			      //"\t"+" , --SAM" +"\n\t\t"+"Output SAM not BAM."+"\n"+
    			      "\t"+"--aligned" +"\t\t"+"Allow reads to be aligned (default "+boolStringify(allowAligned)+")"+"\n"+
    			      "\t"+"-v , --verbose" +"\t\t"+"Turn all messages on (default "+boolStringify(verbose)+")"+"\n"+
    			      "\t"+"--log [log file]" +"\t"+"Print a tally of merged reads to this log file (default only to stderr)"+"\n"+
    			      "\t"+"--phred64" +"\t\t"+"Use PHRED 64 as the offset for QC scores (default : PHRED33)"+"\n"+
			      
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
    			      "\t\t"+"-f , --adapterFirstRead" +"\t\t\t"+"Adapter that is observed after the forward read (def. Multiplex: "+options_adapter_F_BAM.substr(0,30)+")"+"\n"+
    			      "\t\t"+"-s , --adapterSecondRead" +"\t\t"+"Adapter that is observed after the reverse read (def. Multiplex: "+options_adapter_S_BAM.substr(0,30)+")"+"\n"+
    			      "\t\t"+"-c , --FirstReadChimeraFilter" +"\t\t"+"If the forward read looks like this sequence, the cluster is filtered out.\n\t\t\t\t\t\t\tProvide several sequences separated by comma (def. Multiplex: "+options_adapter_chimera_BAM.substr(0,30)+")"+"\n"+
    			      "\t\t"+"-k , --key"+"\t\t\t\t"+"Key sequence with which each sequence starts. Comma separate for forward and reverse reads. (default '"+key+"')"+"\n"+
			      "\t\t"+"--trimkey"     +"\t\t\t\t"+"Trim the key sequence even for untrimmed. (default '"+boolStringify(trimKey)+"')"+"\n"+

			      
    			      "\t\t"+"-i , --allowMissing"+"\t\t\t"+"Allow one base in one key to be missing or wrong. (default "+boolStringify(allowMissing)+")"+"\n"+
    			      "\t\t"+"--trimCutoff"+"\t\t\t"+"Lowest number of adapter bases to be observed for single Read trimming (default "+stringify(trimCutoff)+")");


    if( (argc== 1) ||
    	(argc== 2 && string(argv[1]) == "-h") ||
    	(argc== 2 && string(argv[1]) == "-help") ||
    	(argc== 2 && string(argv[1]) == "--help") ){
    	cout<<"Usage:"<<endl;
    	cout<<""<<endl;
    	cout<<usage<<endl;
    	return 1;
    }

    

    for(int i=1;i<(argc-1);i++){ //all but the last arg

    	if(strcmp(argv[i],"-fq1") == 0 ){
    	    fastqfile1=string(argv[i+1]);
    	    fastqFormat=true;
    	    i++;
    	    continue;
    	}

    	if(strcmp(argv[i],"-fq2") == 0 ){
    	    fastqfile2=string(argv[i+1]);
    	    fastqFormat=true;
    	    singleEndModeFQ=false;
    	    i++;
    	    continue;
    	}

    	if(strcmp(argv[i],"-fqo") == 0 ){
    	    fastqoutfile=string(argv[i+1]);
    	    fastqFormat=true;
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

    	if(string(argv[i]) == "--phred64"  ){
    	    qualOffset=64;
    	    continue;
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
    	    i++;
    	    continue;
    	}


    	if(strcmp(argv[i],"-s") == 0 || strcmp(argv[i],"--adapterSecondRead") == 0 ){
    	    adapter_S =string(argv[i+1]);
    	    i++;
    	    continue;
    	}

    	if(strcmp(argv[i],"-c") == 0 || strcmp(argv[i],"--FirstReadChimeraFilter") == 0 ){
    	    adapter_chimera =string(argv[i+1]);
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

    	if(strcmp(argv[i],"--trimCutoff") == 0 ){
    	    trimCutoff=atoi(argv[i+1]);
    	    i++;
    	    continue;
    	}
	
    	cerr<<"Unknown option "<<argv[i] <<" exiting"<<endl;
    	return 1;	    
    }

    bamFile=argv[argc-1];

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


    MergeTrimReads mtr (adapter_F,adapter_S,adapter_chimera,
    			key1,key2,
    			trimCutoff,allowMissing,ancientDNA,location,scale,useDist,qualOffset);

    fqwriters onereadgroup;

    if(fastqFormat){
	
    	if( bamFileOUT != ""  || produceUnCompressedBAM || allowAligned){
    	    cerr<<"ERROR : Cannot specify options like -o, -u or --allowAligned for fastq"<<endl;
    	    return 1;
    	}

    	if(fastqfile1 == ""){
    	    cerr<<"ERROR : Must specify as least the first file for fastq"<<endl;
    	    return 1;	    
    	}



    	FastQParser * fqp1;
    	FastQParser * fqp2;

    	if(singleEndModeFQ){
    	    fqp1 = new FastQParser (fastqfile1);

    	    string outdirs   = fastqoutfile+".fq.gz";
    	    string outdirsf  = fastqoutfile+".fail.fq.gz";

    	    onereadgroup.single.open(outdirs.c_str(), ios::out);
    	    onereadgroup.singlef.open(outdirsf.c_str(), ios::out);

    	    if(!onereadgroup.single.good()){       cerr<<"Cannot write to file "<<outdirs<<endl; return 1; }
    	    if(!onereadgroup.singlef.good()){      cerr<<"Cannot write to file "<<outdirsf<<endl; return 1; }

	    
    	}else{
    	    fqp1 = new FastQParser (fastqfile1);
    	    fqp2 = new FastQParser (fastqfile2);

    	    string outdirs   = fastqoutfile+".fq.gz";
    	    string outdir1   = fastqoutfile+"_r1.fq.gz";
    	    string outdir2   = fastqoutfile+"_r2.fq.gz";

    	    string outdirsf  = fastqoutfile+".fail.fq.gz";
    	    string outdir1f  = fastqoutfile+"_r1.fail.fq.gz";
    	    string outdir2f  = fastqoutfile+"_r2.fail.fq.gz";

    	    onereadgroup.single.open(outdirs.c_str(), ios::out);
    	    onereadgroup.pairr1.open(outdir1.c_str(), ios::out);
    	    onereadgroup.pairr2.open(outdir2.c_str(), ios::out);

    	    onereadgroup.singlef.open(outdirsf.c_str(), ios::out);
    	    onereadgroup.pairr1f.open(outdir1f.c_str(), ios::out);
    	    onereadgroup.pairr2f.open(outdir2f.c_str(), ios::out);

    	    if(!onereadgroup.single.good()){       cerr<<"Cannot write to file "<<outdirs<<endl; return 1; }
    	    if(!onereadgroup.pairr1.good()){       cerr<<"Cannot write to file "<<outdir1<<endl; return 1; }
    	    if(!onereadgroup.pairr2.good()){       cerr<<"Cannot write to file "<<outdir2<<endl; return 1; }
	    
    	    if(!onereadgroup.singlef.good()){      cerr<<"Cannot write to file "<<outdirsf<<endl; return 1; }
    	    if(!onereadgroup.pairr1f.good()){      cerr<<"Cannot write to file "<<outdir1f<<endl; return 1; }
    	    if(!onereadgroup.pairr2f.good()){      cerr<<"Cannot write to file "<<outdir2f<<endl; return 1; }
	    
    	}


    	unsigned int totalSeqs=0;
    	while(fqp1->hasData()){

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

    		merged result=	mtr.process_PE(*(fo1->getSeq()),*(fo1->getQual()),
    					       *(fo2->getSeq()),*(fo2->getQual()));

    		mtr.incrementCountall();

    		if(result.code != ' '){ //keys or chimeras

    		    if(result.code == 'K'){
    			mtr.incrementCountfkey();
    		    }else{
    			if(result.code == 'D'){
    			    mtr.incrementCountchimera();
    			}else{
    			    cerr << "leehom: Wrong return code =\""<<result.code<<"\""<<endl;
    			    exit(1);
    			}
    		    }
			
    		    onereadgroup.pairr2f<<"@"<<def2s<<"/2" <<endl <<*(fo2->getSeq())<<endl<<"+"<<endl <<*(fo2->getQual())<<endl;
    		    onereadgroup.pairr1f<<"@"<<def1s<<"/1" <<endl <<*(fo1->getSeq())<<endl<<"+"<<endl <<*(fo1->getQual())<<endl;
    		    continue;

    		}else{
    		        if(result.sequence != ""){ //new sequence			    
    			    onereadgroup.single<<"@"<<def1s<<"" <<endl << result.sequence<<endl<<"+"<<endl <<result.quality<<endl;    	    

    			    if( result.sequence.length() > max(fo1->getSeq()->length(),fo2->getSeq()->length()) ){
    				mtr.incrementCountmergedoverlap();
    			    }else{
    				mtr.incrementCountmerged();			  
    			    }

    			}else{ //keep as is
    			    mtr.incrementCountnothing();

			    if(trimKey){
				onereadgroup.pairr2<<"@"<<def2s<<"/2" <<endl << (fo2->getSeq())->substr( key2.length() ) <<endl<<"+"<<endl <<(fo2->getQual())->substr( key2.length() )<<endl;
				onereadgroup.pairr1<<"@"<<def1s<<"/1" <<endl << (fo1->getSeq())->substr( key1.length() ) <<endl<<"+"<<endl <<(fo1->getQual())->substr( key1.length() )<<endl;
			    }else{
				onereadgroup.pairr2<<"@"<<def2s<<"/2" <<endl <<*(fo2->getSeq())                          <<endl<<"+"<<endl <<*(fo2->getQual())                        `<<endl;
				onereadgroup.pairr1<<"@"<<def1s<<"/1" <<endl <<*(fo1->getSeq())                          <<endl<<"+"<<endl <<*(fo1->getQual())                        <<endl;
			    }

    			}
    		}

    	    }else{
		
		
    		merged result=mtr.process_SR(*(fo1->getSeq()),*(fo1->getQual()));
    		mtr.incrementCountall();

    		if(result.code != ' '){ //either chimera or missing key

    		    if(result.code == 'K'){
    			mtr.incrementCountfkey();
    		    }else{
    			if(result.code == 'D'){
    			    mtr.incrementCountchimera();
    			}else{
    			    cerr << "leehom: Wrong return code =\""<<result.code<<"\""<<endl;
    			    exit(1);
    			}
    		    }

    		    onereadgroup.singlef<<""<<*(fo1->getID())<<"" <<endl << *(fo1->getSeq())<<endl<<"+"<<endl <<*(fo1->getQual())<<endl;
    		    continue;
    		}

    		if(result.sequence != ""){ //new sequence
    		    mtr.incrementCounttrimmed();
    		    onereadgroup.single<<""<<*(fo1->getID())<<"" <<endl << result.sequence<<endl<<"+"<<endl <<result.quality<<endl;
    		}else{
    		    mtr.incrementCountnothing();
		    if(trimKey){
			onereadgroup.single<<""<<(fo1->getID())<<"" << endl <<  (fo1->getSeq())->substr( key1.length() )<<endl<<"+"<<endl << (fo1->getQual())->substr( key1.length() )<<endl;
		    }else{
			onereadgroup.single<<""<<*(fo1->getID())<<"" <<endl << *(fo1->getSeq())                         <<endl<<"+"<<endl <<*(fo1->getQual())                         <<endl;
		    }
    		}

    	    }
	


    	    totalSeqs++;
    	}
    
    	delete fqp1;
    	if(!singleEndModeFQ){
    	    delete fqp2;
    	}

    	if(singleEndModeFQ){

    	    onereadgroup.single.close();
    	    onereadgroup.singlef.close();
	    
    	}else{
    	    onereadgroup.single.close();
    	    onereadgroup.pairr1.close();
    	    onereadgroup.pairr2.close();
	    
    	    onereadgroup.singlef.close();
    	    onereadgroup.pairr1f.close();
    	    onereadgroup.pairr2f.close();
    	}
    
    	//fastq
    }else{
    	//else BAM


    	//  initMerge();
    	//     set_adapter_sequences(adapter_F,
    	// 			  adapter_S,
    	// 			  adapter_chimera);
    	//     set_options(trimCutoff,allowMissing,mergeoverlap);






    	if( bamFileOUT == ""  ){
    	    cerr<<"The output must be a be specified, exiting"<<endl;
    	    return 1;
    	}

    	if ( !reader.Open(bamFile) ) {
    	    cerr << "Could not open input BAM file  "<<bamFile << endl;
    	    return 1;
    	}
    	SamHeader header = reader.GetHeader();

    

    	string pID          = "leeHom";   
    	string pName        = "leeHom";   
    	string pCommandLine = "";
    	for(int i=0;i<(argc);i++){
    	    pCommandLine += (string(argv[i])+" ");
    	}
    	putProgramInHeader(&header,pID,pName,pCommandLine,returnGitHubVersion(string(argv[0]),".."));

    	const RefVector references = reader.GetReferenceData();
    	//we will not call bgzip with full compression, good for piping into another program to 
    	//lessen the load on the CPU
    	if(produceUnCompressedBAM) 
    	    writer.SetCompressionMode(BamWriter::Uncompressed);

    	if ( !writer.Open(bamFileOUT,header,references) ) {
    	    cerr << "Could not open output BAM file "<<bamFileOUT << endl;
    	    return 1;
    	}



    	SamHeader sh=reader.GetHeader();
    	//Up to the user to be sure that a sequence is followed by his mate
    	// if(!sh.HasSortOrder() || 
    	//    sh.SortOrder != "queryname"){
    	// 	cerr << "Bamfile must be sorted by queryname" << endl;
    	// 	return 1;
    	// }
    

    	BamAlignment al;
    	BamAlignment orig;
    	BamAlignment orig2;

    	BamAlignment al2;
    	bool al2Null=true;
    
    	while ( reader.GetNextAlignment(al) ) {


    	    if(al.IsMapped() || al.HasTag("NM") || al.HasTag("MD")  ){
    		if(!allowAligned){
    		    cerr << "Reads should not be aligned" << endl;
    		    return 1;
    		}else{
    		    //should we remove tags ?
    		}
    	    }


    	    if(al.IsPaired() && 
    	       al2Null ){
    		al2=al;
    		al2Null=false;
    		continue;
    	    }else{
    		if(al.IsPaired() && 
    		   !al2Null){
    		    if(keepOrig){
    			orig = al;
    		    }

    		    bool  result =  mtr.processPair(al,al2);


    		    if( result ){//was merged

    			if(keepOrig){
    			    orig2 = al2;
			    al.Name = al.Name+"_M";
    			}

    			writer.SaveAlignment(al);

    			if(keepOrig){
    			    orig.SetIsDuplicate(true);
    			    orig2.SetIsDuplicate(true);
    			    writer.SaveAlignment(orig2);
    			    writer.SaveAlignment(orig);
    			}

    			//the second record is empty
    		    }else{
    			//keep the sequences as pairs

    			writer.SaveAlignment(al2);		    
    			writer.SaveAlignment(al);
    		    }

    		    //
    		    //  SINGLE END
    		    //
    		}else{ 

    		    if(keepOrig){
    			orig =al;
			al.Name = al.Name+"_M";
    		    }
    		    mtr.processSingle(al);

    		    if(keepOrig){
    			//write duplicate
    			if(orig.QueryBases.length()  != al.QueryBases.length()){
    			    orig.SetIsDuplicate(true);
    			    writer.SaveAlignment(orig);
    			}
    		    }
    		    writer.SaveAlignment(al);



    		} //end single end
    		al2Null=true;
    	    }//second pair
		    

    	} //while al
    	reader.Close();
    	writer.Close();


    } //else BAM


    cerr <<mtr.reportSingleLine()<<endl;

    if(printLog){
    	ofstream fileLog;
    	fileLog.open(logFileName.c_str());

    	if (fileLog.is_open()){
    	    fileLog <<mtr.reportMultipleLines() <<endl;

    	}else{
    	    cerr << "Unable to print to file "<<logFileName<<endl;
    	}
    	fileLog.close();
    }
    return 0;
}

