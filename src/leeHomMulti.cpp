/*
 * threadsFastq
 * Date: Feb-27-2016 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>
#include <queue>

#include <api/SamHeader.h>
#include <api/BamMultiReader.h>
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <api/BamAux.h>
#include "MergeTrimReads.h"
#include "PutProgramInHeader.h"
#include "FastQParser.h"

#include "utils.h"

//#define DEBUG 


//#define DEBUGMEM

#ifdef DEBUG
#define DEBUGSLEEPING
#endif

//#define DEBUGSLEEPING
using namespace std;
using namespace BamTools;

typedef struct { 
    ogzstream single;
    ogzstream pairr1;
    ogzstream pairr2;

    ogzstream singlef;
    ogzstream pairr1f;
    ogzstream pairr2f;

} fqwriters;


struct argspthread {
    bool  keepOrig;
    bool  singleEndModeFQ;
    unsigned int umif;
    unsigned int umir;
};

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

pthread_mutex_t  mutexCERR    = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t  mutexQueue   = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t  mutexCounter = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t  mutexRank    = PTHREAD_MUTEX_INITIALIZER;


//! Chunk of code to check if a certain thread call failed
/*!
  This block is calls by the pthread


*/				
#define checkResults(string, val) {             \
 if (val) {                                     \
     cerr<<"Failed with "<<val<<" at "<<string<<endl;	\
   exit(1);                                     \
 }                                              \
}
 


inline bool isBamAlignEmpty(const BamAlignment & toTest){
    return ( toTest.Name.empty() &&
	     toTest.Length == 0 );    
}




//BAM
class DataChunk{
private:
    
public:
    vector<BamAlignment>  * dataToProcess;    
    int rank;

    // DataChunk();
    DataChunk(int rank_m);
    DataChunk(const DataChunk & other);
    ~DataChunk();
    DataChunk & operator= (const DataChunk & other);
};


DataChunk::DataChunk(int rank_m){
    rank =  rank_m;
    dataToProcess = new vector<BamAlignment>() ;    

#ifdef DEBUGMEM
    int rc = pthread_mutex_lock(&mutexCERR);
    checkResults("pthread_mutex_lock()\n", rc);

    cerr<<"Constructor addr: "<<this<<" rank "<<rank<<endl;

    rc = pthread_mutex_unlock(&mutexCERR);
    checkResults("pthread_mutex_unlock()\n", rc);
#endif

}


DataChunk::~DataChunk(){
    delete dataToProcess;

#ifdef DEBUGMEM
    int rc = pthread_mutex_lock(&mutexCERR);
    checkResults("pthread_mutex_lock()\n", rc);
    
    cerr<<"Destructor  addr: "<<this<<" rank "<<rank<<endl;

    rc = pthread_mutex_unlock(&mutexCERR);
    checkResults("pthread_mutex_unlock()\n", rc);
#endif

}

class CompareDataChunk {
public:
    bool operator() ( DataChunk * cd1, DataChunk * cd2)  {
        //comparison code here
	return ( cd1->rank > cd2->rank );
    }
};


queue< DataChunk * >                                               queueDataToprocess;
priority_queue<DataChunk *, vector<DataChunk *>, CompareDataChunk > queueDataTowrite;

//fq
class fqrecord{
 public:
    char   code;
    //merged
    string sequence;
    string quality;

    string d1;
    string s1;
    string q1;
    string d2;
    string s2;
    string q2;
    string cmt;
    

    bool paired;
};

class DataChunkFQ{
private:
    
public:
    vector<fqrecord>  * dataToProcess;    
    int rank;

    // DataChunk();
    DataChunkFQ(int rank_m);
    DataChunkFQ(const DataChunkFQ & other);
    ~DataChunkFQ();
    DataChunkFQ & operator= (const DataChunkFQ & other);
};


DataChunkFQ::DataChunkFQ(int rank_m){
    rank =  rank_m;
    dataToProcess = new vector<fqrecord>() ;    

#ifdef DEBUGMEM
    int rc = pthread_mutex_lock(&mutexCERR);
    checkResults("pthread_mutex_lock()\n", rc);

    cerr<<"Constructor addr: "<<this<<" rank "<<rank<<endl;

    rc = pthread_mutex_unlock(&mutexCERR);
    checkResults("pthread_mutex_unlock()\n", rc);
#endif

}


DataChunkFQ::~DataChunkFQ(){
    delete dataToProcess;

#ifdef DEBUGMEM
    int rc = pthread_mutex_lock(&mutexCERR);
    checkResults("pthread_mutex_lock()\n", rc);

    cerr<<"Destructor  addr: "<<this<<" rank "<<rank<<endl;    

    rc = pthread_mutex_unlock(&mutexCERR);
    checkResults("pthread_mutex_unlock()\n", rc);
#endif

}

class CompareDataChunkFQ {
public:
    bool operator() ( DataChunkFQ * cd1, DataChunkFQ * cd2)  {
        //comparison code here
	return ( cd1->rank > cd2->rank );
    }
};


queue< DataChunkFQ * >                                                  queueDataToprocessFQ;
priority_queue<DataChunkFQ *, vector<DataChunkFQ *>, CompareDataChunkFQ > queueDataTowriteFQ;





//GLOBALLY accessed
int maxQueueDataTowriteSize  =10;
int maxQueueDataToprocessSize=10;

int    timeThreadSleep = 1;
bool      readDataDone = false;
unsigned int sizeChunk = 5000;

int numberOfThreads=1;

map<unsigned int, int>       threadID2Rank;






void forceInsertProcedure(DataChunk * chunkToWrite,int rankThread){
    int   rc;
    //---> mutex counter is taken here <---

 forceinsertinit:
    //need to insert chunkToWrite
    
    if(int(queueDataTowrite.size())>(maxQueueDataTowriteSize)){	  // if queue full
	//---> mutex counter is taken here <---

	if(chunkToWrite->rank <  queueDataTowrite.top()->rank){  //if my current chunk out ranks the first one, need to insert
	    
#ifdef DEBUGSLEEPING
	    rc = pthread_mutex_lock(&mutexCERR);
	    checkResults("pthread_mutex_lock()\n", rc);	       
	    cerr<<"Thread #"<<rankThread<<" queue is full with "<<int(queueDataTowrite.size())<<" force insert, current rank "<<chunkToWrite->rank<<" top rank "<<queueDataTowrite.top()->rank<<", sleeping for "<<timeThreadSleep<<"s"<<endl;
	    rc = pthread_mutex_unlock(&mutexCERR);
	    checkResults("pthread_mutex_unlock()\n", rc);     
#endif

	    DataChunk *  tempDataChunk = queueDataTowrite.top();
	    queueDataTowrite.pop();
	    queueDataTowrite.push(chunkToWrite);

	    //release mutex and sleep
	    rc = pthread_mutex_unlock(&mutexCounter);
	    checkResults("pthread_mutex_unlock()\n", rc);	
		
	    sleep(timeThreadSleep);
	    
	    chunkToWrite=tempDataChunk;

	    //take mutex
	    rc = pthread_mutex_lock(&mutexCounter);
	    checkResults("pthread_mutex_lock()\n", rc);
	    goto forceinsertinit;
	}else{//if the currentChunk does not outrank the one on top

#ifdef DEBUGSLEEPING
	    rc = pthread_mutex_lock(&mutexCERR);
	    checkResults("pthread_mutex_lock()\n", rc);
	       
	    cerr<<"Thread #"<<rankThread<<" queue is full with "<<int(queueDataTowrite.size())<<" items, but does not outrank, current rank "<<chunkToWrite->rank<<" top rank "<<queueDataTowrite.top()->rank<<" sleeping for "<<timeThreadSleep<<endl;

	    rc = pthread_mutex_unlock(&mutexCERR);
	    checkResults("pthread_mutex_unlock()\n", rc);
#endif

	    //release mutex
	    rc = pthread_mutex_unlock(&mutexCounter);
	    checkResults("pthread_mutex_unlock()\n", rc);	
		
	    sleep(timeThreadSleep);
	    //take mutex
	    rc = pthread_mutex_lock(&mutexCounter);
	    checkResults("pthread_mutex_lock()\n", rc);
	    goto forceinsertinit;
	}

    }else{

#ifdef DEBUG 
	rc = pthread_mutex_lock(&mutexCERR);
	checkResults("pthread_mutex_lock()\n", rc);
	cerr<<"Thread #"<<rankThread<<" force insert, enough space size= "<<int(queueDataTowrite.size())<<" inserting"<<endl;
	rc = pthread_mutex_unlock(&mutexCERR);
	checkResults("pthread_mutex_unlock()\n", rc);
#endif

	//we now insert it
	
	queueDataTowrite.push(chunkToWrite);

	//release mutex
	rc = pthread_mutex_unlock(&mutexCounter);
	checkResults("pthread_mutex_unlock()\n", rc);	
	return ;
    }
       
}













void forceInsertProcedureFQ(DataChunkFQ * chunkToWrite,int rankThread){
    int   rc;
    //---> mutex counter is taken here <---

 forceinsertinit:
    //need to insert chunkToWrite
    
    if(int(queueDataTowriteFQ.size())>(maxQueueDataTowriteSize)){	  // if queue full
	//---> mutex counter is taken here <---

	if(chunkToWrite->rank <  queueDataTowriteFQ.top()->rank){  //if my current chunk out ranks the first one, need to insert
	    
#ifdef DEBUGSLEEPING
	    rc = pthread_mutex_lock(&mutexCERR);
	    checkResults("pthread_mutex_lock()\n", rc);
	    cerr<<"Thread #"<<rankThread<<" queue is full with "<<int(queueDataTowriteFQ.size())<<" force insert, current rank "<<chunkToWrite->rank<<" top rank "<<queueDataTowriteFQ.top()->rank<<", sleeping for "<<timeThreadSleep<<"s"<<endl;
	    rc = pthread_mutex_unlock(&mutexCERR);
	    checkResults("pthread_mutex_unlock()\n", rc);     
#endif

	    DataChunkFQ *  tempDataChunk = queueDataTowriteFQ.top();
	    queueDataTowriteFQ.pop();
	    queueDataTowriteFQ.push(chunkToWrite);

	    //release mutex and sleep
	    rc = pthread_mutex_unlock(&mutexCounter);
	    checkResults("pthread_mutex_unlock()\n", rc);	
		
	    sleep(timeThreadSleep);
	    
	    chunkToWrite=tempDataChunk;

	    //take mutex
	    rc = pthread_mutex_lock(&mutexCounter);
	    checkResults("pthread_mutex_lock()\n", rc);
	    goto forceinsertinit;
	}else{//if the currentChunk does not outrank the one on top

#ifdef DEBUGSLEEPING
	    rc = pthread_mutex_lock(&mutexCERR);
	    checkResults("pthread_mutex_lock()\n", rc);
	       
	    cerr<<"Thread #"<<rankThread<<" queue is full with "<<int(queueDataTowriteFQ.size())<<" items, but does not outrank, current rank "<<chunkToWrite->rank<<" top rank "<<queueDataTowriteFQ.top()->rank<<" sleeping for "<<timeThreadSleep<<"s"<<endl;

	    rc = pthread_mutex_unlock(&mutexCERR);
	    checkResults("pthread_mutex_unlock()\n", rc);
#endif

	    //release mutex
	    rc = pthread_mutex_unlock(&mutexCounter);
	    checkResults("pthread_mutex_unlock()\n", rc);	
		
	    sleep(timeThreadSleep);

	    //take mutex
	    rc = pthread_mutex_lock(&mutexCounter);
	    checkResults("pthread_mutex_lock()\n", rc);
	    goto forceinsertinit;
	}

    }else{

#ifdef DEBUG 
	rc = pthread_mutex_lock(&mutexCERR);
	checkResults("pthread_mutex_lock()\n", rc);
	cerr<<"Thread #"<<rankThread<<" force insert, enough space size= "<<int(queueDataTowriteFQ.size())<<" inserting"<<endl;
	rc = pthread_mutex_unlock(&mutexCERR);
	checkResults("pthread_mutex_unlock()\n", rc);
#endif

	//we now insert it
	
	queueDataTowriteFQ.push(chunkToWrite);

	//release mutex
	rc = pthread_mutex_unlock(&mutexCounter);
	checkResults("pthread_mutex_unlock()\n", rc);	
	return ;
    }
       
}






//! Method called for each thread BAM
/*!
  

*/				
void *mainComputationThread(void * argc){    
    //bool * keepOrigAddr = (bool *)argc;
    //bool   keepOrig     = *keepOrigAddr;

    bool keepOrig = ((struct argspthread*)argc)->keepOrig;
    unsigned int umif      = ((struct argspthread*)argc)->umif;
    unsigned int umir      = ((struct argspthread*)argc)->umir;
    bool umi      = ((umif+umir)>0);
    
    int   rc;
    int rankThread=0;

    rc = pthread_mutex_lock(&mutexRank);
    checkResults("pthread_mutex_lock()\n", rc);

    threadID2Rank[*(int *)pthread_self()]  = threadID2Rank.size()+1;
    rankThread = threadID2Rank[*(int *)pthread_self()];

    
    rc = pthread_mutex_unlock(&mutexRank);
    checkResults("pthread_mutex_unlock()\n", rc);

 checkqueue:    
    // stackIndex=-1;
    //check stack

    
    rc = pthread_mutex_lock(&mutexQueue);
    checkResults("pthread_mutex_lock()\n", rc);


    bool foundData=false;
    
    //cerr<<"Thread #"<<rankThread <<" started and is requesting data"<<endl;


    DataChunk * currentChunk;


    if(!queueDataToprocess.empty()){    
 	foundData=true;
 	currentChunk = queueDataToprocess.front();
 	queueDataToprocess.pop();
#ifdef DEBUG 

	rc = pthread_mutex_lock(&mutexCERR);
	checkResults("pthread_mutex_lock()\n", rc);
	
	
	cerr<<"Thread #"<<rankThread<<" is reading "<<currentChunk->rank<<endl;
	
	rc = pthread_mutex_unlock(&mutexCERR);
	checkResults("pthread_mutex_unlock()\n", rc);

#endif
	//cout<<"rank "<< &(currentChunk->dataToProcess) <<endl;
    }

    
  

    if(!foundData){
	rc = pthread_mutex_unlock(&mutexQueue);
	checkResults("pthread_mutex_unlock()\n", rc);


	if(readDataDone){
#ifdef DEBUG 
	rc = pthread_mutex_lock(&mutexCERR);
	checkResults("pthread_mutex_lock()\n", rc);
	
	cerr<<"Thread #"<<rankThread<<" is done"<<endl;

	
	rc = pthread_mutex_unlock(&mutexCERR);
	checkResults("pthread_mutex_unlock()\n", rc);

	
#endif
	    return NULL;	
	}else{

#ifdef DEBUGSLEEPING
	rc = pthread_mutex_lock(&mutexCERR);
	checkResults("pthread_mutex_lock()\n", rc);	       
	cerr<<"Thread #"<<rankThread<<" has not found data, sleeping for "<<timeThreadSleep<<endl;	
	rc = pthread_mutex_unlock(&mutexCERR);
	checkResults("pthread_mutex_unlock()\n", rc);
#endif

	    sleep(timeThreadSleep);
	    goto checkqueue;
	}

    }else{
	//release stack
	rc = pthread_mutex_unlock(&mutexQueue);
	checkResults("pthread_mutex_unlock()\n", rc);
    }

    //////////////////////////////////////////////////////////////
    //                BEGIN COMPUTATION                         //
    //////////////////////////////////////////////////////////////
    DataChunk * chunkToWrite = new DataChunk(currentChunk->rank);

    BamAlignment al;
    BamAlignment al2;
    BamAlignment orig;
    BamAlignment orig2;
    bool al2Null=true;
    string umifS;
    string umirS;
    string umifQ;
    string umirQ;
    
    
    for(unsigned i=0;i<currentChunk->dataToProcess->size();i++){
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
	
	al = currentChunk->dataToProcess->at(i);
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

		bool  result =  mtr->processPair(al,al2);

		if( result ){//was merged

		    if(keepOrig){
			orig2 = al2;
			al.Name = al.Name+"_M";
		    }

		    if( umi ){

			if( al.QueryBases.length() > (umif+umir)){
			    if(umif != 0){	    
				umifS         = al.QueryBases.substr(0,umif);
				umifQ         = al.Qualities.substr( 0,umif);

				al.QueryBases = al.QueryBases.substr(umif);
				al.Qualities  = al.Qualities.substr( umif);

				if(umir == 0){	//no reverse UMI, add the tags now
				    if(!al.AddTag("RX","Z",umifS)){ cerr << "Error while editing RX tags for UMI:"<<umifS<<"#"<< endl; exit(1);   }				    
				    if(!al.AddTag("QX","Z",umifQ)){ cerr << "Error while editing QX tags for UMI:"<<umifQ<<"#"<< endl; exit(1);   }
				}

			    }

			    if(umir != 0){
				umirS         = al.QueryBases.substr(al.QueryBases.length()-umir,  umir);
				umirQ         = al.Qualities.substr( al.Qualities.length() -umir,  umir);
			    
				al.QueryBases = al.QueryBases.substr(0,al.QueryBases.length()-umir);
				al.Qualities  = al.Qualities.substr(0, al.Qualities.length() -umir);

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
		    
		    chunkToWrite->dataToProcess->push_back(al);
		    //writer.SaveAlignment(al);
		    
		    if(keepOrig){
			orig.SetIsDuplicate(true);
			orig2.SetIsDuplicate(true);
			// writer.SaveAlignment(orig2);
			// writer.SaveAlignment(orig);
			chunkToWrite->dataToProcess->push_back(orig2);
			chunkToWrite->dataToProcess->push_back(orig);
		    }

		    //the second record is empty
		}else{
		    //keep the sequences as pairs

		    
		    if(umif != 0){	    
			//al2 is the fwd read
			umifS          = al2.QueryBases.substr(0,umif);
			umifQ          = al2.Qualities.substr( 0,umif);
			
			al2.QueryBases = al2.QueryBases.substr(umif);
			al2.Qualities  = al2.Qualities.substr( umif);
			
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
			
			al.QueryBases  = al.QueryBases.substr(umir);
			al.Qualities   = al.Qualities.substr( umir) ;

			
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
		    chunkToWrite->dataToProcess->push_back(al2);
		    chunkToWrite->dataToProcess->push_back(al);

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
		mtr->processSingle(al);
		unsigned int newL=al.QueryBases.length();
				
		if(keepOrig){
		    //write duplicate
		    if(orgL  != newL){
			orig.SetIsDuplicate(true);
			chunkToWrite->dataToProcess->push_back(orig);
		    }
		}
		//if(!al.AddTag("RX","Z",umifS)){ cerr << "Error while editing tags for UMI:"<<umifS<<"#"<< endl;  exit(1); }

		
		if(umif != 0){	    
		    if(al.QueryBases.length() < (umif+umir) ){//the new fragment is shorter than both umis, cannot trim 
			al.SetIsDuplicate(true);  //leave but fail QC
		    }else{				    
			umifS         = al.QueryBases.substr(0,umif);
			umifQ         = al.Qualities.substr(0, umif);

			al.QueryBases = al.QueryBases.substr(umif);
			al.Qualities  = al.Qualities.substr( umif);

			if(umir == 0){	//no reverse UMI, add the tags now
			    if(!al.AddTag("RX","Z",umifS)){  cerr << "Error while editing RX tags for UMI:"<<umifS<<"#"<< endl; exit(1);   }				    
			    if(!al.AddTag("QX","Z",umifQ)){  cerr << "Error while editing QX tags for UMI:"<<umifQ<<"#"<< endl; exit(1);   }
			}	
		    }
		}

		if(umir != 0){	    
		    if(al.QueryBases.length() < (umif+umir) ){//the new fragment is shorter than both umis, cannot trim 
			al.SetIsDuplicate(true);  //leave but fail QC
		    }else{				    
			if(orgL  != newL){//was trimmed
			    umirS         = al.QueryBases.substr(al.QueryBases.length()-umir,  umir);
			    umirQ         = al.Qualities.substr( al.Qualities.length() -umir,  umir);
			
			    al.QueryBases = al.QueryBases.substr(0,al.QueryBases.length()-umir);
			    al.Qualities  = al.Qualities.substr( 0,al.Qualities.length() -umir);
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
		chunkToWrite->dataToProcess->push_back(al);



	    } //end single end
	    al2Null=true;
	}//second pair
		    

    }


	


    ////cerr<<"Thread #"<<rankThread <<" is done with computations"<<endl;

    //////////////////////////////////////////////////////////////
    //                END   COMPUTATION                         //
    //////////////////////////////////////////////////////////////
#ifdef DEBUG    
    int rankLastChunk = currentChunk->rank;
#endif

    delete currentChunk;//input chunk
    //TO REMOVE
    //int fakeSleep=int( 10*(1/double(rankThread)));
    //COUNTERS

#ifdef DEBUG    
    rc = pthread_mutex_lock(&mutexCERR);
    checkResults("pthread_mutex_lock()\n", rc);
	       
    //cerr<<"Thread #"<<rankThread <<" is done with computations, sleep "<<fakeSleep<<" rank produced "<<rankLastChunk<<endl;
    cerr<<"Thread #"<<rankThread <<" is done with computations, rank produced "<<rankLastChunk<<endl;
    //sleep(fakeSleep);
    rc = pthread_mutex_unlock(&mutexCERR);
    checkResults("pthread_mutex_unlock()\n", rc);     
#endif

    //take mutex    
    rc = pthread_mutex_lock(&mutexCounter);
    checkResults("pthread_mutex_lock()\n", rc);


    //---> mutex is taken here <---

    //if(int(queueDataTowrite.size())>(2*numberOfThreads)){	
    //if queue is full
    if(int(queueDataTowrite.size())>(maxQueueDataTowriteSize)){		

	forceInsertProcedure(chunkToWrite,rankThread);

    }else{// just insert it if there is space


#ifdef DEBUG
	rc = pthread_mutex_lock(&mutexCERR);
	checkResults("pthread_mutex_lock()\n", rc);
	       
	cerr<<"Thread #"<<rankThread<<" has space in queue, queue size= "<<int(queueDataTowrite.size())<<" rank of chunk="<<chunkToWrite->rank<<endl;

	rc = pthread_mutex_unlock(&mutexCERR);
	checkResults("pthread_mutex_unlock()\n", rc);     
#endif

	queueDataTowrite.push(chunkToWrite);
	//release mutex
	rc = pthread_mutex_unlock(&mutexCounter);
	checkResults("pthread_mutex_unlock()\n", rc);	


    }
    
#ifdef DEBUG
    rc = pthread_mutex_lock(&mutexCERR);
    checkResults("pthread_mutex_lock()\n", rc);
    
    cerr<<"Thread #"<<rankThread <<" is re-starting"<<endl;

    rc = pthread_mutex_unlock(&mutexCERR);
    checkResults("pthread_mutex_unlock()\n", rc);     
#endif

    goto checkqueue;	   


    

#ifdef DEBUG     
	rc = pthread_mutex_lock(&mutexCERR);
	checkResults("pthread_mutex_lock()\n", rc);
	       
	cerr<<"Thread "<<rankThread<<" ended "<<endl;

	rc = pthread_mutex_unlock(&mutexCERR);
	checkResults("pthread_mutex_unlock()\n", rc);
#endif

    return NULL;

}





//! Method called for each thread FQ
/*!
  

*/				
void *mainComputationThreadFQ(void * argc){    
    // bool * singleEndModeFQAddr = (bool *)argc;
    // bool   singleEndModeFQ     = *singleEndModeFQAddr;

    bool   singleEndModeFQ = ((struct argspthread*)argc)->singleEndModeFQ;
    unsigned int    umif   = ((struct argspthread*)argc)->umif;
    unsigned int    umir   = ((struct argspthread*)argc)->umir;
    bool   umi             = ((umif+umir)>0);
    
    int   rc;
    int rankThread=0;

    rc = pthread_mutex_lock(&mutexRank);
    checkResults("pthread_mutex_lock()\n", rc);

    threadID2Rank[*(int *)pthread_self()]  = threadID2Rank.size()+1;
    rankThread = threadID2Rank[*(int *)pthread_self()];

    
    rc = pthread_mutex_unlock(&mutexRank);
    checkResults("pthread_mutex_unlock()\n", rc);

 checkqueue:    
    // stackIndex=-1;
    //check stack

    
    rc = pthread_mutex_lock(&mutexQueue);
    checkResults("pthread_mutex_lock()\n", rc);


    bool foundData=false;
    
    //cerr<<"Thread #"<<rankThread <<" started and is requesting data"<<endl;

    DataChunkFQ * currentChunk;


    if(!queueDataToprocessFQ.empty()){    
 	foundData=true;
 	currentChunk = queueDataToprocessFQ.front();
 	queueDataToprocessFQ.pop();

#ifdef DEBUG 
	rc = pthread_mutex_lock(&mutexCERR);
	checkResults("pthread_mutex_lock()\n", rc);
		
	cerr<<"Thread #"<<rankThread<<" is reading "<<currentChunk->rank<<" with size "<<currentChunk->dataToProcess->size()<<endl;
	
	rc = pthread_mutex_unlock(&mutexCERR);
	checkResults("pthread_mutex_unlock()\n", rc);
#endif
	//cout<<"rank "<< &(currentChunk->dataToProcess) <<endl;
    }

    
  

    if(!foundData){
	rc = pthread_mutex_unlock(&mutexQueue);
	checkResults("pthread_mutex_unlock()\n", rc);


	if(readDataDone){

#ifdef DEBUG 
	    rc = pthread_mutex_lock(&mutexCERR);
	    checkResults("pthread_mutex_lock()\n", rc);
	    
	    cerr<<"Thread #"<<rankThread<<" is done"<<endl;
	
	    rc = pthread_mutex_unlock(&mutexCERR);
	    checkResults("pthread_mutex_unlock()\n", rc);	
#endif

	    return NULL;
	}else{

#ifdef DEBUGSLEEPING
	    rc = pthread_mutex_lock(&mutexCERR);
	    checkResults("pthread_mutex_lock()\n", rc);	    
	    cerr<<"Thread #"<<rankThread<<" has not found data, sleeping for "<<timeThreadSleep<<"s"<<endl;	    
	    rc = pthread_mutex_unlock(&mutexCERR);
	    checkResults("pthread_mutex_unlock()\n", rc);
#endif

	    sleep(timeThreadSleep);
	    goto checkqueue;
	}

    }else{
	//release stack
	rc = pthread_mutex_unlock(&mutexQueue);
	checkResults("pthread_mutex_unlock()\n", rc);
    }

    //////////////////////////////////////////////////////////////
    //                BEGIN COMPUTATION                         //
    //////////////////////////////////////////////////////////////
    DataChunkFQ * chunkToWrite = new DataChunkFQ(currentChunk->rank);
    string umifS;
    string umirS;
    string umifQ;
    string umirQ;
    
    string cmt;
    // if(dataToWrite->dataToProcess->at(i).sequence != ""){ //new sequence			    
    
    
    for(unsigned i=0;i<currentChunk->dataToProcess->size();i++){
	//cerr<<"process "<<currentChunk->dataToProcess->at(i).d1<<endl;
	fqrecord toadd;


	
	if(!singleEndModeFQ){
	    

	    merged result=	mtr->process_PE(currentChunk->dataToProcess->at(i).s1,currentChunk->dataToProcess->at(i).q1,
						currentChunk->dataToProcess->at(i).s2,currentChunk->dataToProcess->at(i).q2);
	    
	    

	    if(umi){
		cmt.clear();
	        if( result.sequence != ""){//has new sequence

		    if( result.sequence.length() < (umif+umir) ){//the resulting read is too short
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

			    result.sequence = result.sequence.substr(umif);
			    result.quality  = result.quality.substr( umif);
			    
			    if(umir == 0){//add them now
				cmt = cmt+"\tRX:Z:"+umifS;
				cmt = cmt+"\tQX:Z:"+umifQ;
			    }			    
			}
			
			//if UMI reverse
			if(umir != 0){
			    umirS           = result.sequence.substr(  result.sequence.length()-umir,  umir);
			    umirQ           = result.quality.substr(   result.quality.length() -umir,  umir);
			    
			    result.sequence = result.sequence.substr(0,result.sequence.length()-umir);
			    result.quality  = result.quality.substr(0, result.quality.length() -umir);
			    
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
			umifS                                 = currentChunk->dataToProcess->at(i).s1.substr(0,  umif);
			umifQ                                 = currentChunk->dataToProcess->at(i).q1.substr(0,  umif);
			
			currentChunk->dataToProcess->at(i).s1 = currentChunk->dataToProcess->at(i).s1.substr(umif);
			currentChunk->dataToProcess->at(i).q1 = currentChunk->dataToProcess->at(i).q1.substr(umif);	 	

			if(umir == 0){//add them now
			    cmt = cmt+"\tRX:Z:"+umifS;
			    cmt = cmt+"\tQX:Z:"+umifQ;
			}	
		    }
		    

		    if(umir != 0){
			umirS                                 = reverseComplement(currentChunk->dataToProcess->at(i).s2.substr(0,  umir));
			umirQ                                 =                   currentChunk->dataToProcess->at(i).q2.substr(0,  umir);
			reverse(umirQ.begin(),umirQ.end());

			currentChunk->dataToProcess->at(i).s2 = currentChunk->dataToProcess->at(i).s2.substr(0,currentChunk->dataToProcess->at(i).s2.length()-umir);
			currentChunk->dataToProcess->at(i).q2 = currentChunk->dataToProcess->at(i).q2.substr(0,currentChunk->dataToProcess->at(i).q2.length()-umir);
			
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
	    
	    //mtr->incrementCountall();
	    toadd.code         = result.code;
	    
	    toadd.sequence     = result.sequence;
	    toadd.quality      = result.quality;
	    toadd.cmt          = cmt;

	    toadd.d1           = currentChunk->dataToProcess->at(i).d1;
	    toadd.s1           = currentChunk->dataToProcess->at(i).s1;
	    toadd.q1           = currentChunk->dataToProcess->at(i).q1;
	    
	    toadd.d2           = currentChunk->dataToProcess->at(i).d2;
	    toadd.s2           = currentChunk->dataToProcess->at(i).s2;
	    toadd.q2           = currentChunk->dataToProcess->at(i).q2;

	    toadd.paired       = true;

	    // continue;



	}else{//V single mode below V
		
		
	    merged result=mtr->process_SR(currentChunk->dataToProcess->at(i).s1,
					  currentChunk->dataToProcess->at(i).q1); //*(fo1->getSeq()),*(fo1->getQual()));


	    if(umi){
		cmt.clear();
	        if( result.sequence != ""){//has new sequence

		    if( result.sequence.length() < (umif+umir) ){//the resulting read is too short
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

			    result.sequence = result.sequence.substr(umif);
			    result.quality  = result.quality.substr( umif);
			    
			    if(umir == 0){//add them now
				cmt = cmt+"\tRX:Z:"+umifS;
				cmt = cmt+"\tQX:Z:"+umifQ;
			    }			    
			}
			
			//if UMI reverse
			if(umir != 0){
			    umirS           = result.sequence.substr(  result.sequence.length()-umir,  umir);
			    umirQ           = result.quality.substr(   result.quality.length() -umir,  umir);
			    
			    result.sequence = result.sequence.substr(0,result.sequence.length()-umir);
			    result.quality  = result.quality.substr(0, result.quality.length() -umir);
			    
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
			umifS                                 = currentChunk->dataToProcess->at(i).s1.substr(0,  umif);
			umifQ                                 = currentChunk->dataToProcess->at(i).q1.substr(0,  umif);
			
			currentChunk->dataToProcess->at(i).s1 = currentChunk->dataToProcess->at(i).s1.substr(umif);
			currentChunk->dataToProcess->at(i).q1 = currentChunk->dataToProcess->at(i).q1.substr(umif);	 	

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

	   

	    mtr->incrementCountall();
	    toadd.code = result.code;
	    toadd.paired  =  false;


	    
	    toadd.sequence     = result.sequence;
	    toadd.quality      = result.quality;
	    toadd.cmt          = cmt;
	    toadd.d1           = currentChunk->dataToProcess->at(i).d1;
	    toadd.s1           = currentChunk->dataToProcess->at(i).s1;
	    toadd.q1           = currentChunk->dataToProcess->at(i).q1;
	    


	}
	


	//totalSeqs++;
    
	chunkToWrite->dataToProcess->push_back(toadd);


    }


	


    ////cerr<<"Thread #"<<rankThread <<" is done with computations"<<endl;

    //////////////////////////////////////////////////////////////
    //                END   COMPUTATION                         //
    //////////////////////////////////////////////////////////////

#ifdef DEBUG    
    int rankLastChunk = currentChunk->rank;
#endif

    delete currentChunk;//input chunk
    //TO REMOVE
    //int fakeSleep=int( 10*(1/double(rankThread)));
    //COUNTERS

#ifdef DEBUG    
    rc = pthread_mutex_lock(&mutexCERR);
    checkResults("pthread_mutex_lock()\n", rc);
    
    //cerr<<"Thread #"<<rankThread <<" is done with computations, sleep "<<fakeSleep<<" rank produced "<<rankLastChunk<<endl;
    cerr<<"Thread #"<<rankThread <<" is done with computations, rank produced "<<rankLastChunk<<" size of data "<<chunkToWrite->dataToProcess->size()<<endl;
    //sleep(fakeSleep);
    rc = pthread_mutex_unlock(&mutexCERR);
    checkResults("pthread_mutex_unlock()\n", rc);     
#endif

    //take mutex    
    rc = pthread_mutex_lock(&mutexCounter);
    checkResults("pthread_mutex_lock()\n", rc);


    //---> mutex is taken here <---

    //if(int(queueDataTowrite.size())>(2*numberOfThreads)){	
    //if queue is full
    if(int(queueDataTowriteFQ.size())>(maxQueueDataTowriteSize)){		

	forceInsertProcedureFQ(chunkToWrite,rankThread);

    }else{// just insert it if there is space


#ifdef DEBUG
	rc = pthread_mutex_lock(&mutexCERR);
	checkResults("pthread_mutex_lock()\n", rc);
	       
	cerr<<"Thread #"<<rankThread<<" has space in queue, queue size= "<<int(queueDataTowriteFQ.size())<<" rank of chunk="<<chunkToWrite->rank<<endl;

	rc = pthread_mutex_unlock(&mutexCERR);
	checkResults("pthread_mutex_unlock()\n", rc);     
#endif

	queueDataTowriteFQ.push(chunkToWrite);
	//release mutex
	rc = pthread_mutex_unlock(&mutexCounter);
	checkResults("pthread_mutex_unlock()\n", rc);	


    }
    
#ifdef DEBUG
    rc = pthread_mutex_lock(&mutexCERR);
    checkResults("pthread_mutex_lock()\n", rc);
    
    cerr<<"Thread #"<<rankThread <<" is re-starting"<<endl;

    rc = pthread_mutex_unlock(&mutexCERR);
    checkResults("pthread_mutex_unlock()\n", rc);     
#endif

    goto checkqueue;	   


    

#ifdef DEBUG     

	rc = pthread_mutex_lock(&mutexCERR);
	checkResults("pthread_mutex_lock()\n", rc);
	       
	cerr<<"Thread "<<rankThread<<" ended "<<endl;

	rc = pthread_mutex_unlock(&mutexCERR);
	checkResults("pthread_mutex_unlock()\n", rc);

#endif

    return NULL;

}



bool checkForWritingLoop(BamWriter * writer,int * lastWrittenChunk){


    //check something to write
    int rc = pthread_mutex_lock(&mutexCounter);
    checkResults("pthread_mutex_lock()\n", rc);
	    
    if(!queueDataTowrite.empty()){
	DataChunk *  dataToWrite= queueDataTowrite.top();
	if( *lastWrittenChunk == (dataToWrite->rank-1) ){

#if defined(DEBUG) ||defined(DEBUGSLEEPING)
	    int qsize = queueDataTowrite.size(); 
#endif

	    queueDataTowrite.pop();
	    rc = pthread_mutex_unlock(&mutexCounter);
	    checkResults("pthread_mutex_unlock()\n", rc);

#ifdef DEBUG 
	    rc = pthread_mutex_lock(&mutexCERR);
	    checkResults("pthread_mutex_lock()\n", rc);			
	    cerr<<"loop: writing "<<dataToWrite->rank<<" queue size= "<<qsize<<endl;			
	    rc = pthread_mutex_unlock(&mutexCERR);
	    checkResults("pthread_mutex_unlock()\n", rc);
#endif

	    //writing dataToWrite
	    for(unsigned int i=0;i<dataToWrite->dataToProcess->size();i++){
		//cout<<"loop: writing "<<dataToWrite->dataToProcess->at(i).Name<<endl;
		writer->SaveAlignment( dataToWrite->dataToProcess->at(i) );
		// outfile<< dataToWrite->dataToProcess->at(i).Name 
		//        << endl  
		//        << dataToWrite->dataToProcess->at(i).QueryBases << endl
		//        << "+" <<endl 
		//        << dataToWrite->dataToProcess->at(i).Qualities << endl;
	    }
	    *lastWrittenChunk=dataToWrite->rank;

#ifdef DEBUG 
	    rc = pthread_mutex_lock(&mutexCERR);
	    checkResults("pthread_mutex_lock()\n", rc);			
	    cerr<<"loop: written "<<dataToWrite->rank<<" lastWrittenChunk= "<<*lastWrittenChunk<<endl;			
	    rc = pthread_mutex_unlock(&mutexCERR);
	    checkResults("pthread_mutex_unlock()\n", rc);
#endif

	    delete dataToWrite;
	    return true;
	}else{
			
#ifdef DEBUG 
	    rc = pthread_mutex_lock(&mutexCERR);
	    checkResults("pthread_mutex_lock()\n", rc);
			
	    cerr<<"loop: is waiting for chunk "<<*lastWrittenChunk<<" but got  "<<(dataToWrite->rank)<<endl;
			
	    rc = pthread_mutex_unlock(&mutexCERR);
	    checkResults("pthread_mutex_unlock()\n", rc);
#endif


	    rc = pthread_mutex_unlock(&mutexCounter);
	    checkResults("pthread_mutex_unlock()\n", rc);

	    return false;
	}
    }else{//if queue is empty
	//cout<<"queue to write is empty "<<endl;

		    			
#ifdef DEBUG 
	rc = pthread_mutex_lock(&mutexCERR);
	checkResults("pthread_mutex_lock()\n", rc);
		    
	cerr<<"loop: writing queue is empty"<<endl;
			
	rc = pthread_mutex_unlock(&mutexCERR);
	checkResults("pthread_mutex_unlock()\n", rc);
#endif

	rc = pthread_mutex_unlock(&mutexCounter);
	checkResults("pthread_mutex_unlock()\n", rc);

	return false;
    }
    return false;

}



bool checkForWritingLoopFQ(fqwriters * onereadgroup,int * lastWrittenChunk,MergeTrimReads *     mtr,    unsigned int umif,    unsigned int umir){


    //check something to write
    int rc = pthread_mutex_lock(&mutexCounter);
    checkResults("pthread_mutex_lock()\n", rc);
	    
    if(!queueDataTowriteFQ.empty()){
	DataChunkFQ *  dataToWrite = queueDataTowriteFQ.top();
	if( *lastWrittenChunk == (dataToWrite->rank-1) ){

#if defined(DEBUG) ||defined(DEBUGSLEEPING)
	    int qsize = queueDataTowriteFQ.size(); 
#endif

	    queueDataTowriteFQ.pop();
	    rc = pthread_mutex_unlock(&mutexCounter);
	    checkResults("pthread_mutex_unlock()\n", rc);

#ifdef DEBUG 
	    rc = pthread_mutex_lock(&mutexCERR);
	    checkResults("pthread_mutex_lock()\n", rc);			
	    cerr<<"loop: writing "<<dataToWrite->rank<<" queue size= "<<qsize<<endl;			
	    rc = pthread_mutex_unlock(&mutexCERR);
	    checkResults("pthread_mutex_unlock()\n", rc);
#endif

	    //writing dataToWrite
	    for(unsigned int i=0;i<dataToWrite->dataToProcess->size();i++){
		if(dataToWrite->dataToProcess->at(i).paired){//PAIRED-END

		    mtr->incrementCountall();
		    
		    if(dataToWrite->dataToProcess->at(i).code != ' '){ //keys or chimeras
			
			if(dataToWrite->dataToProcess->at(i).code == 'K'){
			    mtr->incrementCountfkey();
			}else{
			    if(dataToWrite->dataToProcess->at(i).code == 'D'){
				mtr->incrementCountchimera();
			    }else{
				if(dataToWrite->dataToProcess->at(i).code == 'U'){//UMI problem, fragment is too short, do not increment anything
				    mtr->incrementCountUMIp();
				}else{
				    cerr << "leehom: Wrong return code =\""<<dataToWrite->dataToProcess->at(i).code<<"\""<<endl;
				    exit(1);
				}
			    }
			}
			
			onereadgroup->pairr2f<<"@"<<dataToWrite->dataToProcess->at(i).d2<<"/2" <<dataToWrite->dataToProcess->at(i).cmt<<endl <<dataToWrite->dataToProcess->at(i).s2<<endl<<"+"<<endl <<dataToWrite->dataToProcess->at(i).q2<<endl;
			onereadgroup->pairr1f<<"@"<<dataToWrite->dataToProcess->at(i).d1<<"/1" <<dataToWrite->dataToProcess->at(i).cmt<<endl <<dataToWrite->dataToProcess->at(i).s1<<endl<<"+"<<endl <<dataToWrite->dataToProcess->at(i).q1<<endl;

			continue;
		    
		    }else{
		        if(dataToWrite->dataToProcess->at(i).sequence != ""){ //new sequence			    
			    onereadgroup->single<<"@"<<dataToWrite->dataToProcess->at(i).d1<<" " <<dataToWrite->dataToProcess->at(i).cmt<<endl << dataToWrite->dataToProcess->at(i).sequence<<endl<<"+"<<endl <<dataToWrite->dataToProcess->at(i).quality<<endl;    	    
			    
			    if( (dataToWrite->dataToProcess->at(i).sequence.length()+umif+umir) > max(dataToWrite->dataToProcess->at(i).s1.length(),dataToWrite->dataToProcess->at(i).s2.length()) ){
				mtr->incrementCountmergedoverlap();
			    }else{
				mtr->incrementCountmerged();			  
			    }

			}else{ //keep as is
			    mtr->incrementCountnothing();
			    
			    onereadgroup->pairr2<<"@"<<dataToWrite->dataToProcess->at(i).d2<<"/2" <<dataToWrite->dataToProcess->at(i).cmt<<  endl <<dataToWrite->dataToProcess->at(i).s2<<endl<<"+"<<endl <<dataToWrite->dataToProcess->at(i).q2<<endl;
			    onereadgroup->pairr1<<"@"<<dataToWrite->dataToProcess->at(i).d1<<"/1" <<dataToWrite->dataToProcess->at(i).cmt<<  endl <<dataToWrite->dataToProcess->at(i).s1<<endl<<"+"<<endl <<dataToWrite->dataToProcess->at(i).q1<<endl;
			    
			}
		    }
		    
		}else{//SINGLE-END

		    
		    if(dataToWrite->dataToProcess->at(i).code != ' '){ //either chimera or missing key
			
			if(dataToWrite->dataToProcess->at(i).code == 'K'){
			    mtr->incrementCountfkey();
			}else{
			    if(dataToWrite->dataToProcess->at(i).code == 'D'){
				mtr->incrementCountchimera();
			    }else{
				if(dataToWrite->dataToProcess->at(i).code == 'U'){//UMI problem
				    mtr->incrementCountUMIp();
				}else{
				    cerr << "leehom: Wrong return code =\""<<dataToWrite->dataToProcess->at(i).code<<"\""<<endl;
				    exit(1);
				}
			    }
			}
			
			onereadgroup->singlef<<"@"<<dataToWrite->dataToProcess->at(i).d1<<" " <<dataToWrite->dataToProcess->at(i).cmt<<endl << dataToWrite->dataToProcess->at(i).s1 <<endl<<"+"<<endl << dataToWrite->dataToProcess->at(i).q1 <<endl;
			continue;
		    }
		    
		    if(dataToWrite->dataToProcess->at(i).sequence != ""){ //new sequence
			mtr->incrementCounttrimmed();
			onereadgroup->single<<"@"<<dataToWrite->dataToProcess->at(i).d1<<" "<<dataToWrite->dataToProcess->at(i).cmt <<endl << dataToWrite->dataToProcess->at(i).sequence <<endl<<"+"<<endl << dataToWrite->dataToProcess->at(i).quality<<endl;
		    }else{
			mtr->incrementCountnothing();
			onereadgroup->single<<"@"<<dataToWrite->dataToProcess->at(i).d1<<" "<<dataToWrite->dataToProcess->at(i).cmt <<endl << dataToWrite->dataToProcess->at(i).s1       <<endl<<"+"<<endl << dataToWrite->dataToProcess->at(i).q1     <<endl;
		    }
		}
			
	    }
	    *lastWrittenChunk=dataToWrite->rank;

#ifdef DEBUG 
	    rc = pthread_mutex_lock(&mutexCERR);
	    checkResults("pthread_mutex_lock()\n", rc);			
	    cerr<<"loop: written "<<dataToWrite->rank<<" lastWrittenChunk= "<<*lastWrittenChunk<<endl;			
	    rc = pthread_mutex_unlock(&mutexCERR);
	    checkResults("pthread_mutex_unlock()\n", rc);
#endif

	    delete dataToWrite;
	    return true;
	}else{
			
#ifdef DEBUG 
	    rc = pthread_mutex_lock(&mutexCERR);
	    checkResults("pthread_mutex_lock()\n", rc);
			
	    cerr<<"loop: is waiting for chunk "<<*lastWrittenChunk<<" but got  "<<(dataToWrite->rank)<<endl;
			
	    rc = pthread_mutex_unlock(&mutexCERR);
	    checkResults("pthread_mutex_unlock()\n", rc);
#endif


	    rc = pthread_mutex_unlock(&mutexCounter);
	    checkResults("pthread_mutex_unlock()\n", rc);
	    return false;
	}
    }else{//if queue is empty
	//cout<<"queue to write is empty "<<endl;

		    			
#ifdef DEBUG 
	rc = pthread_mutex_lock(&mutexCERR);
	checkResults("pthread_mutex_lock()\n", rc);
		    
	cerr<<"loopfq: writing queue is empty"<<endl;
			
	rc = pthread_mutex_unlock(&mutexCERR);
	checkResults("pthread_mutex_unlock()\n", rc);
#endif

	rc = pthread_mutex_unlock(&mutexCounter);
	checkResults("pthread_mutex_unlock()\n", rc);
	return false;
    }
    return false;
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
    bool   trimKey=false;
    
    bool   allowMissing=false;
    int    trimCutoff=1;

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

    unsigned int umif=0;
    unsigned int umir=0;
    
    //int    numberOfThreads   = 1;

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
			      "\t\t"+"-f , --adapterFirstRead" +"\t\t\t"+"Adapter that is observed after the forward read (def. Multiplex: "+options_adapter_F_BAM.substr(0,30)+")"+"\n"+
			      "\t\t"+"-s , --adapterSecondRead" +"\t\t"+"Adapter that is observed after the reverse read (def. Multiplex: "+options_adapter_S_BAM.substr(0,30)+")"+"\n"+
			      "\t\t"+"-c , --FirstReadChimeraFilter" +"\t\t"+"If the forward read looks like this sequence, the cluster is filtered out.\n\t\t\t\t\t\t\tProvide several sequences separated by comma (def. Multiplex: "+options_adapter_chimera_BAM.substr(0,30)+")"+"\n"+
			      "\t\t"+"-k , --key"+"\t\t\t\t"+"Key sequence with which each sequence starts. Comma separate for forward and reverse reads. (default '"+key+"')"+"\n"+
			      "\t\t"+"--trimkey"     +"\t\t\t\t"+"Trim the key sequence even for untrimmed. (default '"+boolStringify(trimKey)+"')"+"\n"+

			      "\t\t"+"-i , --allowMissing"+"\t\t\t"+"Allow one base in one key to be missing or wrong. (default "+boolStringify(allowMissing)+")"+"\n"+

			      "\t\t"+"--umif [bp]"+"\t\t\t"+"Extract bp for the UMI for the forward read. (default "+stringify(umif)+")"+"\n"+
			      "\t\t"+"--umir [bp]"+"\t\t\t"+"Extract bp for the UMI for the reverse read. (default "+stringify(umir)+")"+"\n"+

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
	//int    numberOfThreads   = 1;

	if(strcmp(argv[i],"-t") == 0 ){
	    numberOfThreads =destringify<int>(argv[i+1]);
	    sizeChunk       =5000*numberOfThreads;
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

	
	if(strcmp(argv[i],"--trimCutoff") == 0 ){
	    trimCutoff=atoi(argv[i+1]);
	    i++;
	    continue;
	}
	
	cerr<<"Unknown option "<<argv[i] <<" exiting"<<endl;
	return 1;	    
    }

    maxQueueDataTowriteSize   = 3*numberOfThreads;
    maxQueueDataToprocessSize = 3*numberOfThreads;


#ifdef DEBUG 
    cerr<<"sizeChunk "<<sizeChunk<<" queue size  "<<maxQueueDataTowriteSize<<endl;
#endif

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

    mtr = new MergeTrimReads(adapter_F,adapter_S,adapter_chimera,
			     key1,key2,trimKey,
			     trimCutoff,allowMissing,ancientDNA,location,scale,useDist,qualOffset);
    
    fqwriters onereadgroup;
    if(fastqFormat){
	//TODO

	
	pthread_mutex_init(&mutexCERR,    NULL);
	pthread_mutex_init(&mutexQueue,   NULL);
	pthread_mutex_init(&mutexCounter, NULL);
	pthread_mutex_init(&mutexRank ,   NULL);
    
	pthread_t             thread[numberOfThreads];
	int                   rc=0;


	 
	struct argspthread *argstopass = (struct argspthread *)malloc(sizeof(struct argspthread));
	argstopass->singleEndModeFQ = singleEndModeFQ;
	argstopass->umif            = umif;
	argstopass->umir            = umir;
	    
	for(int i=0;i<numberOfThreads;i++){
	    //rc = pthread_create(&thread[i], NULL, mainComputationThreadFQ, &singleEndModeFQ);
	    rc = pthread_create(&thread[i], NULL, mainComputationThreadFQ, (void *)argstopass );
	    checkResults("pthread_create()\n", rc);
	}

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

	    if(!onereadgroup.single.good()){       cerr<<"Cannot write to file "<<outdirs<<endl;  return 1; }
	    if(!onereadgroup.pairr1.good()){       cerr<<"Cannot write to file "<<outdir1<<endl;  return 1; }
	    if(!onereadgroup.pairr2.good()){       cerr<<"Cannot write to file "<<outdir2<<endl;  return 1; }
	    
	    if(!onereadgroup.singlef.good()){      cerr<<"Cannot write to file "<<outdirsf<<endl; return 1; }
	    if(!onereadgroup.pairr1f.good()){      cerr<<"Cannot write to file "<<outdir1f<<endl; return 1; }
	    if(!onereadgroup.pairr2f.good()){      cerr<<"Cannot write to file "<<outdir2f<<endl; return 1; }	    
	}

	unsigned int rank=0;

	DataChunkFQ * currentChunkfq = new DataChunkFQ(rank);
	unsigned int counter=0;
	int lastWrittenChunk=-1;

	//unsigned int totalSeqs=0;
	while(fqp1->hasData()){

	    FastQObj * fo1=fqp1->getData();
	    vector<string> def1=allTokens( *(fo1->getID()), ' '  );
	    string def1s=def1[0];
	

	    FastQObj * fo2;
	    string def2s;
	    string ext2s;
	    
	    fqrecord fqtoadd;
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
	

	    if(counter== (sizeChunk-1)){
		     
		currentChunkfq->dataToProcess->push_back(fqtoadd);
	       
		rc = pthread_mutex_lock(&mutexQueue);
		checkResults("pthread_mutex_lock()\n", rc);
		//cout<<"got it"<<endl;
		//cout<<"main addr"<< &(currentChunk->dataToProcess) <<endl;
		//
		while(true){
		    //mutex queue is taken here
		    if(int(queueDataToprocessFQ.size())>(maxQueueDataToprocessSize)){	  // if queue full
#ifdef DEBUGSLEEPING
			int queueSize = int(queueDataToprocessFQ.size());
#endif
			
			//release mutex queue
			rc = pthread_mutex_unlock(&mutexQueue);
			checkResults("pthread_mutex_unlock()\n", rc);

			bool hasWritten = checkForWritingLoopFQ(&onereadgroup,&lastWrittenChunk,mtr,umif,umir);
			if(!hasWritten){
#ifdef DEBUGSLEEPING
			    rc = pthread_mutex_lock(&mutexCERR);
			    checkResults("pthread_mutex_lock()\n", rc);			
			    cerr<<"main loop: queue to process is full, capacity: "<<queueSize<<" sleeping for "<<timeThreadSleep<<"s"<<endl;			
			    rc = pthread_mutex_unlock(&mutexCERR);
			    checkResults("pthread_mutex_unlock()\n", rc);
#endif

			    sleep(timeThreadSleep);
			}
			

			//take queue mutex
			rc = pthread_mutex_lock(&mutexQueue);
			checkResults("pthread_mutex_lock()\n", rc);
		    }else{
#ifdef DEBUG 
			rc = pthread_mutex_lock(&mutexCERR);
			checkResults("pthread_mutex_lock()\n", rc);			
			cerr<<"main loop: queue to process has space, adding rank: "<<currentChunkfq->rank<<endl;			
			rc = pthread_mutex_unlock(&mutexCERR);
			checkResults("pthread_mutex_unlock()\n", rc);
#endif

			queueDataToprocessFQ.push(currentChunkfq);

			rc = pthread_mutex_unlock(&mutexQueue);
			checkResults("pthread_mutex_unlock()\n", rc);
			break;
		    }
		}



		rc = pthread_mutex_unlock(&mutexQueue);
		checkResults("pthread_mutex_unlock()\n", rc);

		rank++;

		checkForWritingLoopFQ(&onereadgroup,&lastWrittenChunk,mtr,umif,umir);
	 
		currentChunkfq = new DataChunkFQ(rank);	    
		//currentChunk->rank=rank;
		counter=0;
	    }else{
		currentChunkfq->dataToProcess->push_back(fqtoadd);

		counter++;
	    }



	    //totalSeqs++;
	}//for each data record
    
	delete fqp1;
	if(!singleEndModeFQ){
	    delete fqp2;
	}
		    			
#ifdef DEBUG 
	rc = pthread_mutex_lock(&mutexCERR);
	checkResults("pthread_mutex_lock()\n", rc);
	
	cerr<<"main loop: is done"<<endl;
	
	rc = pthread_mutex_unlock(&mutexCERR);
	checkResults("pthread_mutex_unlock()\n", rc);
#endif




	//add last chunk
	int lastRank=currentChunkfq->rank;
	queueDataToprocessFQ.push(currentChunkfq);
	bool wroteEverything=false;

	readDataDone=true;
	//cout<<"main thread is done reading"<<endl;

	while(!wroteEverything){
	    //cout<<"write queue size = "<< queueDataTowrite.size() <<endl;
	    //threads are running here
	    rc = pthread_mutex_lock(&mutexCounter);
	    checkResults("pthread_mutex_lock()\n", rc);

	    //#ifdef DEBUG    
#if defined(DEBUG) ||defined(DEBUGSLEEPING)
	    int qsize = queueDataTowriteFQ.size(); 
#endif
	    //#endif

	    bool wroteData=false;
	    if(!queueDataTowriteFQ.empty()){
	
		DataChunkFQ *  dataToWrite= queueDataTowriteFQ.top();
		
		if( lastWrittenChunk == (dataToWrite->rank-1) ){
		    queueDataTowriteFQ.pop();
		    
		    rc = pthread_mutex_unlock(&mutexCounter);
		    checkResults("pthread_mutex_unlock()\n", rc);

#ifdef DEBUG 
		    rc = pthread_mutex_lock(&mutexCERR);
		    checkResults("pthread_mutex_lock()\n", rc);	       
		    cerr<<"final: writing "<<dataToWrite->rank<<" queue size= "<<qsize<<" size of data to write "<<dataToWrite->dataToProcess->size()<<endl;
		    rc = pthread_mutex_unlock(&mutexCERR);
		    checkResults("pthread_mutex_unlock()\n", rc);
#endif

		    for(unsigned int i=0;i<dataToWrite->dataToProcess->size();i++){
			// cerr<<dataToWrite->dataToProcess->at(i).d1<<endl;
			if(dataToWrite->dataToProcess->at(i).paired){
			    
			    mtr->incrementCountall();
		    
			    if(dataToWrite->dataToProcess->at(i).code != ' '){ //keys or chimeras
			
				if(dataToWrite->dataToProcess->at(i).code == 'K'){
				    mtr->incrementCountfkey();
				}else{
				    if(dataToWrite->dataToProcess->at(i).code == 'D'){
					mtr->incrementCountchimera();
				    }else{
					if(dataToWrite->dataToProcess->at(i).code == 'U'){//UMI problem, sequence is too short
					    mtr->incrementCountUMIp();
					}else{					    
					    cerr << "leehom: Wrong return code =\""<<dataToWrite->dataToProcess->at(i).code<<"\""<<endl;
					    exit(1);
					}
				    }
				}
				
				onereadgroup.pairr2f<<"@"<<dataToWrite->dataToProcess->at(i).d2<<"/2" <<dataToWrite->dataToProcess->at(i).cmt<<endl <<dataToWrite->dataToProcess->at(i).s2<<endl<<"+"<<endl <<dataToWrite->dataToProcess->at(i).q2<<endl;
				onereadgroup.pairr1f<<"@"<<dataToWrite->dataToProcess->at(i).d1<<"/1" <<dataToWrite->dataToProcess->at(i).cmt<<endl <<dataToWrite->dataToProcess->at(i).s1<<endl<<"+"<<endl <<dataToWrite->dataToProcess->at(i).q1<<endl;

				continue;
		    
			    }else{
				if(dataToWrite->dataToProcess->at(i).sequence != ""){ //new sequence			    
				    onereadgroup.single<<"@"<<dataToWrite->dataToProcess->at(i).d1 << dataToWrite->dataToProcess->at(i).cmt<<endl << dataToWrite->dataToProcess->at(i).sequence<<endl<<"+"<<endl <<dataToWrite->dataToProcess->at(i).quality<<endl;    	    
			    
				    if( (dataToWrite->dataToProcess->at(i).sequence.length()+umif+umir) > max(dataToWrite->dataToProcess->at(i).s1.length(),dataToWrite->dataToProcess->at(i).s2.length()) ){
					mtr->incrementCountmergedoverlap();
				    }else{
					mtr->incrementCountmerged();			  
				    }

				}else{ //keep as is
				    mtr->incrementCountnothing();
				    if(trimKey){
					onereadgroup.pairr2<<"@"<<dataToWrite->dataToProcess->at(i).d2<<"/2"<<dataToWrite->dataToProcess->at(i).cmt <<endl <<dataToWrite->dataToProcess->at(i).s2.substr( key2.length() )<<endl<<"+"<<endl <<dataToWrite->dataToProcess->at(i).q2.substr( key2.length() )<<endl;
					onereadgroup.pairr1<<"@"<<dataToWrite->dataToProcess->at(i).d1<<"/1"<<dataToWrite->dataToProcess->at(i).cmt <<endl <<dataToWrite->dataToProcess->at(i).s1.substr( key1.length() )<<endl<<"+"<<endl <<dataToWrite->dataToProcess->at(i).q1.substr( key1.length() )<<endl;
				    }else{
					onereadgroup.pairr2<<"@"<<dataToWrite->dataToProcess->at(i).d2<<"/2"<<dataToWrite->dataToProcess->at(i).cmt <<endl <<dataToWrite->dataToProcess->at(i).s2                        <<endl<<"+"<<endl <<dataToWrite->dataToProcess->at(i).q2                        <<endl;
					onereadgroup.pairr1<<"@"<<dataToWrite->dataToProcess->at(i).d1<<"/1"<<dataToWrite->dataToProcess->at(i).cmt <<endl <<dataToWrite->dataToProcess->at(i).s1                        <<endl<<"+"<<endl <<dataToWrite->dataToProcess->at(i).q1                        <<endl;
				    }
			    
				}
			    }
		    
			}else{

		    
			    if(dataToWrite->dataToProcess->at(i).code != ' '){ //either chimera or missing key
			
				if(dataToWrite->dataToProcess->at(i).code == 'K'){
				    mtr->incrementCountfkey();
				}else{
				    if(dataToWrite->dataToProcess->at(i).code == 'D'){
					mtr->incrementCountchimera();
				    }else{
					if(dataToWrite->dataToProcess->at(i).code == 'U'){//UMI problem, sequence is too short
					    mtr->incrementCountUMIp();
					}else{					    					    
					    cerr << "leehom: Wrong return code =\""<<dataToWrite->dataToProcess->at(i).code<<"\""<<endl;
					    exit(1);
					}
				    }
				}
			
				onereadgroup.singlef<<"@"<<dataToWrite->dataToProcess->at(i).d1<<""<<dataToWrite->dataToProcess->at(i).cmt <<endl << dataToWrite->dataToProcess->at(i).s1 <<endl<<"+"<<endl << dataToWrite->dataToProcess->at(i).q1 <<endl;
				continue;
			    }
		    
			    if(dataToWrite->dataToProcess->at(i).sequence != ""){ //new sequence
				mtr->incrementCounttrimmed();
				onereadgroup.single<<"@"<<dataToWrite->dataToProcess->at(i).d1<<""<<dataToWrite->dataToProcess->at(i).cmt <<endl << dataToWrite->dataToProcess->at(i).sequence <<endl<<"+"<<endl << dataToWrite->dataToProcess->at(i).quality<<endl;
			    }else{
				mtr->incrementCountnothing();
				if(trimKey){
				    onereadgroup.single<<"@"<<dataToWrite->dataToProcess->at(i).d1<<""<<dataToWrite->dataToProcess->at(i).cmt <<endl << dataToWrite->dataToProcess->at(i).s1.substr( key1.length() )<<endl<<"+"<<endl << dataToWrite->dataToProcess->at(i).q1.substr( key1.length() )<<endl;
				}else{
				    onereadgroup.single<<"@"<<dataToWrite->dataToProcess->at(i).d1<<""<<dataToWrite->dataToProcess->at(i).cmt <<endl << dataToWrite->dataToProcess->at(i).s1                        <<endl<<"+"<<endl << dataToWrite->dataToProcess->at(i).q1                        <<endl;
				}
			    }
			}
			
		    }
		
		    wroteData=true;		
		    lastWrittenChunk=dataToWrite->rank;
		
		    if(dataToWrite->rank == lastRank)
			wroteEverything=true;	
		    delete dataToWrite;
		}else{//not the correct rank					
		    //cout<<"loop: is waiting for chunk "<<lastWrittenChunk<<" but got  "<<(dataToWrite->rank)<<endl;		    
		    rc = pthread_mutex_unlock(&mutexCounter);
		    checkResults("pthread_mutex_unlock()\n", rc);
		}

	    }else{
		rc = pthread_mutex_unlock(&mutexCounter);
		    checkResults("pthread_mutex_unlock()\n", rc);
	    }

	 
	    if(!wroteData){//if did not write all data, sleep

#ifdef DEBUGSLEEPING
		rc = pthread_mutex_lock(&mutexCERR);
		checkResults("pthread_mutex_lock()\n", rc);
		cerr<<"final: waiting to write, queue size= "<<qsize<<", sleeping for "<<timeThreadSleep<<"s"<<endl;
		rc = pthread_mutex_unlock(&mutexCERR);
		checkResults("pthread_mutex_unlock()\n", rc);
#endif

		sleep(timeThreadSleep);
	    }

	}//while did not write everything


	//waiting for threads to finish
	for (int i=0; i <numberOfThreads; ++i) {
	    rc = pthread_join(thread[i], NULL);
	    checkResults("pthread_join()\n", rc);
	}

	pthread_mutex_destroy(&mutexRank);
	pthread_mutex_destroy(&mutexQueue);
	pthread_mutex_destroy(&mutexCounter);
	pthread_mutex_destroy(&mutexCERR);



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
    
	
	cerr <<mtr->reportSingleLine()<<endl;

	if(printLog){
	    ofstream fileLog;
	    fileLog.open(logFileName.c_str());

	    if (fileLog.is_open()){
		fileLog <<mtr->reportMultipleLines() <<endl;

	    }else{
		cerr << "Unable to print to file "<<logFileName<<endl;
	    }
	    fileLog.close();
	}
    

    	//end fastq
    }else{
	//else BAM

	if( bamFileOUT == ""  ){
	    cerr<<"The output must be a be specified, exiting"<<endl;
	    return 1;
	}

	pthread_mutex_init(&mutexCERR,    NULL);
	pthread_mutex_init(&mutexQueue,   NULL);
	pthread_mutex_init(&mutexCounter, NULL);
	pthread_mutex_init(&mutexRank ,   NULL);
    
	pthread_t             thread[numberOfThreads];
	int                   rc=0;


	struct argspthread *argstopass = (struct argspthread *)malloc(sizeof(struct argspthread));
	argstopass->keepOrig        = keepOrig;
	argstopass->umif            = umif;
	argstopass->umir            = umir;

	for(int i=0;i<numberOfThreads;i++){
	    //	    rc = pthread_create(&thread[i], NULL, mainComputationThread, &keepOrig);
	    rc = pthread_create(&thread[i], NULL, mainComputationThread, (void *)argstopass);
	    checkResults("pthread_create()\n", rc);
	}


	// FastQParser * fqp1;
	// fqp1 = new FastQParser (fastqfile);
	BamReader reader;
	BamWriter writer;

	if ( !reader.Open(bamFile) ) {
	    cerr << "Could not open input BAM file" << bamFile<<endl;
	    return 1;
	}


	SamHeader       header     = reader.GetHeader();
	const RefVector references = reader.GetReferenceData();
    
	string pID          = "leeHom";   
	string pName        = "leeHom";   
	string pCommandLine = "";
	for(int i=0;i<(argc);i++){
	    pCommandLine += (string(argv[i])+" ");
	}
	putProgramInHeader(&header,pID,pName,pCommandLine,returnGitHubVersion(string(argv[0]),".."));

	//const RefVector references = reader.GetReferenceData();
	//we will not call bgzip with full compression, good for piping into another program to 
	//lessen the load on the CPU
	if(produceUnCompressedBAM) 
	    writer.SetCompressionMode(BamWriter::Uncompressed);

	if ( !writer.Open(bamFileOUT,header,references) ) {
	    cerr << "Could not open output BAM file "<<bamFileOUT << endl;
	    return 1;
	}

	// if ( !writer.Open(bamfilew,header,references) ) {
	// 	cerr << "Could not open output BAM file "<<bamfilew << endl;
	// 	return 1;
	// }


	// ogzstream outfile;
	// outfile.open(fastqfilew.c_str(), ios::out);
	// if(!outfile.good()){       cerr<<"Cannot write to file "<<fastqfilew<<endl; return 1; }
    
	unsigned int rank=0;

	DataChunk * currentChunk = new DataChunk(rank);
	unsigned int counter=0;
	int lastWrittenChunk=-1;
	// currentChunk->rank=rank;

	BamAlignment al;

	while ( reader.GetNextAlignment(al) ) {
	    
	    //FastQObj * fo1=fqp1->getData();
	    //cout<<counter<<"\t"<<currentChunk->dataToProcess->size()<<"\t"<<al.Name<<endl;
			    
	    if(al.IsMapped() || al.HasTag("NM") || al.HasTag("MD")  ){
		if(!allowAligned){
		    cerr << "Reads should not be aligned" << endl;
		    return 1;
		}else{
		    //should we remove tags ?
		}
	    }

	    if(counter== (sizeChunk-1)){
		
		//store old one
		currentChunk->dataToProcess->push_back(al);
	    

		rc = pthread_mutex_lock(&mutexQueue);
		checkResults("pthread_mutex_lock()\n", rc);
		
		while(true){
		    //mutex queue is taken here
		    if(int(queueDataToprocess.size())>(maxQueueDataToprocessSize)){	  // if queue full
#ifdef DEBUGSLEEPING
			int queueSize = int(queueDataToprocess.size());
#endif
			//release mutex queue
			rc = pthread_mutex_unlock(&mutexQueue);
			checkResults("pthread_mutex_unlock()\n", rc);


			bool hasWritten = checkForWritingLoop(&writer,&lastWrittenChunk);
			if(!hasWritten){
#ifdef DEBUGSLEEPING
			    rc = pthread_mutex_lock(&mutexCERR);
			    checkResults("pthread_mutex_lock()\n", rc);
			    cerr<<"main loop: queue to process is full, capacity: "<<queueSize<<", sleeping for "<<timeThreadSleep<<"s"<<endl;		
			    rc = pthread_mutex_unlock(&mutexCERR);
			    checkResults("pthread_mutex_unlock()\n", rc);
#endif

			    sleep(timeThreadSleep);

			}

			//take queue mutex
			rc = pthread_mutex_lock(&mutexQueue);
			checkResults("pthread_mutex_lock()\n", rc);
		    }else{
#ifdef DEBUG 
			rc = pthread_mutex_lock(&mutexCERR);
			checkResults("pthread_mutex_lock()\n", rc);			
			cerr<<"main loop: queue to process has space, adding rank: "<<currentChunk->rank<<endl;			
			rc = pthread_mutex_unlock(&mutexCERR);
			checkResults("pthread_mutex_unlock()\n", rc);
#endif

			queueDataToprocess.push(currentChunk);
			
			rc = pthread_mutex_unlock(&mutexQueue);
			checkResults("pthread_mutex_unlock()\n", rc);
			break;
		    }
		}

		rank++;

		checkForWritingLoop(&writer,&lastWrittenChunk);

		//new one
		currentChunk = new DataChunk(rank);	    
		//currentChunk->rank=rank;
		counter=0;
	    }else{
		currentChunk->dataToProcess->push_back(al);

		counter++;
	    }

	    // fqrecord toadd;
	    // toadd.ids  = *(fo1->getID());
	    // toadd.seqs = *(fo1->getSeq());
	    // toadd.qual = *(fo1->getQual());
	}//end bam file has data
	reader.Close();
		    			
#ifdef DEBUG 
	rc = pthread_mutex_lock(&mutexCERR);
	checkResults("pthread_mutex_lock()\n", rc);
	
	cerr<<"main loop: is done"<<endl;
	
	rc = pthread_mutex_unlock(&mutexCERR);
	checkResults("pthread_mutex_unlock()\n", rc);
#endif


	//add last chunk
	int lastRank=currentChunk->rank;
	queueDataToprocess.push(currentChunk);
	bool wroteEverything=false;

	readDataDone=true;
	//cout<<"main thread is done reading"<<endl;

	while(!wroteEverything){
	    //cout<<"write queue size = "<< queueDataTowrite.size() <<endl;
	    //threads are running here
	    rc = pthread_mutex_lock(&mutexCounter);
	    checkResults("pthread_mutex_lock()\n", rc);

#if defined(DEBUG) ||defined(DEBUGSLEEPING)
	    int qsize = queueDataTowrite.size(); 
#endif

	    bool wroteData=false;
	    if(!queueDataTowrite.empty()){
	
		DataChunk *  dataToWrite= queueDataTowrite.top();

		if( lastWrittenChunk == (dataToWrite->rank-1) ){
		    queueDataTowrite.pop();
		    rc = pthread_mutex_unlock(&mutexCounter);
		    checkResults("pthread_mutex_unlock()\n", rc);

#ifdef DEBUG 
		    rc = pthread_mutex_lock(&mutexCERR);
		    checkResults("pthread_mutex_lock()\n", rc);	       
		    cerr<<"final: writing "<<dataToWrite->rank<<" queue size= "<<qsize<<endl;
		    rc = pthread_mutex_unlock(&mutexCERR);
		    checkResults("pthread_mutex_unlock()\n", rc);
#endif

		    //writing dataToWrite
		    for(unsigned int i=0;i<dataToWrite->dataToProcess->size();i++){
			//cout<<"final: writing "<<dataToWrite->dataToProcess->at(i).Name<<endl;
			writer.SaveAlignment( dataToWrite->dataToProcess->at(i) );
			// outfile<< dataToWrite->dataToProcess->at(i).Name
			// 	   << endl  
			// 	   << dataToWrite->dataToProcess->at(i).QueryBases << endl
			// 	   << "+" <<endl 
			// 	   << dataToWrite->dataToProcess->at(i).Qualities << endl;
		    }
		
		    wroteData=true;		
		    lastWrittenChunk=dataToWrite->rank;
		
		    if(dataToWrite->rank == lastRank)
			wroteEverything=true;	
		    delete dataToWrite;
		}else{//not the correct rank					
		    //cout<<"loop: is waiting for chunk "<<lastWrittenChunk<<" but got  "<<(dataToWrite->rank)<<endl;		    
		    rc = pthread_mutex_unlock(&mutexCounter);
		    checkResults("pthread_mutex_unlock()\n", rc);
		}

	    }else{
		rc = pthread_mutex_unlock(&mutexCounter);
		    checkResults("pthread_mutex_unlock()\n", rc);
	    }

	 
	    if(!wroteData){//if did not write all data, sleep

#ifdef DEBUGSLEEPING
		rc = pthread_mutex_lock(&mutexCERR);
		checkResults("pthread_mutex_lock()\n", rc);     
		cerr<<"final: waiting to write, queue size= "<<qsize<<", sleeping for "<<timeThreadSleep<<"s"<<endl;
		rc = pthread_mutex_unlock(&mutexCERR);
		checkResults("pthread_mutex_unlock()\n", rc);
#endif

		sleep(timeThreadSleep);
	    }

	}//while did not write everything


	//waiting for threads to finish
	for (int i=0; i <numberOfThreads; ++i) {
	    rc = pthread_join(thread[i], NULL);
	    checkResults("pthread_join()\n", rc);
	}

	pthread_mutex_destroy(&mutexRank);
	pthread_mutex_destroy(&mutexQueue);
	pthread_mutex_destroy(&mutexCounter);
	pthread_mutex_destroy(&mutexCERR);
	writer.Close();

	cerr <<mtr->reportSingleLine()<<endl;

	if(printLog){
	    ofstream fileLog;
	    fileLog.open(logFileName.c_str());

	    if (fileLog.is_open()){
		fileLog <<mtr->reportMultipleLines() <<endl;

	    }else{
		cerr << "Unable to print to file "<<logFileName<<endl;
	    }
	    fileLog.close();
	}
    
	pthread_exit(NULL);
    }


    delete mtr;
    return 0;
}

