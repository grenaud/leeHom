
#ifndef MergeTrimReads_h
#define MergeTrimReads_h

#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <vector>
#include <cfloat>
#include <math.h>
#include <sys/time.h>

#include <stdint.h>
/* #include <boost/math/distributions/lognormal.hpp> */

#include "libgab.h"
#include "Kernels.h"
#include "BamRecord.h"

using namespace std;
/* using namespace boost::math; */

#define MAXLENGTHSEQUENCE 1000 //should be good for a few years for Illumina at least ...
#define MAXQCSCORE 64


typedef struct{
    char   code;
    string sequence;
    string quality;
} merged;

typedef struct {
    char   base;
    int    qual;
    double prob;
} baseQual;

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

class MergeTrimReads{
 private:
    //VARIABLES
    double maxLikelihoodRatio  ;
    double log10maxLikelihoodRatio  ;

    double fastModeProbError  ;
    double log10fastModeProbError  ;

    double likelihoodChimera  ;
    double likelihoodAdapterSR ;
    
    double likelihoodAdapterPR ;
    /* static const double likelihoodAdapterPR = -1.0; */
    bool initialized ;

    const size_t min_length   ;
    const int    qualOffset    ;
    bool trimKey;

    /* static const double max_prob_N = 0.25; */
    /* extern double cutoff_merge_trim; */
    size_t maxadapter_comp;

    size_t min_overlap_seqs;

    unsigned int seed;
    bool seedSpecified;

    //random numbers (used to break ties), see setRandomContext()
    uint64_t baseSeed;
    uint64_t rngState;
    bool     legacyRng;     //use the global rand(), results then depend on the order in which reads are processed

    /* //  Key variables /// */
    bool handle_key;
    bool options_allowMissing;
    string keys0;
    string keys1;
    string ikey;
    int len_key1;
    int len_key2;
    int len_ikey;

    size_t options_trimCutoff;
    //bool options_mergeoverlap;
    bool ancientDNA ;
    double max_prob_N ;
    /* extern size_t min_length ; */

    //Chimera options and adapter
    //char*  chimInit[];/* = { */
     
    vector<string> adapter_chimeras ;
    string options_adapter_F;
    string options_adapter_S;
    


    // //  Key variables ///
    /*     bool handle_key; */
    /*     string keys0; */
    /*     string keys1; */
    /*     int len_key1 */
    /*     int len_key2; */


    //likelihood variables
    //the quality offset is either 33 or 64
    double likeMatch[64+MAXQCSCORE];
    double likeMismatch[64+MAXQCSCORE];
    double likeMatchProb[64+MAXQCSCORE];
    double likeMismatchProb[64+MAXQCSCORE];

    /* double likeMatch33[33+64]; */
    /* double likeMismatch33[33+64]; */
    /* double likeMatchProb33[33+64]; */
    /* double likeMismatchProb33[33+64]; */

    /* double likeMatch64[64+64]; */
    /* double likeMismatch64[64+64]; */
    /* double likeMatchProb64[64+64]; */
    /* double likeMismatchProb64[64+64]; */

    double likeMatchPair[64+MAXQCSCORE][64+MAXQCSCORE];
    double likeMismatchPair[64+MAXQCSCORE][64+MAXQCSCORE];
    /* double likeMatchPair33[33+64][33+64]; */
    /* double likeMismatchPair33[33+64][33+64]; */
    /* double likeMatchPair64[64+64][64+64]; */
    /* double likeMismatchPair64[64+64][64+64]; */

    double probForQual[64+MAXQCSCORE];
    /* double probForQual33[33+64]; */
    /* double probForQual64[64+64]; */
    
    double likeRandomMatch;    // 1/4
    double likeRandomMisMatch; // 3/4
    double likeRandomMatchProb;    // 1/4
    double likeRandomMisMatchProb; // 3/4

    double newprob[5][5][64+MAXQCSCORE][64+MAXQCSCORE];

    //combined tables for the vectorized kernels, see Kernels.h
    double adapterTab[2*(64+MAXQCSCORE)];
    double overlapTab[(64+MAXQCSCORE)*2*(64+MAXQCSCORE)];
    int    quickBases;
    double quickFraction;
    bool   referenceKernels;      //use the original one-candidate-at-a-time code (for verification)
    

    //prior dist
    long double pdfDist[MAXLENGTHSEQUENCE];    
    long double cdfDist[MAXLENGTHSEQUENCE];
    
    //vector<string> adapter_chimeras;

    //FUNCTIONS
    string returnFirstToken(string * toparse,string delim);
    char revComp(char c);
    //inline string convert_logprob_quality(vector<int> logScores);
    inline double randomGen();
    inline baseQual cons_base_prob(baseQual  base1,baseQual base2);
    inline baseQual cons_base_probInit(baseQual  base1,baseQual base2);



    /* double computePDF(const double x); */
    /* double computeCDF(const double x); */
    
    void    setLikelihoodScores(double likelihoodChimera_,
				double likelihoodAdapterSR_,
				double likelihoodAdapterPR_);

    void set_options(int trimcutoff=1,bool allowMissing=false,bool mergeoverlap=false);
    void set_adapter_sequences(const string& forward, const string& reverse, const string& chimera);
    void set_keys(const string& key1="", 
		  const string& key2="");

    void initMerge();
    double detectChimera(const string      & read,
			 //const vector<int> & qual,
			 const string      & qual,
			 const string      & chimeraString,
			 unsigned int        offsetChimera=0);
    double measureOverlap(const string      & read1,
			  const string      & qual1,
			  const string      & read2,
			  const string      & qual2,
			  const int         & maxLengthForPair,
			  unsigned int      offsetRead=0,				    
			  //double *          iterations =0 ,
			  int  *            matches=0);
    double measureOverlapQuick(const string & read1,
			       const string & qual1,
			       const string & read2,
			       const string & qual2,
			       unsigned int   offsetRead);
    double detectAdapter(const string      & read,
			 const string      & qual,
			 const string      & adapterString,
			 unsigned int        offsetRead=0,
			 int              *  matches=0);

    long double logcomppdf(long double mu,long double sigma,long double x);
    long double logcompcdf(long double mu,long double sigma,long double x);


    int edits(const string & seq1,const string & seq2);
    void sanityCheckLength(const string & seq,const string & qual);
    bool checkKeySingleEnd(string & read1,string & qual1,merged & toReturn);
    bool checkKeyPairedEnd(string & read1,string & qual1,
			   string & read2,string & qual2,
			   merged & toReturn);
    bool checkChimera(const string & read1,
		      const string & qualv1,
		      merged & toReturn, 
		      const double & logLikelihoodTotal);
    //void string2NumericalQualScores(const string & qual,vector<int> & qualv);
    void computeBestLikelihoodSingle(const string      & read1,
				     const string      & qualv1,
				     double & logLikelihoodTotal,
				     int &    logLikelihoodTotalIdx,
				     double & sndlogLikelihoodTotal,
				     int &    sndlogLikelihoodTotalIdx);
    void computeBestLikelihoodSingleReference(const string      & read1,
					      const string      & qualv1,
					      double & logLikelihoodTotal,
					      int &    logLikelihoodTotalIdx,
					      double & sndlogLikelihoodTotal,
					      int &    sndlogLikelihoodTotalIdx);
    leehom::KernelTables kernelTables() const;
    leehom::ChimeraLanes chimeraLanes;     //adapter_chimeras laid out for the vectorized kernel

    void computeBestLikelihoodPairedEnd(const string &      read1,
					const string &      qualv1,
					    
					const string &      read2,
					const string &      qualv2,
					    
					const string &      read2_rev,
					const string &      qualv2_rev,
					    
					const int & lengthRead1,
					const int & lengthRead2,
					const int & maxLengthForPair,

					double & logLikelihoodTotal,
					int    & logLikelihoodTotalIdx,
					int    & logLikelihoodTotalMatches,
					    
					double & sndlogLikelihoodTotal,
					int    & sndlogLikelihoodTotalIdx,
					int    & sndlogLikelihoodTotalMatches);


    int countMatchesPE(const string & read1, const string & qualv1,
		       const string & read2, const string & qualv2,
		       const string & read2_rev, const string & qualv2_rev,
		       const int & lengthRead1, const int & lengthRead2, const int & maxLengthForPair,
		       int indexAdapter);

    void computeBestLikelihoodPairedEndReference(const string &      read1,
						 const string &      qualv1,
						 const string &      read2,
						 const string &      qualv2,
						 const string &      read2_rev,
						 const string &      qualv2_rev,
						 const int & lengthRead1,
						 const int & lengthRead2,
						 const int & maxLengthForPair,
						 double & logLikelihoodTotal,
						 int    & logLikelihoodTotalIdx,
						 int    & logLikelihoodTotalMatches,
						 double & sndlogLikelihoodTotal,
						 int    & sndlogLikelihoodTotalIdx,
						 int    & sndlogLikelihoodTotalMatches);

    void computeConsensusPairedEnd( const string & read1,
				    const string & qualv1,
					
				    const string & read2_rev,
				    const string & qualv2_rev,
					
							      
				    const double & logLikelihoodTotal,
				    const int    & logLikelihoodTotalIdx,
				    const int    & logLikelihoodTotalMatches,
					
				    const double & sndlogLikelihoodTotal,
				    const int    & sndlogLikelihoodTotalIdx,
				    const int    & sndlogLikelihoodTotalMatches,
					
				    const int & maxLengthForPair,
				    merged & toReturn);


    string sortUniqueChar(string v);
    bool set_extra_flag( BamRecord &al, int32_t f );

    const string MERGEDBAMFLAG ;
    const int32_t TRIMMEDFLAG       ;
    const int32_t MERGEDFLAG        ;
    const int32_t TRIMMEDMERGEDFLAG ;
 
    unsigned int count_all ;
    unsigned int count_fkey ;
    unsigned int count_merged ;
    unsigned int count_merged_overlap ;
    unsigned int count_trimmed ;
    unsigned int count_nothing ;
    unsigned int count_chimera ;
    unsigned int count_UMIp ;//UMI problems

    vector<string> checkedTags;

    //lnnorm
    double location;
    double scale;
    bool   useDist;
    //lognormal_distribution<> p;
    void bamAlgnToString( BamRecord & al ,
			  BamRecord & al2,
			  string & read1,
			  string & read2,
			  string & qual1,
			  string & qual2);
    void inferadaptPair(   string & read1,string & read2,string & qual1,string & qual2,vector<string> & adapterSeqs_fwd, vector< string > & qualadaptSeq_fwd,vector<string> & adapterSeqs_rev, vector< string  > & qualadaptSeq_rev);

 public:
    MergeTrimReads (const string& forward_, const string& reverse_, const string& chimera_,
		    const string& key1_="", const string& key2_="",const bool trimKey_=false,const string& ikey_="",
		    const int trimcutoff_=1,
		    const bool allowMissing_=false,
		    const bool ancientDNA_=false,
		    const double location_=-1.0,
		    const double scale_=-1.0,
		    const bool useDist_=false,
		    const int    qualOffset=33,
		    const unsigned int seed=0,const bool seedSpecified=false);

    ~MergeTrimReads();
    
    merged process_PE(string  read1,string  qual1,string read2,string qual2);
    merged process_SR(string  read1, string qual1);

    /*!
      Sets the state of the random number generator from the index of the cluster about to be processed.
      The random numbers drawn for a given cluster then only depend on the seed and on that index, so the
      results are the same whatever the number of threads or the way reads are split into chunks.
    */
    void setRandomContext(uint64_t clusterIndex);

    //! Adds the counts of another object (used to combine the counts of the threads)
    void addCounts(const MergeTrimReads & other);

    /*!
      For verification only: the log-likelihood (before the prior) and the number of matches for every
      candidate position of the adapter, with either the original scalar code (reference=true) or the
      vectorized kernels. Both must give exactly the same values.
    */
    void debugLikelihoodsPE(const string & read1, const string & qual1,
			    const string & read2, const string & qual2,
			    bool reference, vector<double> & ll);
    void debugLikelihoodsSR(const string & read1, const string & qual1,
			    bool reference, vector<double> & ll);
    /*!
      For verification only: what the search for the most likely adapter index returns for a read pair (best and second
      best log-likelihood after the prior, their index and number of matches), with the original or the vectorized code.
    */
    void debugChimera(const string & read1, const string & qual1, bool reference, vector<double> & out);
    void debugBestPE(const string & read1, const string & qual1,
		     const string & read2, const string & qual2,
		     bool reference,
		     double & best, int & bestIdx, int & bestMatches,
		     double & second, int & secondIdx, int & secondMatches);

    /*!
      Quick mode: a candidate adapter position whose overlap likelihood, after bases compared bases, is below what
      random sequences would give is not computed any further (0 turns it off). See Kernels.h.
    */
    void setQuick(int bases, double fraction = 0.4){ quickBases = (bases > 0) ? bases : 0; quickFraction = fraction; }

    //! Use the original scalar likelihood code instead of the vectorized kernels (same results, slower)
    void setReferenceKernels(bool b){ referenceKernels = b; }

    void setFadapter(const string& forward_);
    void setRadapter(const string& reverse_);
		    
    /* pair<BamRecord,BamRecord> processPair(const BamRecord & al,const BamRecord & al2); */
    /* BamRecord                    processSingle(const BamRecord & al); */
    string mlConsensus(vector<string> & adapterSeqs, vector< string  > & qualadaptSeq,const float threshold=0.99);
    void inferadaptPairFQ(fqrecord & fr, vector<string> & adapterSeqs_fwd, vector< string  > & qualadaptSeq_fwd,vector<string> & adapterSeqs_rev, vector< string > & qualadaptSeq_rev);

    void inferadaptPairBAM(   BamRecord & al , BamRecord & al2,vector<string> & adapterSeqs_fwd, vector< string  > & qualadaptSeq_fwd,vector<string> & adapterSeqs_rev, vector< string  > & qualadaptSeq_rev);
    bool processPair(   BamRecord & al , BamRecord & al2);
    void processSingle( BamRecord & al );

    string reportSingleLine();
    string reportMultipleLines();

    string revcompl(const string seq);


    int getCountall();
    int getCountfkey();
    int getCountmerged();
    int getCountmergedoverlap();
    int getCounttrimmed();
    int getCountnothing();
    int getCountchimera();
    int getCountUMIp();


    void incrementCountall();
    void incrementCountfkey();
    void incrementCountmerged();
    void incrementCountmergedoverlap();
    void incrementCounttrimmed();
    void incrementCountnothing();
    void incrementCountchimera();
    void incrementCountUMIp();

};
#endif
