#include "MergeTrimReads.h"


//#define DEBUG2
// #define DEBUGSR

// #define DEBUGADAPT
//#define DEBUGOVERLAP
// // //#define DEBUGTOTALOV
// // #define DEBUGPARTIALOV
// // #define CONSBASEPROB

//#define DEBUGPR
// #define DEBUGSCORE



const long double PI  = atanl(1.0L)*4;   


// vector<string> adapter_chimeras (chimInit,chimInit+13);



// double cutoff_merge_seqs_early = 0.95;
// double cutoff_merge_seqs = 0.90;
// bool   initialized       = false;
// double cutoff_merge_trim = 0.80;




//! Computes log of probability density function
/*!
  Computes the log base 10 of pdf(x)
  
  \param mu The mean (location)
  \param sigma The variance (scale)
  \param x The value for which we want the pdf
  \return The values of log base 10 of (pdf(x))
*/
long double MergeTrimReads::logcomppdf(long double mu,long double sigma,long double x){
    if(x==0){
	x=1.0;
    }

    long double two = 2.0;   
    long double exponent = logl(x) - mu;
    exponent *= -exponent;
    exponent /= two * sigma * sigma;
    
    long double result = exponent/logl(10);
    result -= logl(sigma * sqrtl(two * PI) * x)/logl(10);

    return result;
}


//! Computes log of cummulative distribution function
/*!
  Computes the log base 10 of 1 - cdf(x)
  
  \param mu The mean (location)
  \param sigma The variance (scale)
  \param x The value for which we want the cdf
  \return The values of log base 10 of (1 - cdf(x))
*/
long double MergeTrimReads::logcompcdf(long double mu,long  double sigma,long  double x){
    if(x==0){
	x=1.0;
    }

    long double two = 2.0;
    long double result = logl(x) - mu;
    result /= sigma * sqrtl( two );
    result =  erfcl(result); //erf(x)  = 1 - erfc(x)
    result *= 0.5;
    return logl(result)/logl(10);    
}


//! Computes reverse complement (char)
/*!
  Computes reverse complement of a char and returns it, returns 'N' for 'N'
  
  \param c the char to reverse complement
  \return The reverse complement
*/
inline char MergeTrimReads::revComp(char c){
    if(c ==    'A')
	return 'T';

    if(c ==    'C')
	return 'G';

    if(c ==    'G')
	return 'C';

    if(c ==    'T')
	return 'A';

    if(c ==    'N')
	return 'N';
    
    cerr<<"Wrong base in revComp function = "<<c<<endl;
    exit(1);
}


//! Computes reverse complement (string)
/*!
  Computes reverse complement of a string by calling revComp() on each character and returns it
  
  \param seq The string to reverse complement
  \return The reverse complement
*/
inline string MergeTrimReads::revcompl(const string seq){
    string toReturn="";
    int seqL=seq.length();
    int i=0;
    while(i<seqL){
	toReturn+=revComp(seq[seqL-i-1]);
	i++;
    }
    return toReturn;
}



//! Initializes likelihood scores
/*!
  Since quality scores can only take discrete values, we pre-compute
  the likelihood of a match for quality scores 0 to 64.
  For all scores, we can the log_10 of it.

  The likelihood of observing any given base is 1/4
  
  We first compute the likelihood scores for quality scores
  less than 2 (since they should not exist, we approximate them)
  We then compute the likelihood for the remaing quality scores

  For a given quality score i, we use the following formula
  	likeMatch(i)        =   1.0-10.0^(i/-10.0)
	    
        likeMismatch(i)     =      (10.0^(i/-10.0) )/3.0  

  For a pair of nucleotide, we must compute the likelihood of
  stemming from the same template. If the nucleotides are identical (say AA)
  then there are 4 possibilities :
         template    prob
         A           (right)(right)
         C           (wrong)(wrong)
         G           (wrong)(wrong)
         T           (wrong)(wrong)

  If the nucleotides differ (say AC) then we have the following:
         template    prob
         A           (right)(wrong)
         C           (wrong)(right)
         G           (wrong)(wrong)
         T           (wrong)(wrong)
  
  We pre-compute the likelihood for each. In the measureOverlap() subroutine, the prior for the nucleotide is added
  We also pre-compute the likelihood for the pdf and cdf for the lognormal distribution.
  

*/
void MergeTrimReads::initMerge(){

    likeRandomMatchProb          =  1.0/4.0 ;
    likeRandomMisMatchProb       =  3.0/4.0 ;

    likeRandomMatch       = log( 1.0/4.0 )/log(10);
    likeRandomMisMatch    = log( 3.0/4.0 )/log(10);

#ifdef DEBUG2
    cout<<"likeRandomMisMatch "<<likeRandomMisMatch<<endl;
    cout<<"likeRandomMatch    "<<likeRandomMatch<<endl;
#endif

    //probability scores for qscores less than 2 make no sense
    //since random DNA has by default a probability of error of 0.75
    //anyway, we set the min for qual scores at 2
    for(int i=0;i<2;i++){

	likeMatch[i]        = log1p(    -pow(10.0,2.0/-10.0) )    /log(10);	    
        likeMismatch[i]     = log  (     pow(10.0,2.0/-10.0)/3.0 )/log(10);

	likeMatchProb[i]           = 1.0-pow(10.0,2.0/-10.0) ;
        likeMismatchProb[i]        =     pow(10.0,2.0/-10.0)/3.0 ;
 
	probForQual[i]      = max(double(1.0)-pow(double(10.0),double(2.0)/double(-10.0)),
				  max_prob_N);

    }


    //Computing for quality scores 2 and up
    for(int i=2;i<64;i++){
	likeMatch[i]        = log1p(    -pow(10.0,i/-10.0) )     /log(10);	    
        likeMismatch[i]     = log  (     pow(10.0,i/-10.0)/3.0  )/log(10); //i/(-10.0*1.0);  


	likeMatchProb[i]           = 1.0-pow(10.0,i/-10.0);
	likeMismatchProb[i]        =     pow(10.0,i/-10.0)/3.0;

	probForQual[i] = max(double(1.0)-pow(double(10.0),double(i)/double(-10.0)),
			     max_prob_N);

#ifdef DEBUG2
        cout<<endl<<"qual = "<<i<<endl;
        cout<<"match    " <<likeMatch[i]<<endl;
        cout<<"mismatch " <<likeMismatch[i]<<endl;
	cout<<"match    " <<likeMatchProb[i]<<endl;
        cout<<"mismatch " <<likeMismatchProb[i]<<endl;
#endif

    }

    //pre-computing the likelihood for pairs of nucleotides
    for(int i=0;i<64;i++){
	for(int j=0;j<64;j++){	    
	    likeMatchPair[i][j]    = 
		log  ( 
		      (1.0*   likeMatchProb[i] *    likeMatchProb[j] +     //mm
		       3.0*likeMismatchProb[i] * likeMismatchProb[j] )     //ww
		       )/log(10);
	    likeMismatchPair[i][j] = 
		log  ( 
		      (1.0*   likeMatchProb[i] * likeMismatchProb[j] +     //mw
		       1.0*likeMismatchProb[i] *    likeMatchProb[j] +     //wm
		       2.0*likeMismatchProb[i] * likeMismatchProb[j] )     //ww
		       )/log(10);  
	}
    }

#ifdef DEBUG2
    cout<<"pair qual match"<<endl;
    cout<<"\t";
    for(int i=0;i<64;i++)
	cout<<i<<"\t";

    for(int i=0;i<64;i++){
	cout<<i<<"\t";
	for(int j=0;j<64;j++){
	    cout<<pow(10.0,likeMatchPair[i][j])<<"\t";	    
	}
	cout<<endl;
    }

    cout<<"pair qual mismatch"<<endl;
    cout<<"\t";
    for(int i=0;i<64;i++)
	cout<<i<<"\t";

    for(int i=0;i<64;i++){
	cout<<i<<"\t";
	for(int j=0;j<64;j++){
	    cout<<pow(10.0,likeMismatchPair[i][j])<<"\t";	    
	}
	cout<<endl;
    }


#endif


    if(useDist){
	//compute pdf and cdf for sequences	
	for(int i=0;i<MAXLENGTHSEQUENCE;i++){
	    pdfDist[i] = logcomppdf(location,scale, i ); 
	    cdfDist[i] = logcompcdf(location,scale, i ); 

#ifdef DEBUG2
	    cout<<"pdfDist["<<i<<"]="<<pow(10.0,pdfDist[i])<<"\t"<<"cdfDist["<<i<<"]="<<pow(10.0,cdfDist[i])<<endl;
#endif
	}
    }



}



//! Computes the ASCII string using the values of the PHRED quality scores
/*!
  Computes the ASCII string using the values of the PHRED quality scores supplied as vector of ints
  
  \param  logScores vector of ints of log scores
  \return The string ready to be written as a BAM
*/
inline string MergeTrimReads::convert_logprob_quality(vector<int> logScores){
    string toReturn="";
    for(unsigned int i=0;i<logScores.size();i++){
	toReturn+=char(max(qualOffset,min(126,logScores[i]+qualOffset)));
    }
    return toReturn;
}





//! Returns random double
/*!
  Initializes the seed if it wasn't initialized and 
  returns double(rand())/double(RAND_MAX).
  
  This is done because when there is a base with two different 
  potential bases and a tie in quality score, we must take one 
  at random
  
 
  \return Random double
*/
inline double MergeTrimReads::randomGen(){
    //    static bool initialized = false;

    if(initialized == false){
	struct timeval time;
	gettimeofday(&time,NULL);
	srand(time.tv_usec);
	initialized=true;
	return  (double(rand())/double(RAND_MAX));       
    }else{
	return  (double(rand())/double(RAND_MAX));       
    }
}


//! Computes the consensus for two individual bases with their quality scores
/*!
  
  Computes a new base and a new quality for two (somewhat)
  independent observations of a single base

  \param  base1 first base (as baseQual struct)
  \param  base2 second base (as baseQual struct)

  \return consensus base (as baseQual struct)
*/
inline baseQual MergeTrimReads::cons_base_prob(baseQual  base1,baseQual base2){
    const string dnaAlphabet  = "ACGT";
    
    vector<char> bases;
    if(base1.base == base2.base){
	bases.push_back(base1.base);
    }else{
	bases.push_back(base1.base);
	bases.push_back(base2.base);
    }

    double aprob1 = log10(1.0-base1.prob)-log10(3.0); 
    double lprob1 = log10(base1.prob);
    double lprob2 = log10(base2.prob);
    double aprob2 = log10(1.0-base2.prob)-log10(3.0);

#ifdef CONSBASEPROB
    cerr<<base1.base<<"\t"<<base1.prob<<"\t"<<base2.base<<"\t"<<base2.prob<<"\t"<<endl;
#endif

    double total_prob = 0.0;
    for(int i=0;i<4;i++){
	double help = 0.0;
	if(base1.base == dnaAlphabet[i])
	    help+=lprob1;
	else 
	    help+=aprob1;
	if(base2.base == dnaAlphabet[i])
	    help+=lprob2;
	else
	    help+=aprob2;
	total_prob+=pow(10.0,help);
    }
    total_prob = log10(total_prob);

    baseQual toReturn;
    toReturn.base='N' ;
    toReturn.qual= -1;

    for(unsigned int i=0;i<bases.size();i++){
	double thelp=0.0;
	for(int k=0;k<4;k++){
	    if(bases[i] != dnaAlphabet[k]){
		double help=0.0;
		if(base1.base == dnaAlphabet[k])
		    help+=lprob1;
		else 
		    help+=aprob1;

		if(base2.base == dnaAlphabet[k])
		    help+=lprob2;
		else
		    help+=aprob2;
		thelp+=pow(double(10.0),help);
	    }
	}
	thelp=log10(thelp);
	int hqual = int( floor( min(60.0,-10.0*(thelp-total_prob)) + 0.5) );
       	if( (hqual > toReturn.qual) || ((hqual == toReturn.qual) && (randomGen() >=0.5)) ){
	    toReturn.base = bases[i];
	    toReturn.qual = hqual;
	}
    }

    return toReturn;
}




//! Sets the adapters sequences
/*!
  
  Used by mergeTrimReadsBAM.cpp to set the adapter
  and chimeric sequences.

  \param forward Forward adapter
  \param reverse Reverse adapter
  \param chimera Potential chimeras (comma separated)
*/
void MergeTrimReads::set_adapter_sequences(const string& forward, const string& reverse, const string& chimera){
    options_adapter_F.assign( forward.length(), 0 );
    options_adapter_S.assign( reverse.length(), 0 );
    transform(forward.begin(), forward.end(), options_adapter_F.begin(), ::toupper);
    transform(reverse.begin(), reverse.end(), options_adapter_S.begin(), ::toupper);

    string tempChimera( chimera.length(), 0 );
    transform(chimera.begin(), chimera.end(), tempChimera.begin(), ::toupper);

    adapter_chimeras.clear();
    string token=returnFirstToken(&tempChimera,",");
    bool foundAdapterF=false;
    while(token.length()!=0){
	if(options_adapter_F == token)
	    foundAdapterF=true;
	adapter_chimeras.push_back(token);
	token=returnFirstToken(&tempChimera,",");	    
    }
    if(!foundAdapterF)
	adapter_chimeras.push_back(options_adapter_F);
    vector<string>::iterator it;
    it=adapter_chimeras.begin() ; 
    while(it < adapter_chimeras.end()  ){

	it->erase( remove( it->begin(), it->end(), ' ' ), it->end() ); //remove white spaces

	if( ( it->length() > 0) && (it->length()  <= maxadapter_comp)){
	    while ( it->length() <= maxadapter_comp ){
		*it += "I";
	    }
	    it++;
	}else{
	    if( it->empty() ){
		adapter_chimeras.erase(it);
	    }else{
		it++;
	    }
	}
    }
  
    
    while ( options_adapter_F.length() <= maxadapter_comp){
	options_adapter_F += "I";
    }

    while ( options_adapter_S.length() <= maxadapter_comp){
	options_adapter_S += "I";
    }


}
   


//! Tokenize string
/*!
  
  Tokenize string, returns the first token and 
  deletes it from the main string

  \param toparse main string to tokenize, will be altered by this subroutine
  \param delim   Field delimited for main string
  \return First token
*/
string MergeTrimReads::returnFirstToken(string * toparse,string delim){
    size_t found;
    string toreturn;
    found=toparse->find(delim);
    if (found!=string::npos){
	toreturn=toparse->substr(0,found);
	toparse->erase(0,found+delim.length());
	return toreturn;
    }
    toreturn=(*toparse);
    toparse->erase(0,toparse->length());    
    return toreturn;
}



//! Sets the key sequences
/*!
  
  Used by mergeTrimReadsBAM.cpp to set the key sequences


  \param key1 Forward key
  \param key2 Reverse key
*/
void MergeTrimReads::set_keys(const string& key1, const string& key2){

    keys0=key1;
    if(key2==""){
	keys1=key1;
    }else{
	keys1=key2;
    }

    len_key1=keys0.length();
    len_key2=keys1.length();

    transform(keys0.begin(), keys0.end(),keys0.begin(), ::toupper);
    transform(keys1.begin(), keys1.end(),keys1.begin(), ::toupper);

    if( (len_key1 > 0) || (len_key2 > 0))
	handle_key=true;
    else
	handle_key=false;
	
}

//! Sets options by driver program
/*!
  
  Used by mergeTrimReadsBAM.cpp to set options


  \param trimcutoff_     After how many bases do we trim a read ?
  \param allowMissing_   Whether we allow some missing bases in key
  \param ancientDNA_     Whether we merge partially overlapping reads
  
*/
void MergeTrimReads::set_options(int trimcutoff_,bool allowMissing_,  bool ancientDNA_ ){ //,bool mergeoverlap){
    options_trimCutoff   = trimcutoff_;
    options_allowMissing = allowMissing_;
    //options_mergeoverlap = mergeoverlap;
    ancientDNA = ancientDNA_;
}


//! Computes the likelihood of a given match for a read to a chimera
/*!
  This subroutine computes the likelihood of the following:

  L(observing first bases matching the chimeric string) 
  x 
  L(observing the remainder of the first bases)
  
  If this exceeds  L(observing the bases of the reads)
  we will flag it as a chimeric sequence

  \param read Actual bases of the read in which to detect the chimeric string
  \param qual Quality scores associated with the read
  \param chimeraString the actual chimeric string
  \param offsetChimera Whether we will allow a small offset in the chimeric string (e.g. 1)
  \return The likelihood of the observation 
*/
inline double MergeTrimReads::detectChimera(const string      & read,
					    const vector<int> & qual,
					    const string      & chimeraString,
					    unsigned int        offsetChimera){

    double likelihoodMatch=0.0;
    unsigned int maxidx=min(read.length(),chimeraString.length()-offsetChimera);

    if(maxidx <= 0)
	return -DBL_MAX;

    for(unsigned i=0;i<maxidx;i++){
	// cerr<<"i= "<<i<<endl;
	if(          read[i]                        == chimeraString[i+offsetChimera] || 
		     chimeraString[i+offsetChimera] == 'I' ){
	    likelihoodMatch  +=    likeMatch[ qual[i] ];
	}else{
	    likelihoodMatch  += likeMismatch[ qual[i] ];
	}
	// cerr<<"i= "<<likelihoodMatch<<endl;
    }
    //add likelihood of remaining bases
    //    cout<<double(read.length() - maxidx )<<endl;
    likelihoodMatch  += double(read.length() - maxidx ) * likeRandomMatch;
    
    return likelihoodMatch;
}





//! Computes the likelihood of a observing an overlap between two sequences
/*!

  This subroutine computes the likelihood of the following:

  L(observing the same bases in "read1" and "read2" for "offsetRead" bases)  
  = L(base 1 in read 1 matching  base 1 in read 2) * L(base 2 in read 1 matching  base 2 in read 2) * ... 
  
  We model the likelihood of matching as:
  L(base 1 in read 1 matching  base 1 in read 2)  =  P(A) * L(base 1 in read 1 matching  base 1 in read 2|original seq. had an A)
                                                     P(C) * L(base 1 in read 1 matching  base 1 in read 2|original seq. had an C)
                                                     P(G) * L(base 1 in read 1 matching  base 1 in read 2|original seq. had an G)
                                                     P(T) * L(base 1 in read 1 matching  base 1 in read 2|original seq. had an T)

  The L(base 1 in read 1 matching  base 1 in read 2|original seq. had an X) is pre-computed in likeMatchPair and likeMismatchPair
  for a match or mismatch, The prior P(X) of the base is simply 1/4.
  
						     
  \param read1 Actual bases of the first read 
  \param qual1 Quality scores associated with the first read
  \param read2 Actual bases of the second read that were previously reverse complemented
  \param qual2 Quality scores associated with the second read

  \param maxLengthForPair pre-computed maximum length of read1 and read2
  \param offsetRead Number of bases to compute
  \param matches Increase for each match found

  \return The likelihood of observing an overlap between two sequences
*/
inline double MergeTrimReads::measureOverlap(const string      & read1,
					     const vector<int> & qual1,
					     const string      & read2,
					     const vector<int> & qual2,
					     const int         & maxLengthForPair,
					     unsigned int      offsetRead,				    
					     //double *          iterations =0 ,
					     int  *            matches){
    
#ifdef DEBUGOVERLAP
    string comparedread1;
    string comparedread2;
#endif

    int i1=0; //start index read 1
    int i2=int(read2.size())-int(offsetRead);
    double likelihoodMatch=0.0;
    
    for(int i=0;i<int(offsetRead);i++){


	//Both have characters
	if(     i1<int(read1.size()) && 
		i2<int(read2.size()) && 
		i2>=0){
		
#ifdef DEBUGOVERLAP
	    // cerr<<"overlap "<<" read1["<<(i1)<<"] "<<read1[i1]<<"\tread2["<<i2<< "] "<<read2[i2]<<endl;
	    comparedread1+=read1[i1];
	    comparedread2+=read2[i2];
#endif

	    if(read1[i1] == read2[i2] ){	    
		(*matches)++;
		//likelihoodMatch  +=    likeMatch[    min(qual1[i1],qual2[i2])  ];	    
		likelihoodMatch  +=    likeMatchPair[   qual1[i1] ][ qual2[i2] ];
	    }else{
		//likelihoodMatch  +=    likeMismatch[   min(qual1[i1],qual2[i2])  ];
		likelihoodMatch  +=   likeMismatchPair[ qual1[i1] ][ qual2[i2] ];	     //[ min(qual1[i1],qual2[i2])  ];	    
	    }
	    likelihoodMatch  += likeRandomMatch; // this is added as a prior for all possible matches, since it is equal for all possible nucleotides
	    //                                      it can be factored out and added as a log

	    
	    i1++;
	    i2++;
	    continue;
	}


	//Only read1 has characters
	if( i1<int(read1.size()) ){

#ifdef DEBUGOVERLAP
	    //cerr<<"overlap "<<" read1["<<(i1)<<"] "<<read1[i1]<<"\tread2["<<i2<< "] "<<read2[i2]<<endl;
	    comparedread1+=read1[i1];
	    comparedread2+="-";
#endif

	    likelihoodMatch  += likeRandomMatch; 

	
	    i1++;
	    i2++;
	    continue;
	}

	//only read2 has characters
	if( i2<int(read2.size()) && 
	    i2>=0 ){
		
#ifdef DEBUGOVERLAP
	    //cerr<<"overlap "<<" read1["<<(i1)<<"] "<<read1[i1]<<"\tread2["<<i2<< "] "<<read2[i2]<<endl;
	    comparedread1+="-";
	    comparedread2+=read2[i2];
#endif
	    likelihoodMatch  += likeRandomMatch; 
	    
	    i1++;
	    i2++;
	    continue;
	}


	//none have characters, can occur while comparing strings for partial overlap for strings with different lengths
#ifdef DEBUGOVERLAP
	comparedread1+="-";
	comparedread2+="-";
#endif
	likelihoodMatch  += likeRandomMatch; 
	    
	i1++;
	i2++;
	continue;
	    
	cerr<<"Wrong state measureOverlap()"<<i<<"\t"<<i1<<"\t"<<i2<<endl;
	exit(1);
	    
    }	
	

#ifdef DEBUGOVERLAP   
    cerr<<"comparedread1 "<<comparedread1<<endl;
    cerr<<"comparedread2 "<<comparedread2<<endl;
    cerr<<"result        "<<likelihoodMatch<<endl;
#endif

    return likelihoodMatch;
}




    


//! Computes the likelihood of a observing the adapter in a read
/*!
  This subroutine computes the likelihood of observing the adapter in "read"
    at index "offsetRead".  

  L(observing last bases matching the adapter) 
  x 
  L(observing the bases in the first part of the sequence)
  
  If this exceeds  L(observing the bases of the reads)


  \param read Actual bases of the read in which to detect the chimeric string
  \param qual Quality scores associated with the read
  \param adapterString String of the adapter 
  \param offsetRead The index at which to start the matching
  \param matches Variable used as accumulator for the # of matches
  \return The likelihood of the observation 
*/
inline double MergeTrimReads::detectAdapter(const string      & read,
					    const vector<int> & qual,
					    const string      & adapterString,
					    unsigned int        offsetRead,
					    int              *  matches){
    

    double likelihoodMatch=0.0;
    
    unsigned maxIterations=min ( min( (read.length()-offsetRead),
				      adapterString.length())  , maxadapter_comp  );


#ifdef DEBUGADAPT
    string comparedread;
    cerr<<"offsetRead "<<offsetRead<<endl<<"read "<<read<<"\n"<<adapterString<<endl;
#endif

    unsigned  i=0;
    unsigned  indexRead=offsetRead;

    while(i<maxIterations){
	indexRead = i+offsetRead;


#ifdef DEBUGADAPT
	comparedread+=read[indexRead];
	cerr<<"apt0 "<<indexRead<<endl;
 	cerr<<"apt1 "<<i<<"\t"<<endl;
	cerr<<"apt2 "<<likelihoodMatch<<endl;
	cerr<<"apt3 "<<read[indexRead]<<endl;
	cerr<<"apt3 "<<qual[indexRead]<<endl;
	cerr<<"apt4 "<<likeMatch[ qual[indexRead] ]<<endl;
	cerr<<"apt5 "<<likeMismatch[ qual[indexRead] ]<<endl;
	cerr<<"apt6 "<<adapterString[i]<<endl;
#endif

	if(         read[indexRead]           == adapterString[i] || 
		    adapterString[i]          == 'I' ){ //match
	    (*matches)++;
	    likelihoodMatch  +=    likeMatch[ qual[indexRead] ];
	}else{ //mismatch
	    likelihoodMatch  += likeMismatch[ qual[indexRead] ];
	}

#ifdef DEBUGADAPT
	cerr<<"apt "<<i<<"\t"<<likelihoodMatch<<"\tq:"<<qual[indexRead]<<"\t"<<likeMatch[ qual[indexRead] ]<<"\t"<<likeMismatch[ qual[indexRead] ]<<endl;
#endif

	//(*iterations)++;
	i++;

    }

    int extraBases = max(0,int(read.length())-int(indexRead)-1);

#ifdef DEBUGADAPT
    cerr<<"detectAdapter1 "<<read<<"\n              "<<comparedread <<"\n              "<<adapterString<<"\tlm:"<<likelihoodMatch<<"\trl:"<<read.length()<<"\tor:"<<offsetRead<<"\tadd:"<<"\tmi:"<<maxIterations<<"\teb:"<<extraBases<<endl;
#endif

   
    //Adding the likelihood of the remaining bases
    likelihoodMatch += double( extraBases ) * likeRandomMatch;


#ifdef DEBUGADAPT
    cerr<<"detectAdapter2 "<< comparedread <<"\n              "<<adapterString<<"\tlm:"<<likelihoodMatch<<"\trl:"<<read.length()<<"\tor:"<<offsetRead<<endl;
    //cerr<<double(read.length()-indexRead)<<endl;
#endif

    return likelihoodMatch;
}



//! Computes the edit distance between 2 sequences
/*!
  This subroutine computes the Levenshtein distance
  between two strings

  \param seq1 First string to evaluate
  \param seq2 Second string to evaluate
  \return The # of mismatches
*/
inline int MergeTrimReads::edits(const string & seq1,const string & seq2){
    int lmin = min(seq1.length(),seq2.length());
    int lmax = max(seq1.length(),seq2.length());
    int dist = lmax-lmin;

    for(int pos=0;pos<lmin;pos++){
	if (seq1[pos] != seq2[pos]) 
	    dist+=1;	
    }

    return dist;
}


//! Checks if the length of the sequence is the same as the length of the quality string
/*!
  Checks if the length of the sequence is the same as the length of the quality string
  will exit() if this is the case
  
  \param seq DNA string to evaluate
  \param qual Quality score string to check
*/
inline void MergeTrimReads::sanityCheckLength(const string & seq,const string & qual){

    if( (seq.length() != qual.length()) ){
	cerr<<"MergeTrimReads: The reads and qualities must have equal lengths"<<endl;
	exit(1);
    }

}



//! Checks if the DNA string starts with the key for single-end
/*!
  Checks if a single-end DNA string starts with the key
  if so, it will trim it
  otherwise, it will fail it
  
  \param read1 DNA string to check for the key
  \param qual1 quality string of read1
  \param toReturn merged struct to hold the results
  \return true if the key was not observed, false otherwise
*/
inline bool MergeTrimReads::checkKeySingleEnd(string & read1,string & qual1,merged & toReturn){


    if( handle_key ){
	if ((read1.substr(0,len_key1) == keys0)  || //perfect match
	    (options_allowMissing && (edits(read1.substr(0,len_key1),keys0) == 1))  //1mm in first key	     
	    ){ //1mm in second key
	    read1 = read1.substr(len_key1,read1.size()-len_key1);
	    qual1 = qual1.substr(len_key1,read1.size()-len_key1);
	}else{
	    if(options_allowMissing && 
	       (read1.substr(0,len_key1-1) == keys0.substr(1,len_key1-1))  ){
		read1 = read1.substr(len_key1-1,read1.size()-(len_key1-1));
		qual1 = qual1.substr(len_key1-1,read1.size()-(len_key1-1));
	    }else{
		toReturn.code    ='K';
		toReturn.sequence="";
		toReturn.quality ="";		    
		return true;	       
	    }
	}
    }

    return false;
}



//! Checks if DNA string pairs starts with the key for paired-end
/*!
  Checks if the paired-end DNA strings start with their respective keys
  if so, it will trim them
  otherwise, it will fail them
  
  \param read1 DNA string of first mate to check for the key
  \param qual1 quality string of read1
  \param read2 DNA string of second mate to check for the key
  \param qual2 quality string of read2
  \param toReturn merged struct to hold the results
  \return true if the key was not observed, false otherwise
*/
inline bool MergeTrimReads::checkKeyPairedEnd(string & read1,
					      string & qual1,
					      string & read2,
					      string & qual2,
					      merged & toReturn){


    if( handle_key && read1.length() > 0){

	if (((read1.substr(0,len_key1) == keys0) && (read2.substr(0,len_key2) == keys1)) || //perfect match
	    (options_allowMissing && (edits(read1.substr(0,len_key1),keys0) == 1) && (read2.substr(0,len_key2) == keys1)) || //1mm in first key
	    (options_allowMissing && (read1.substr(0,len_key1) == keys0) && (edits(read2.substr(0,len_key2),keys1) == 1))
	    ){ //1mm in second key
	    read1 = read1.substr(len_key1,read1.size()-len_key1);
	    qual1 = qual1.substr(len_key1,read1.size()-len_key1);
	    read2 = read2.substr(len_key2,read2.size()-len_key2);
	    qual2 = qual2.substr(len_key2,read2.size()-len_key2);
	}else{
	    if(options_allowMissing && 
	       ( (len_key1>0?(read1.substr(0,len_key1-1) == keys0.substr(1,len_key1-1)):true) && (read2.substr(0,len_key2) == keys1)) ){
		read1 = read1.substr(len_key1-1,read1.size()-(len_key1-1));
		qual1 = qual1.substr(len_key1-1,read1.size()-(len_key1-1));
		read2 = read2.substr(len_key2,read2.size()-len_key2);
		qual2 = qual2.substr(len_key2,read2.size()-len_key2);
	    }else{
		if(options_allowMissing && 
		   (read1.substr(0,len_key1) == keys0) && (len_key2>0?(read2.substr(0,len_key2-1) == keys1.substr(1,len_key2-1)):true)){
		    read1 = read1.substr(len_key1,read1.size()-len_key1);
		    qual1 = qual1.substr(len_key1,read1.size()-len_key1);
		    read2 = read2.substr(len_key2-1,read2.size()-(len_key2-1));
		    qual2 = qual2.substr(len_key2-1,read2.size()-(len_key2-1));		    
		}else{
		    toReturn.code    ='K';
		    toReturn.sequence="";
		    toReturn.quality ="";		    
		    return true;
		}

	    }
	}
    }

    return false;

}


//! Checks if the first mate seems to be a chimera
/*!
  Calls the detectChimera() to compute the likelihood of all potential chimeric
  sequences, if it exceeds 

  \param read1 DNA string of first mate to check for the key
  \param qual1 quality string of read1
  \param read2 DNA string of second mate to check for the key
  \param qual2 quality string of read2
  \param toReturn merged struct to hold the results
  \return true if the key was not observed, false otherwise
*/
inline bool MergeTrimReads::checkChimera(const string & read1,
					 const vector<int> & qualv1,
					 merged & toReturn, 
					 const double & logLikelihoodTotal){


    for(unsigned int indexChimera=0;indexChimera<adapter_chimeras.size();indexChimera++){
	
	
	if( detectChimera( read1 , qualv1 , adapter_chimeras[indexChimera], 0 )  > logLikelihoodTotal ){
	    toReturn.code    ='D';
	    toReturn.sequence="";
	    toReturn.quality ="";		    
	    return true;
	}


	if( detectChimera( read1 , qualv1 , adapter_chimeras[indexChimera], 1 )  > logLikelihoodTotal ){
	    toReturn.code    ='D';
	    toReturn.sequence="";
	    toReturn.quality ="";
	    return true;
	}

    }

    return false;
}


//! Transforms the quality string into a vector of ints
/*!
  Transforms the quality string into a vector of integers
  where each ints is the numerical value of the quality score (with respect 
  to the qualOffset) 
  

  \param qual  The string of quality scores (input)
  \param qualv The vector of quality scores as integers (output)
*/
inline void MergeTrimReads::string2NumericalQualScores(const string & qual,vector<int> & qualv){

    for(unsigned int i=0;i<qual.length();i++){
	//this nonsense is to make old compilers happy
	int t2  = int(char( qual[i] ))-qualOffset;
	int t = max( t2,2);
	qualv.push_back( t  );
    }
    
}



//! Computes all likelihoods for single-end reads
/*!
  This subroutine computes all possible likelihoods
  for every possibility of adapters index.

  It calls detectAdapter() for every possible index 
  of the adapter and retains the most likely one.


  For a given putative index, we compute the likelihood of observing the following:
    
  
                                  L(adapter in read1)
                                             |--------> (adapter)
  read1                         ---------------------------------
  let :
  logLike1 = L(adapter in read1)

  We compute the following using this :
  logLike1 = L(seeing adapter bases in read1) + L(remaining bases) (computed using detectAdapter())
  
  The total likelihood for that index given by logLike1 
  
  
  \param read1   The single read
  \param qualv1  The quality scores for read1
  \param logLikelihoodTotal        Best likelihood found for all indices
  \param logLikelihoodTotalIdx     Index of best likelihood found for all indices
  \param sndlogLikelihoodTotal     Second Best likelihood found for all indices
  \param sndlogLikelihoodTotalIdx  Index of second best likelihood found for all indices

*/
inline void MergeTrimReads::computeBestLikelihoodSingle(const string      & read1,
							const vector<int> & qualv1,
							double & logLikelihoodTotal,
							int    & logLikelihoodTotalIdx,
							double & sndlogLikelihoodTotal,
							int    & sndlogLikelihoodTotalIdx){
    // sndlogLikelihoodTotalMatches = logLikelihoodTotalMatches;
	    

    for(unsigned int indexAdapter=0;
	indexAdapter<(read1.length()-options_trimCutoff);
	indexAdapter++){
	//double iterations=0;
	int    matches   =0;

#ifdef DEBUGSR
	cerr<<indexAdapter<<"\t"<<options_adapter_F<<endl;
#endif

	double logLike=
	    (double( indexAdapter ) * likeRandomMatch ) //likelihood of remaining bases
	    +
	    detectAdapter( read1 , qualv1 , options_adapter_F,indexAdapter , &matches );

	if(useDist){
	    //logLike += logcomppdf(location,scale,indexAdapter);
	    logLike += pdfDist[indexAdapter];
	}


	if(logLikelihoodTotal < logLike){
	    sndlogLikelihoodTotal        = logLikelihoodTotal;
	    sndlogLikelihoodTotalIdx     = logLikelihoodTotalIdx;
	    // sndlogLikelihoodTotalMatches = logLikelihoodTotalMatches;
	    
	    logLikelihoodTotal           = logLike;
	    logLikelihoodTotalIdx        = indexAdapter;
	    // logLikelihoodTotalMatches    = matches;

	}else{ 
	    if(sndlogLikelihoodTotal < logLike){ //not more likely than first but more likely than second
		sndlogLikelihoodTotal        = logLike;
		sndlogLikelihoodTotalIdx     = indexAdapter;
		// sndlogLikelihoodTotalMatches = matches;
	    }
	}

#ifdef DEBUGSR
	cerr<<indexAdapter<<"\t"<<logLike<<endl;
#endif
    }
}




//! Computes all likelihoods for paired-end reads
/*!
  This subroutine computes all possible likelihoods
  for every possibility of adapters index.


  For a given putative index, we compute the likelihood of observing the following:
    

          L(adapter in read2)+  L(likelihood overlap) + L(adapter in read1)
                                             |--------> (adapter)
  read1                         ---------------------------------
  read2      ---------------------------------
                      <--------| (adapter 2)
  let :
  logLike1 = L(adapter in read1)
  logLike2 = L(adapter in read2)
  logLike3 = L(adapter in overlapping region)

  We compute the following using this :
  logLike1 = L(seeing adapter bases in read1) + L(remaining bases) (detectAdapter())
  logLike2 = L(seeing adapter bases in read2) + L(remaining bases) (detectAdapter())
  logLike3 = L(bases matching read1/read2 overlap)   (computed using measureOverlap())
  
  The total likelihood for that index is:
  logLike1 + logLike2 + logLike3

  The overall likelihood is given by the sum of :
    - detectAdapter() for the first mate
    - detectAdapter() for the second mate
    - measureOverlap() for the overlapping portion of read1 and read2_rev
  
  \param read1       The first mate
  \param qualv1      The quality scores for the first mate
  \param read2       The second mate
  \param qualv2      The quality scores for the second mate
  \param read2_rev   The second mate
  \param qualv2_rev  The quality scores for the second mate

  \paramlengthRead1  Length of first mate
  \paramlengthRead2  Length of second mate
  \parammaxLengthForPair Max of both length

  \param logLikelihoodTotal        Best likelihood found for all indices
  \param logLikelihoodTotalIdx     Index of best likelihood found for all indices
  \param sndlogLikelihoodTotal     Second Best likelihood found for all indices
  \param sndlogLikelihoodTotalIdx  Index of second best likelihood found for all indices

*/

inline void MergeTrimReads::computeBestLikelihoodPairedEnd(const string &      read1,
							   const vector<int> & qualv1,

							   const string &      read2,
							   const vector<int> & qualv2,

							   const string &      read2_rev,
							   const vector<int> & qualv2_rev,

							   const int & lengthRead1,
							   const int & lengthRead2,
							   const int & maxLengthForPair,

							   double & logLikelihoodTotal,
							   int    & logLikelihoodTotalIdx,
							   int    & logLikelihoodTotalMatches,

							   double & sndlogLikelihoodTotal,
							   int    & sndlogLikelihoodTotalIdx,
							   int    & sndlogLikelihoodTotalMatches){
    


    //int lengthDiffR1_R2 =  (lengthRead1-lengthRead2);

    for(int indexAdapter=0; //let index adapters be the index of the potential P7/P5 adapter
	indexAdapter<(2*maxLengthForPair-min_overlap_seqs);
	indexAdapter++){
	double logLike=0.0;

 	// double iterations=0;
	int    matches   =0;

	if(indexAdapter > maxLengthForPair) //partial overlap
	    //if(!options_mergeoverlap) //no point in continuing 
	    if(!ancientDNA)
		break;
	    
	//adapters
 	if(indexAdapter<lengthRead1)
 	    logLike  +=  detectAdapter(    read1 ,     qualv1     , options_adapter_F,indexAdapter, &matches);

 	if(indexAdapter<lengthRead2)
 	    logLike  +=  detectAdapter(    read2 ,     qualv2     , options_adapter_S,indexAdapter, &matches);

	
	//overlap
	logLike  +=	 measureOverlap(   read1      ,
					   qualv1     ,
					   read2_rev  ,
					   qualv2_rev ,  
					   maxLengthForPair,					      
					   indexAdapter,
					   &matches);
	//adding prior
	if(useDist){
	    //logLike += logcomppdf(location,scale,indexAdapter);
	    logLike +=  pdfDist[indexAdapter];
	}

	

#ifdef DEBUGPR
    	cerr<<"idx: "<<indexAdapter<<"\ttl:\t"<<logLike<<"\tlr:"<<"\tit:"<<"\tma:"<<matches<<"\tlprior:\t"<<( logcomppdf( location,scale,indexAdapter ))<<endl;
    	//cerr<<logLike1<<"\t"<<logLike2<<"\t"<<logLike3<<endl;
#endif
	
    	if(logLikelihoodTotal < logLike){
    	    sndlogLikelihoodTotal        = logLikelihoodTotal;
    	    sndlogLikelihoodTotalIdx     = logLikelihoodTotalIdx;
	    sndlogLikelihoodTotalMatches = logLikelihoodTotalMatches;
	    
    	    logLikelihoodTotal           = logLike;
    	    logLikelihoodTotalIdx        = indexAdapter;
	    logLikelihoodTotalMatches    = matches;
    	}else{ 
    	    if(sndlogLikelihoodTotal < logLike){ //not more likely than first but more likely than second
    		sndlogLikelihoodTotal        = logLike;
    		sndlogLikelihoodTotalIdx     = indexAdapter;
		sndlogLikelihoodTotalMatches = matches;
    	    }
    	}
        

    }
    
}


//! Computes the consensus sequence for paired-end reads
/*!

  After having called computeBestLikelihoodPairedEnd(), the
  most likely position of the adapter is known. Using this
  information, this subroutine computes the consensus if reads
  are to be merged



  \param read1       The first mate
  \param qualv1      The quality scores for the first mate
  \param read2_rev   The second mate
  \param qualv2_rev  The quality scores for the second mate

  \param logLikelihoodTotal         Best likelihood found for all indices
  \param logLikelihoodTotalIdx      Index of best likelihood found for all indices
  \param logLikelihoodTotalMatches  # of matches for the best likelihood 

  \param sndlogLikelihoodTotal         Second Best likelihood found for all indices
  \param sndlogLikelihoodTotalIdx      Index of second best likelihood found for all indices
  \param sndlogLikelihoodTotalMatches  # of matches for the second best likelihood

  \param maxLengthForPair Max of both length
  \param toReturn merged struct (this will be returned to the main program)

*/

inline void MergeTrimReads::computeConsensusPairedEnd( const string & read1,
						       const vector<int> &   qualv1,
							      
						       const string & read2_rev,
						       const vector<int> & qualv2_rev,

							      
						       const double & logLikelihoodTotal,
						       const int    & logLikelihoodTotalIdx,
						       const int    & logLikelihoodTotalMatches,
					      
						       const double & sndlogLikelihoodTotal,
						       const int    & sndlogLikelihoodTotalIdx,
						       const int    & sndlogLikelihoodTotalMatches,
					      
						       const int & maxLengthForPair,
						       merged & toReturn){

    baseQual b1;
    baseQual b2;
    //test for merging:
    if( logLikelihoodTotalIdx != -1    &&                       //the status quo is not the most likely
	//(logLikelihoodTotal/sndlogLikelihoodTotal)      < maxLikelihoodRatio &&
	( (sndlogLikelihoodTotal - logLikelihoodTotal) < log10(maxLikelihoodRatio) ) &&
	logLikelihoodTotalMatches >= int(min_overlap_seqs)  ){       //sufficient # of mismatches for partial overlap (artifically to min_overlap_seqs for complete overlap)

	
	int i1=0; //start index read 1
	int i2=int(read2_rev.size())-int(logLikelihoodTotalIdx);

	string        newSeq  = "";
	vector<int>   newQual ;

	for(int i=0;i<int(logLikelihoodTotalIdx);i++){



	    if(     i1<int(read1.size()) && 
		    i2<int(read2_rev.size()) && 
		    i2>=0){
		
		b1.base = read1[i1];
		b1.prob = probForQual[ qualv1[i1] ];
		b1.qual = -1;

		b2.base = read2_rev[i2];
		b2.prob = probForQual[ qualv2_rev[i2] ]; 
		b2.qual = -1;

		baseQual RT    = cons_base_prob(b1,b2);
		newSeq         +=  RT.base ;
		newQual.push_back( RT.qual );

		i1++;
		i2++;
		continue;
	    }

	
	    if( i1<int(read1.size()) ){

		newSeq              +=	 read1[i1] ;
		newQual.push_back(      qualv1[i1] );

		i1++;
		i2++;
		continue;
	    }

	    if( i2<int(read2_rev.size()) && 
		i2>=0 ){
		

		newSeq              +=	read2_rev[i2] ;
		newQual.push_back(     qualv2_rev[i2] );
	    
		i1++;
		i2++;
		continue;
	    }

	    cerr<<"Wrong state in computeBestLikelihoodPairedEnd()"<<endl;
	    exit(1);
	    
	}



	toReturn.code    =' ';
	toReturn.sequence=newSeq;	   
	toReturn.quality =convert_logprob_quality(newQual);	

	// cout<<toReturn.sequence<<endl;
	//return toReturn;

	
    }else{

	toReturn.code    =' ';
	toReturn.sequence="";	   
	toReturn.quality ="";
	//return toReturn;    
    }
    
}




//! Processes a single read
/*!

  This subroutine :
    - Checks the key (if needed)
    - Checked for chimera
    - Computes the likelihood of various  adapter indices using computeBestLikelihoodSingle()
    - Trims if sufficient evidence is found
    
  \param read1       The single-end read
  \param qualv1      The quality scores for the read

  \return A merged struct, detailing the result

*/
merged MergeTrimReads::process_SR(string  read1, 
				  string  qual1){

    merged toReturn;


    //check for the key (if needed)
    if(checkKeySingleEnd(read1,qual1,toReturn)){
	return toReturn;
    }
    //end check for the key (if needed)

    sanityCheckLength(read1,qual1);



    vector<int> qualv1;
    string2NumericalQualScores(qual1,qualv1);
    

    //start detecting chimera //
    double logLikelihoodTotal     = double(read1.length())*likeRandomMatch;
    int logLikelihoodTotalIdx     = -1;


    if(checkChimera(read1,
		    qualv1,
		    toReturn, 
		    logLikelihoodTotal)){
	return toReturn;	
    }

    // end detecting chimera //


    // cerr<<"read1 2 "<<read1<<endl;

    //adding distance prior
    if(useDist){	
	//logLikelihoodTotal +=  logcompcdf(location,scale, (read1.length()-options_trimCutoff) );
	logLikelihoodTotal +=  cdfDist[ (read1.length()-options_trimCutoff) ];
    }

  


    //second best hit
    double sndlogLikelihoodTotal      = -DBL_MAX;
    int sndlogLikelihoodTotalIdx      = -1;

    // cerr<<"read1 3 "<<read1<<endl;

    computeBestLikelihoodSingle(read1,
				qualv1,
				logLikelihoodTotal,
				logLikelihoodTotalIdx,
				sndlogLikelihoodTotal,
				sndlogLikelihoodTotalIdx);

#ifdef DEBUGSR
     cerr<<logLikelihoodTotal<<"\t"<<endl;
#endif






     if( logLikelihoodTotalIdx != -1 &&
	 ( (sndlogLikelihoodTotal - logLikelihoodTotal) < log10(maxLikelihoodRatio) )
	 ){
	
        read1 = read1.substr(0,logLikelihoodTotalIdx);
        qual1 = qual1.substr(0,logLikelihoodTotalIdx);


	if( read1.length() < min_length){
	    toReturn.code    ='D';
	    toReturn.sequence="";	   
	    toReturn.quality ="";
	    return toReturn;	
	}else{
	    toReturn.code    =' ';
	    toReturn.sequence=read1;	   
	    toReturn.quality =qual1;	
	    return toReturn;
	}


    }


    toReturn.code    =' ';
    toReturn.sequence="";	   
    toReturn.quality ="";
    return toReturn;
  
}


//! Processes a paired-end read
/*!

  This subroutine :
    - Checks the key (if needed)
    - Checked for chimera
    - Computes the likelihood of various  adapter indices for both reads using computeBestLikelihoodPairedEnd()
    - Merges both reads if there is sufficient evidence using computeConsensusPairedEnd()

    


    \param read1       The first mate read
    \param qualv1      The quality scores for the first mate

    \param read2       The second mate read
    \param qualv2      The quality scores for the second mate

    \return A merged struct, detailing the type of result obtained
*/
merged MergeTrimReads::process_PE( string  read1,  string  qual1,
				   string  read2,  string  qual2){

    // //read1 is larger, N pad the second
    // if(read1.size() > read2.size() ){

    // 	for(int i=0;i< ( int(read1.size()) - int(read2.size()));i++){
    // 	    read2+="N";
    // 	    qual2+="#";
    // 	}
    // }

    // //read2 is larger, N pad the first
    // if(read1.size() < read2.size() ){

    // 	for(int i=0;i< ( int(read2.size()) - int(read1.size()));i++){
    // 	    read1+="N";
    // 	    qual1+="#";
    // 	}

    // }



    merged toReturn;



    sanityCheckLength(read1,qual1);
    sanityCheckLength(read2,qual2);

    if(checkKeyPairedEnd(read1,qual1,read2,qual2,toReturn)){
	return toReturn;
    }


    vector<int> qualv1;
    vector<int> qualv2;

    string2NumericalQualScores(qual1,qualv1);
    string2NumericalQualScores(qual2,qualv2);







    
    int lengthRead1 = int(read1.size());
    int lengthRead2 = int(read2.size());

    double logLikeSecondRead      = double(lengthRead2)*likeRandomMatch;
    double logLikeFirstRead       = double(lengthRead1)*likeRandomMatch;
    double logLikelihoodTotal     = logLikeFirstRead + logLikeSecondRead;


    int logLikelihoodTotalIdx     = -1;
    int logLikelihoodTotalMatches =0;

    //start detecting chimera //
    //    double lowChimeraLike=-DBL_MAX;


    if(checkChimera(read1,
		    qualv1,
		    toReturn, 
		    logLikeFirstRead)){ //we just look at the first read hence we need to only use the likelihood of the first read
	return toReturn;	
    }   
    // end detecting chimera //




    //computing rev compl for read 2
    string  read2_rev=  revcompl(read2);
    string  qual2_rev=           qual2;

    reverse(qual2_rev.begin(), 
	    qual2_rev.end());
    
    vector<int> qualv2_rev (qualv2);

    reverse(qualv2_rev.begin(), 
	    qualv2_rev.end());

   

    //COMPUTE 
    int maxLengthForPair=max(lengthRead1,lengthRead2);

    if(useDist){
	//double priorDist    = 1.0 - cdf(p, (2*maxLengthForPair-min_overlap_seqs)  );
	//logLikelihoodTotal += logcompcdf(location,scale, (2*maxLengthForPair-min_overlap_seqs) ); //log( priorDist )/log(10);
	logLikelihoodTotal += cdfDist[ (2*maxLengthForPair-min_overlap_seqs) ];
    }

  
#ifdef DEBUGPR
    // cerr<< (2*maxLengthForPair-min_overlap_seqs)<<endl;
    // cerr<< cdf(p, (2*maxLengthForPair-min_overlap_seqs)) <<endl;
    // cerr<< (1.0-cdf(p, (2*maxLengthForPair-min_overlap_seqs))) <<endl;

    cerr<<"idx: "<<-1<<"\ttl:\t"<<logLikelihoodTotal<<"\tlr:"<<"\tit:"<<"\tma:"<<0<<"\tlprior:\t"<<( (2*maxLengthForPair-min_overlap_seqs)   )<<endl;
#endif

    //second best hit
    double sndlogLikelihoodTotal      = -DBL_MAX;
    int sndlogLikelihoodTotalIdx      = -1;
    int sndlogLikelihoodTotalMatches  =  0;

#ifdef DEBUGPR
    cerr<<"fst: "<<read1<<endl<<"raw: "<<read2<<endl<<"rev: "<<read2_rev<<endl<<endl;
#endif
  

    
    // #ifdef DEBUGSCORE
    //     cerr<<"mergeo "<<"\t"<<logLikelihoodTotal<<"\t"<<logLikelihoodTotalIdx<<"\t"<<logLikelihoodTotalMatches<<endl;
    //     //cout<<logLikelihoodTotal<<"\t"<<logLikelihoodTotalMatches<<endl;
    // #endif
    


    
    
    computeBestLikelihoodPairedEnd(read1,
				   qualv1,

				   read2,
				   qualv2,

				   read2_rev,
				   qualv2_rev,
				   
				   lengthRead1,
				   lengthRead2,
				   maxLengthForPair,

				   logLikelihoodTotal,
				   logLikelihoodTotalIdx,
				   logLikelihoodTotalMatches,
				   
				   sndlogLikelihoodTotal,
				   sndlogLikelihoodTotalIdx,
				   sndlogLikelihoodTotalMatches);
	
    
    
    
#ifdef DEBUGSCORE
    cerr<<"mergeo "<<"\t"<<logLikelihoodTotal<<"\t"<<logLikelihoodTotalIdx<<"\t"<<logLikelihoodTotalMatches<<"\trt:\t"<<(sndlogLikelihoodTotal-logLikelihoodTotal)<<"\t"<<sndlogLikelihoodTotal<<"\t"<<sndlogLikelihoodTotalIdx<<"\t"<<(logLikelihoodTotal/sndlogLikelihoodTotal)<<endl;
    //cout<<logLikelihoodTotal<<"\t"<<logLikelihoodTotalMatches<<endl;
#endif
    

    
    computeConsensusPairedEnd( read1,
			       qualv1,

			       read2_rev,
			       qualv2_rev,


			       logLikelihoodTotal,
			       logLikelihoodTotalIdx,
			       logLikelihoodTotalMatches,
					      
			       sndlogLikelihoodTotal,
			       sndlogLikelihoodTotalIdx,
			       sndlogLikelihoodTotalMatches,
					      
			       maxLengthForPair,
			       toReturn);
    
    return toReturn;    


}





//! Subroutine to return the unique characters as a sorted string 
/*!
Subroutine to return the unique characters as a sorted string 

\param v input string
\return The sorted string of unique characters
*/
inline string MergeTrimReads::sortUniqueChar(string v){
    // vector<char> v( input.begin(), input.end() );
    string::iterator it;
    sort (v.begin(), v.end());
    char   previousChar='!';
    string toReturn    ="";
    for (it=v.begin(); it!=v.end();it++ ){
	if(previousChar == '!' ||
	   previousChar != *it){
	    previousChar=*it;
	    toReturn+=previousChar;
	}	
    }
    return toReturn;
}


//! Subroutine to add extra bam flags
/*!
  This subroutine sets the MERGEDBAMFLAG bam flag in for the BamAlignment.
  If the flag is already present, it will remove it

  \param al The BamAlignment object where the flags will be added
  \param f an int for the type of flag to add
  \return true iff successful
*/
bool MergeTrimReads::set_extra_flag( BamAlignment &al, int32_t f )
{
    char tp=0;
    if(!al.GetTagType(MERGEDBAMFLAG,tp)) {
        if(al.AddTag(MERGEDBAMFLAG,"i",f)) return true;
    }
    else {
        int v=0;
        switch(tp) {
            case 'i': case 'I': case 's': case 'S': case 'c': case 'C': 
                al.GetTag(MERGEDBAMFLAG,v);
                al.RemoveTag(MERGEDBAMFLAG);
                if(al.AddTag(MERGEDBAMFLAG,"i",v|f)) return true;
                break ;
            default:
                cerr << "ERROR: " << MERGEDBAMFLAG << " has type " << tp << endl;
        }
    }
    cerr << "Unable to add tag " << MERGEDBAMFLAG << endl;
    return false;
}




//! General contructor
/*!
  This constructor initialize basic variables to begin using proceePair or processSingle()

  \param forward_ Sequence of the adapters for the adapter seen in the first read
  \param reverse_ Sequence of the adapters for the adapter seen in the second read
  \param chimera_ A string of comma separated substring representing putative chimeras
  \param key1_ If keying was used, sequence of the key for the first read
  \param key2_ If keying was used, sequence of the key for the seonc read
  \param trimcutoff The minimum number of bases we should have observed for single reads
  \param allowMissing boolean, if set to true, it will allow imperfect matches to the key
  \param mergeoverlap boolean, if set to true, it will allow merging of reads where the adapter was never seen and the reconstructed sequence is longer than the initial reads
*/
MergeTrimReads::MergeTrimReads (const string& forward_, const string& reverse_, const string& chimera_,
				const string& key1_, const string& key2_,
				int trimcutoff_,bool allowMissing_,
				bool ancientDNA_,double location_,double scale_,bool useDist_):
    //bool mergeoverlap_) : 
    min_length (5),
    qualOffset (33),



    //FLAGs for merged reads
    MERGEDBAMFLAG     ( "FF" ),
    TRIMMEDFLAG       ( 1 ),
    MERGEDFLAG        ( 2 ),
    TRIMMEDMERGEDFLAG ( 3 ),
    
    location(location_),
    scale(scale_),
    useDist(useDist_)
 {
     
     initialized = false;

    
     //TODO: the default values should be stored in the config.json file, not here
     max_prob_N = 0.25;
    
     maxLikelihoodRatio = 1.0/20.0;
    
     maxadapter_comp  = 30; /**< maximum number of bases to be compared in the adapter */
     min_overlap_seqs = 10; /**< maximum number that have to match in the case of partial overlap */ 
    
     //     min_length =5;
     //     qualOffset =33;

     count_all = 0;
     count_fkey = 0;
     count_merged = 0;
     count_merged_overlap = 0;
     count_trimmed = 0;
     count_nothing = 0;
     count_chimera = 0;

     checkedTags.push_back("RG");
     checkedTags.push_back("XI");
     checkedTags.push_back("YI");
     checkedTags.push_back("XJ");
     checkedTags.push_back("YJ");


    

     //  chimInit[]= {
     // 	"ACACTCTTTCCCTACACGTCTGAACTCCAG",
     // 	"ACACTCTTTCCCACACGTCTGAACTCCAGT",
     // 	"ACACTCTTTCCCTACACACGTCTGAACTCC",
     // 	"CTCTTTCCCTACACGTCTGAACTCCAGTCA",
     // 	"GAAGAGCACACGTCTGAACTCCAGTCACII",
     // 	"GAGCACACGTCTGAACTCCAGTCACIIIII",
     // 	"GATCGGAAGAGCACACGTCTGAACTCCAGT",
     // 	"AGATCGGAAGAGCACACGTCTGAACTCCAG",
     // 	"AGAGCACACGTCTGAACTCCAGTCACIIII",
     // 	"ACACGTCTGAACTCCAGTCACIIIIIIIAT",
     // 	"GTGCACACGTCTGAACTCCAGTCACIIIII",
     // 	"AGCACACGTCTGAACTCCAGTCACIIIIII",
     // 	"CGTATGCCGTCTTCTGCTTGAAAAAAAAAA"};

     //     options_adapter_F="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACIIIIIIIATCTCGTATGCCGTCTTCTGCTTG";
     //     options_adapter_S="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT";

     set_adapter_sequences(forward_,reverse_,chimera_);

     //adapter_chimeras = vector<string>(chimInit,chimInit+13);

     set_keys(key1_,key2_);

     // //  Key variables ///
     //     handle_key           = false;
     //     keys0="";
     //     keys1="";
     //     len_key1=0;
     //     len_key2=0;
     set_options(trimcutoff_,allowMissing_,ancientDNA_);

     if(useDist){
	 //p = lognormal_distribution<>(location,scale);	
	 //getting rid of arbitrary cutoffs
	 min_overlap_seqs   = 1; 
	 options_trimCutoff = 1;
	 ancientDNA         = true;
     }

     // size_t  options_trimCutoff   = 
     // bool    options_mergeoverlap = false;
     // bool    options_allowMissing = false;

     initMerge();
    
}


//! Destructor
/*!
  Destructor, no memory is allocated on the heap 
*/
MergeTrimReads::~MergeTrimReads (){

}


//! Processes the reads for paired-end reads 
/*!
  This subroutine will call process_PE() on the sequence and return a pair<> of BamAlignment objects
  corresponding to the pair of reads. If the reads have been merged, the second object will be empty.
  It will also fix the BAM tags accordingly.

  \param al BamAlignment object as input representing the first pair
  \param al2 BamAlignment object as input representing the first pair
  \return True if the pair was merged and is contained within al only
*/
//pair<BamAlignment,BamAlignment> MergeTrimReads::processPair(const BamAlignment & al,const BamAlignment & al2){
bool MergeTrimReads::processPair( BamAlignment & al, BamAlignment & al2){
    
    string read1;
    string read2;
    string qual1;
    string qual2;

    if(al.Name != al2.Name ){
	cerr << "Seq#1 has a different id than seq #2, exiting " << endl;
	exit(1);
    }
    count_all ++;
    if(al.IsFirstMate()  &&
       al2.IsSecondMate() ){
	read1 =string(al.QueryBases);
	read2 =string(al2.QueryBases);
	qual1 =string(al.Qualities);
	qual2 =string(al2.Qualities);		    
    }else{
	if(al2.IsFirstMate()  &&
	   al.IsSecondMate() ){
	    read1 =string(al2.QueryBases);
	    qual1 =string(al2.Qualities);
	    read2 =string(al.QueryBases);
	    qual2 =string(al.Qualities);		    
	}else{
	    cerr << "Seq#1 must be the first mate for seq #2, exiting " << endl;
	    exit(1);
	}
    }
    
    if(qual1 == "*"){
	qual1=string(read1.length(),'0');
    }
    if(qual2 == "*"){
	qual2=string(read1.length(),'0');
    }
    
    // cerr<<"read1 "<<read1<<endl;
    // cerr<<"read2 "<<read2<<endl;

    merged result=process_PE(read1,qual1,
			     read2,qual2);

    //     if(al.Name == "M00518_0167_000000000-A3Y1J_CH_A2109:1:1101:12584:1649"){
    // 	cerr<<"SEQ #"<<result.sequence<<"#"<<endl;
    // 	exit(1);
    //     }

    if(result.code != ' '){ //keys or chimeras
	string prevZQ1="";
	string prevZQ2="";
	//	pair<BamAlignment,BamAlignment> toReturn (al,al2);

	al.SetIsFailedQC(true);
	al.GetTag("ZQ",prevZQ1);		    
	prevZQ1+=result.code;
	if(al.HasTag("ZQ") ){ //this is done because bamtools was not intelligent enough to understand that "ZQ:A" becomes "ZQ:Z" when you add a char, oh well.. 
	    al.RemoveTag("ZQ");
	    if(al.HasTag("ZQ") ){
		cerr << "Failed to remove tag for "<< al.Name<<endl;
		exit(1);
	    }
	}

	if(prevZQ1 != "")
	    if(!al.AddTag("ZQ","Z",sortUniqueChar(prevZQ1))){
		cerr << "Error while editing tags new tag11:"<<prevZQ1 <<"#"<< endl;
		exit(1);
	    }

		   

	al2.SetIsFailedQC(true);
	al2.GetTag("ZQ",prevZQ2);
	prevZQ2+=result.code;
	if(al2.HasTag("ZQ") ){ 
	    al2.RemoveTag("ZQ");
	    if(al2.HasTag("ZQ") ){ 
		cerr << "Failed to remove tag for "<< al2.Name<< endl;
		exit(1);
	    }
	}

	if(prevZQ2 != "")
	    if(!al2.AddTag("ZQ","Z",sortUniqueChar(prevZQ2))){
		cerr << "Error while editing tags new tag21:" << prevZQ2<<"#"<<endl;
		exit(1);
	    }

		  
		    
	if( result.code == 'K'){
	    count_fkey ++;
	}else{
	    if( result.code  == 'D'){
		count_chimera++;
	    }else{
		cerr << "mergeTrimReadsBAM: Wrong return code =\""<<result.code<<"\""<<endl;
		exit(1);
	    }
	}
	return false;
    }





    //    cerr<<"new #"<<result.sequence<<"#"<<endl;

    if(result.sequence != ""){ //new sequence
	//BamAlignment toWrite (al);//build from the previous one
	string towriteZQ="";
	al.GetTag("ZQ",towriteZQ);   //get from the first one
	// string towriteZQ2="";
	// al2.GetTag("ZQ",towriteZQ2);  //get from the first one


	al.AlignmentFlag=4; 	//not failed
	
	al.SetIsFailedQC( al.IsFailedQC() && al2.IsFailedQC() ); //fail the new one if both fail 
	
	al.Position    = -1;
	al.MapQuality  =  0;
	//RNEXT
	al.MatePosition=-1;
	

	al.QueryBases = result.sequence;
	al.Qualities  = result.quality;
	// toWrite.SetIsMapped(false);
		

	//copy tag info
	// toWrite.TagData=al.TagData; 	
	if( result.code == ' '){
	    //cerr<<"TEST2 #"<<result.code<<" "<<endl<<result.sequence<<endl<<read1<<endl;
	    if( result.sequence.length() > max(read1.length(),read2.length())){
		//cerr<<"test3"<<endl;
		count_merged_overlap ++;			  
		if( !set_extra_flag(al,  MERGEDFLAG) ) exit(1);
	    }else{
 		//cerr<<"test4"<<endl;
		count_merged++;
		if( !set_extra_flag(al,  TRIMMEDMERGEDFLAG) ) exit(1);
	    }
	}
	//  	cerr<<"TEST #"<<result.code<<"#"<<endl;
	// 	exit(1);

	//toWrite.RemoveTag("ZQ");
	if(al.HasTag("ZQ") ){ 
	    al.RemoveTag("ZQ");
	    if(al.HasTag("ZQ") ){ 
		cerr << "Failed to remove tag for new "<< al.Name<< endl;
		exit(1);
	    }
	}

	if( al.QueryBases.length()  < min_length){
	    al.SetIsFailedQC( true );		   
	    // if(!al.EditTag("ZQ","Z",string("L"))){
	    // 	  cerr << "Error while editing tags" << endl;
	    // 	  exit(1);
	    // }
	    towriteZQ+="L";
	}

	/////////////////////////////
	//       Fixing tags       //
	/////////////////////////////
	string dummy1;
	string dummy2;
		  

	//paranoid check to make sure our tags are identical
	for(size_t idx=0;idx<checkedTags.size();idx++)
	    if(al.GetTag(checkedTags[idx],dummy1)){
		if(al2.GetTag(checkedTags[idx],dummy2)){
		    if(dummy1 != dummy2){
			cerr << "Value for "<<checkedTags[idx]<<" cannot differ between mates " << endl;
			exit(1);
		    }
		    //fine otherwise
		}else{
		    cerr << "One read has been assigned a  "<<checkedTags[idx]<<" tag but not the other " << endl;
		    exit(1);

		}
	    }

	//The new read has the same ZQ tag as al. at this point
	if(al.HasTag("ZQ") || al2.HasTag("ZQ") ){
	    if( al2.HasTag("ZQ") ){
		if(!al.GetTag("ZQ",dummy1)) {
		    cerr << "Failed to get ZQ field from read 1" << endl;
		    exit(1);
		}
		if(!al2.GetTag("ZQ",dummy2)) {
		    cerr << "Failed to get ZQ field from read 2" << endl;
		    exit(1);
		}
		//we then need to add the ZQ from the second read
		if(dummy1 != dummy2){
		    towriteZQ+=dummy2;
		}	       	    
	    }
	}


	if(towriteZQ != ""){
	    if(!al.AddTag("ZQ","Z",sortUniqueChar(towriteZQ))){
		cerr << "Error while editing tags new tag20 :"<<towriteZQ <<"#"<< endl;
		exit(1);
	    }
	}	
	//we keep the original reads
	// 	if(keepOrig){
	// 	    al.SetIsDuplicate(true);
	// 	    al2.SetIsDuplicate(true);
	// 	    writer.SaveAlignment(al2);
	// 	    writer.SaveAlignment(al);
	// 	}
	// BamAlignment empty;
	// return 	pair<BamAlignment,BamAlignment>(toWrite,empty);
	return true;
	    
    }else{ //keep as is
	
	if( result.code == ' ')
	    count_nothing++;
	//if we use the Dup flag ourselves, clear it
	// 	if(keepOrig){ 
	// 	    al.SetIsDuplicate(false);
	// 	    al2.SetIsDuplicate(false);
	// 	}
	//keep the sequences as pairs
	
	// return 	pair<BamAlignment,BamAlignment>(al,al2);

	//writer.SaveAlignment(al2);
	//writer.SaveAlignment(al);
	return false;

    }

}



/*!
  This subroutine will call process_SR() on the sequence and return a new BamALignment object
  corresponding to the potentially trimmed sequence.
  It will fix the BAM tags accordingly.
  \param al BamAlignment object as input and that will be written out
  \return The BamAlignment to be written
*/
void MergeTrimReads::processSingle(BamAlignment & al){

    string read1;
    string qual1;

    count_all ++;

    read1 =string(al.QueryBases);
    qual1 =string(al.Qualities);
    if(qual1 == "*"){
	qual1=string(read1.length(),'0');
    }
	
    // cerr<<read1<<endl;
    // cerr<<"-----"<<endl;
    // cerr<<qual1<<endl;

    merged result=process_SR(read1,qual1);

    if(result.code != ' '){ //either chimera or missing key
	string prevZQ1="";
	//BamAlignment toWrite (al);//build from the previous one

	al.SetIsFailedQC(true);
	al.GetTag("ZQ",prevZQ1);
	prevZQ1+=result.code;
	if(al.HasTag("ZQ") ){ 
	    al.RemoveTag("ZQ");
	    if(al.HasTag("ZQ") ){ 
		cerr << "Failed to remove tag for "<< al.Name<< endl;
		exit(1);
	    }
	}

	if(prevZQ1 != ""){
	    if(!al.EditTag("ZQ","Z",sortUniqueChar(prevZQ1))){
		cerr << "Error while editing tags new tag11:"<<prevZQ1 <<"#"<< endl;
		exit(1);
	    }
	}	
	    
	if( result.code == 'K'){
	    count_fkey ++;
	}else{
	    if( result.code  == 'D'){
		count_chimera++;
	    }else{
		cerr << "mergeTrimReadsBAM: Wrong return code =\""<<result.code<<"\""<<endl;
		exit(1);
	    }

	}
	//return toWrite;
    }

    if(result.sequence != ""){ //new sequence

	///BamAlignment toWrite (al);//build from the previous al
	al.MapQuality=0;

	al.QueryBases = result.sequence;
	al.Qualities  = result.quality;
	al.SetIsMapped(false);
	if(!set_extra_flag( al, TRIMMEDFLAG ))
	    exit(1); 

	al.Position    =-1;
	al.MatePosition=-1;		    
	if( result.code == ' ')
	    count_trimmed++;

	// 	if(keepOrig){
	// 	    al.SetIsDuplicate(true);
	// 	    writer.SaveAlignment(al);
	// 	}
	//return toWrite;
	
    }else{
	if( result.code == ' ')
	    count_nothing++;

	//return al;
    }


}


//! Returns a summary report as a string
/*!
  Returns a summary report as a string of all operations performed by this object since its creation
  
  \return A string reprenting the tally for the reads
*/
string MergeTrimReads::reportSingleLine(){
    return "Total " + stringify( count_all ) +"; Merged (trimming) "+ stringify(count_merged  ) +"; Merged (overlap) "+ stringify(count_merged_overlap  ) +"; Kept PE/SR "+ stringify( count_nothing ) +"; Trimmed SR "+ stringify(count_trimmed  ) +"; Adapter dimers/chimeras "+ stringify(count_chimera  ) +"; Failed Key "+ stringify(count_fkey  ) ;

}


//! Returns a full report as a string
/*!
  Returns a full report as a string of all operations performed by this object since its creation
  
  \return A string reprenting the tally for the reads
*/
string MergeTrimReads::reportMultipleLines(){
    return  "Total reads :"            +stringify(count_all)+           "\t"+stringify(100.0*double(count_all)           /double(count_all))+"%\n"+
      	    "Merged (trimming) "       +stringify(count_merged)+        "\t"+stringify(100.0*double(count_merged)        /double(count_all))+"%\n"+
	    "Merged (overlap) "        +stringify(count_merged_overlap)+"\t"+stringify(100.0*double(count_merged_overlap)/double(count_all))+"%\n"+
	    "Kept PE/SR "              +stringify(count_nothing)+       "\t"+stringify(100.0*double(count_nothing)       /double(count_all))+"%\n"+
	    "Trimmed SR "              +stringify(count_trimmed)+       "\t"+stringify(100.0*double(count_trimmed)       /double(count_all))+"%\n"+
	    "Adapter dimers/chimeras " +stringify(count_chimera)+       "\t"+stringify(100.0*double(count_chimera)       /double(count_all))+"%\n"+
	"Failed Key "              +stringify(count_fkey)+          "\t"+stringify(100.0*double(count_fkey)          /double(count_all))+"%\n";
}





int MergeTrimReads::getCountall(){            return count_all; }
int MergeTrimReads::getCountfkey(){           return count_fkey; }
int MergeTrimReads::getCountmerged(){         return count_merged; }
int MergeTrimReads::getCountmergedoverlap(){  return count_merged_overlap; }
int MergeTrimReads::getCounttrimmed(){        return count_trimmed; }
int MergeTrimReads::getCountnothing(){        return count_nothing; }
int MergeTrimReads::getCountchimera(){        return count_chimera; }


void MergeTrimReads::incrementCountall(){             count_all++;            }
void MergeTrimReads::incrementCountfkey(){            count_fkey++;           }
void MergeTrimReads::incrementCountmerged(){          count_merged++;         }
void MergeTrimReads::incrementCountmergedoverlap(){   count_merged_overlap++; }
void MergeTrimReads::incrementCounttrimmed(){         count_trimmed++;        }
void MergeTrimReads::incrementCountnothing(){         count_nothing++;        }
void MergeTrimReads::incrementCountchimera(){         count_chimera++;        }
