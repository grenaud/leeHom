/*
  Checks that the vectorized likelihood kernels give exactly the results of the original scalar code.

  Random reads (including very short ones, mates of different lengths, adapters with undefined bases,
  extreme qualities, ...) are processed with a range of options and a checksum of everything that
  comes out of process_PE() / process_SR() is printed. The checksum is the same whichever code is
  used:

      fuzz_kernels                       vectorized (AVX2 if the CPU has it)
      LEEHOM_NOAVX2=1 fuzz_kernels       vectorized, portable version
      LEEHOM_REFERENCE=1 fuzz_kernels    original scalar code

  test.sh runs the three and compares.
*/

#include <iostream>
#include <stdint.h>
#include <string>
#include <vector>
#include <cstdio>
#include <cstring>

#include "MergeTrimReads.h"

using namespace std;

static uint64_t rngState = 0x1234567887654321ULL;
static uint64_t rnd(){
    rngState ^= rngState >> 12;
    rngState ^= rngState << 25;
    rngState ^= rngState >> 27;
    return rngState * 2685821657736338717ULL;
}
static int rint_(int lo, int hi){ return lo + int(rnd() % uint64_t(hi-lo+1)); }

static uint64_t checksum = 14695981039346656037ULL;
static void mix(const string & s){
    for(size_t i=0;i<s.size();i++){
	checksum ^= (unsigned char) s[i];
	checksum *= 1099511628211ULL;
    }
    checksum ^= 0xff;
    checksum *= 1099511628211ULL;
}

static const char * BASES = "ACGT";

static string randomSeq(int len){
    string s;
    for(int i=0;i<len;i++)
	s += (rint_(0,40) == 0) ? 'N' : BASES[rint_(0,3)];
    return s;
}

static string randomQual(int len){
    string q;
    const int mode = rint_(0,3);
    for(int i=0;i<len;i++){
	switch(mode){
	case 0:  q += char(33 + rint_(2,41)); break;      //typical
	case 1:  q += char(33 + 40);          break;      //all high
	case 2:  q += char(33 + rint_(0,63)); break;      //whole supported range
	default: q += char(33 + ((rint_(0,5)==0) ? rint_(0,10) : rint_(30,41))); break;
	}
    }
    return q;
}

static string revcomp(const string & s){
    string r(s.rbegin(), s.rend());
    for(size_t i=0;i<r.size();i++){
	switch(r[i]){
	case 'A': r[i]='T'; break;
	case 'C': r[i]='G'; break;
	case 'G': r[i]='C'; break;
	case 'T': r[i]='A'; break;
	default: break;
	}
    }
    return r;
}

//a copy of a with a few errors
static string mutate(const string & a, int errors){
    string s = a;
    for(int e=0;e<errors && !s.empty();e++)
	s[rint_(0,int(s.size())-1)] = BASES[rint_(0,3)];
    return s;
}

int main(void){

    const string fwd[] = { "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACIIIIIIIATCTCGTATGCCGTCTTCTGCTTG",
			   "AGATCGGAAGAGCACACGTCTGAACTCCAG",
			   "AGATCGGAAGAGC",
			   "AGATCGGAAG",
			   "AGAT",
			   "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACIIIIIIIATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAAAAAAAAAAA" };
    const string rev[] = { "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT",
			   "GGAAGAGCGTCGTGTAGGGAAAGAGTGTAG",
			   "GGAAGAGCGTCG",
			   "GGAAGAGCGT",
			   "GGAA",
			   "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTTAAAAAAAAAAAAAAAAAA" };
    const string chimera = "ACACTCTTTCCCTACACGTCTGAACTCCAG,ACACTCTTTCCCACACGTCTGAACTCCAGT,CGTATGCCGTCTTCTGCTTGAAAAAAAAAA";

    long nPE = 0, nSR = 0, mergedPE = 0, trimmedSR = 0, flaggedPE = 0;
    long compared = 0, mismatches = 0;
    vector<double> llA, llB;

    const int quickChoices[8] = {0, 4, 8, 16, 0, 24, 12, 3};

    for(int config=0; config<48; config++){
	const int  ad          = config % 6;
	const bool ancient     = (config/6) % 2 == 1;
	const bool prior       = ((config/12) % 2 == 1) && !ancient;     // --loc/--scale excludes --ancientdna
	const int  trimCutoff  = 1 + (config % 3);

	MergeTrimReads m(fwd[ad], rev[ad], chimera,
			 "", "", false, "",
			 trimCutoff, false, ancient,
			 prior ? 5.14 : -1.0, prior ? 0.29 : -1.0, prior,
			 33, 42, true);
	m.setQuick(quickChoices[config % 8], (config % 3 == 0) ? 0.0 : ((config % 3 == 1) ? 0.5 : 0.9));

	for(int it=0; it<800; it++){
	    m.setRandomContext(uint64_t(config)*100000 + it);

	    //an insert, with the adapters (possibly cut) after it
	    const int insertLen = rint_(1, 220);
	    const string insert = randomSeq(insertLen);
	    const string a1     = mutate(fwd[ad].substr(0, rint_(0, int(fwd[ad].size()))), rint_(0,2));
	    const string a2     = mutate(rev[ad].substr(0, rint_(0, int(rev[ad].size()))), rint_(0,2));

	    //single-end read
	    {
		int len = rint_(2, 170);
		string s = (insert + a1 + randomSeq(170));
		s = s.substr(0, len);
		if(rint_(0,9) == 0) s = randomSeq(len);
		string q = randomQual(s.size());
		if(int(s.size()) > trimCutoff){
		    //the vectorized and the original code must give the very same values
		    m.debugLikelihoodsSR(s, q, true,  llA);
		    m.debugLikelihoodsSR(s, q, false, llB);
		    compared++;
		    if(llA.size() != llB.size() || (!llA.empty() && memcmp(&llA[0], &llB[0], llA.size()*sizeof(double)) != 0)){
			mismatches++;
			if(mismatches < 5) fprintf(stderr, "MISMATCH single-end read, config %d length %d\n", config, int(s.size()));
		    }
		    merged r = m.process_SR(s, q);
		    mix(string(1,r.code)); mix(r.sequence); mix(r.quality);
		    nSR++;
		    if(!r.sequence.empty()) trimmedSR++;
		}
	    }

	    //paired-end reads
	    {
		int len1 = rint_(1, 160);
		int len2 = (rint_(0,3) == 0) ? rint_(1,160) : len1;
		string frag = insert;
		string r1 = (frag + a1 + randomSeq(170)).substr(0, len1);
		string r2 = (revcomp(frag) + a2 + randomSeq(170)).substr(0, len2);
		if(rint_(0,15) == 0){ r1 = randomSeq(len1); r2 = randomSeq(len2); }
		string q1 = randomQual(r1.size());
		string q2 = randomQual(r2.size());
		m.debugLikelihoodsPE(r1, q1, r2, q2, true,  llA);
		m.debugLikelihoodsPE(r1, q1, r2, q2, false, llB);
		compared++;
		if(llA.size() != llB.size() || (!llA.empty() && memcmp(&llA[0], &llB[0], llA.size()*sizeof(double)) != 0)){
		    mismatches++;
		    if(mismatches < 5) fprintf(stderr, "MISMATCH paired-end reads, config %d lengths %d %d\n", config, int(r1.size()), int(r2.size()));
		}
		//the likelihood of each chimeric sequence
		m.debugChimera(r1, q1, true,  llA);
		m.debugChimera(r1, q1, false, llB);
		compared++;
		if(llA.size() != llB.size() || (!llA.empty() && memcmp(&llA[0], &llB[0], llA.size()*sizeof(double)) != 0)){
		    mismatches++;
		    if(mismatches < 5) fprintf(stderr, "MISMATCH chimeras, config %d length %d\n", config, int(r1.size()));
		}
		//the best and second best candidates, with their number of matches
		{
		    double bA, bB, sA, sB; int biA, biB, bmA, bmB, siA, siB, smA, smB;
		    m.debugBestPE(r1, q1, r2, q2, true,  bA, biA, bmA, sA, siA, smA);
		    m.debugBestPE(r1, q1, r2, q2, false, bB, biB, bmB, sB, siB, smB);
		    compared++;
		    if(memcmp(&bA,&bB,sizeof(double)) || memcmp(&sA,&sB,sizeof(double)) || biA!=biB || bmA!=bmB || siA!=siB || smA!=smB){
			mismatches++;
			if(mismatches < 5) fprintf(stderr, "MISMATCH in the best candidates, config %d lengths %d %d (idx %d/%d matches %d/%d)\n", config, int(r1.size()), int(r2.size()), biA, biB, bmA, bmB);
		    }
		}
		merged r = m.process_PE(r1, q1, r2, q2);
		mix(string(1,r.code)); mix(r.sequence); mix(r.quality);
		nPE++;
		if(!r.sequence.empty()) mergedPE++;
		if(r.code != ' ') flaggedPE++;
	    }
	}
    }

    printf("pairs %ld (merged %ld, flagged %ld) single %ld (trimmed %ld) checksum %016llx  likelihood vectors compared %ld, differing %ld\n", nPE, mergedPE, flaggedPE, nSR, trimmedSR, (unsigned long long) checksum, compared, mismatches);
    return (mismatches == 0) ? 0 : 1;
}
