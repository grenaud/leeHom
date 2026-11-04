#ifndef leeHom_Kernels_h
#define leeHom_Kernels_h

#include <cstddef>
#include <string>
#include <vector>

/*
  Likelihood kernels used by MergeTrimReads.

  These compute, for every candidate position at once, exactly the values that
  MergeTrimReads::detectAdapter(), MergeTrimReads::measureOverlap() and MergeTrimReads::detectChimera()
  return one candidate at a time. The floating-point additions for a given candidate are done
  in the same order as in those functions, so the results are bit-identical; the
  speedup comes from evaluating several independent candidates side by side (SIMD, with AVX2
  when the CPU has it).
*/

namespace leehom {

//! Number of candidates evaluated together by the overlap scan; the quick mode looks at a group at once
const int QUICK_GROUP = 8;

//! Pre-combined lookup tables, all indexed by the raw quality character (0..127)
struct KernelTables {
    //! [2][128] : [0*128+q] = likeMismatch[q], [1*128+q] = likeMatch[q]
    const double * adapterTab;
    //! [128][2][128] : [(q1*2+eq)*128+q2] = eq ? likeMatchPair[q1][q2] : likeMismatchPair[q1][q2]
    const double * overlapTab;
    //! likeRandomMatch
    double         randMatch;
    //! log10fastModeProbError
    double         fastModeThreshold;
    //! Maximum number of adapter bases compared
    size_t         maxAdapterComp;
    /*!
      --quick: 0 means off. Otherwise, a candidate overlap whose log-likelihood after quickBases compared bases is
      below what random sequences would give is not computed any further: the remaining terms are added as
      their expected value under random sequences (see quickExpected()), the number of terms is unchanged.
    */
    int            quickBases;
    //! how far, between the random expectation (0) and a perfect overlap (1), the likelihood must be after quickBases to go on
    double         quickFraction;
};

/*!
  Expected value of the log-likelihood of one pair of bases of qualities qa1 and qa2 if the sequences were random:
  3/4 of mismatches and 1/4 of matches. tab is KernelTables::overlapTab.
*/
inline double quickExpected(const double * tab, int qa1, int qa2);
inline double quickMatchTerm(const double * tab, int qa1, int qa2){
    return tab[size_t(qa1 & 127) * 256 + 128 + (qa2 & 127)];
}

/*!
  The likelihood (relative to the value before the first compared base) that a candidate must have reached after
  N compared bases to be computed any further: that of random sequences, plus fraction of the way to a perfect overlap.
*/
inline double quickThreshold(const double * tab, int qa1, int qa2, int N, double r, double fraction){
    const double E  = quickExpected(tab, qa1, qa2);
    const double Mt = quickMatchTerm(tab, qa1, qa2);
    double thr = double(N) * (E + r);
    thr += fraction * (double(N) * (Mt - E));
    return thr;
}

inline double quickExpected(const double * tab, int qa1, int qa2){
    const size_t row = size_t(qa1 & 127) * 256;
    return 0.75 * tab[row +       (qa2 & 127)] +
	   0.25 * tab[row + 128 + (qa2 & 127)];
}

/*!
  For each offset in [0,count), adds to ll[offset] the log-likelihood of observing
  the adapter at that offset of the read.
  Equivalent to ll[o] += detectAdapter(read,qual,adapter,o);
  Requires count <= len.
*/
void adapterScan(const KernelTables & kt,
		 const char * read,
		 const char * qual,
		 int          len,
		 const char * adapter,
		 int          adapterLen,
		 int          count,
		 double *     ll);

/*!
  For each offset in [0,count), adds to ll[offset] the log-likelihood of the overlap
  of read1 and the reverse-complemented read2 if the adapter starts at that offset.
  Equivalent to ll[o] += measureOverlap(read1,qual1,read2,qual2,.,o);
*/
void overlapScan(const KernelTables & kt,
		 const char * read1,
		 const char * qual1,
		 int          len1,
		 const char * read2rev,
		 const char * qual2rev,
		 int          len2,
		 int          count,
		 double *     ll);


//! The chimeric sequences, laid out so that all of them (with the offsets 0 and 1) can be compared to a read at once
struct ChimeraLanes {
    int                         nlanes;     //2 per sequence: offsets 0 and 1
    int                         stride;     //nlanes rounded up to a multiple of 8
    int                         steps;      //length of the longest sequence
    std::vector<unsigned char>  chars;      //[steps][stride] : base of the lane at that step, 0 if none
    std::vector<int>            length;     //[stride] : number of bases of the lane (sequence length - offset)

    ChimeraLanes():nlanes(0),stride(0),steps(0){}
    //! Lays out the sequences
    void build(const std::vector<std::string> & chimeras);
};

/*!
  out[lane] = detectChimera(read,qual,sequence,offset) for lane = 2*sequence + offset (-DBL_MAX when there are no bases to compare).
  out must have room for cl.nlanes values.
*/
void chimeraScan(const KernelTables & kt,
		 const ChimeraLanes & cl,
		 const char * read,
		 const char * qual,
		 int          len,
		 double *     out);

}

#endif
