#include <algorithm>
#include <vector>
#include <string>
#include <cstring>
#include <cstdint>
#include <cstdlib>
#include <cfloat>

#include "Kernels.h"

using namespace std;

namespace leehom {

static const int G   = 8;                // candidates evaluated together by the portable code
static const int PAD = 64;               // padding on both sides of the scratch copies

static inline unsigned qidx(char c){ return unsigned((unsigned char)c) & 127u; }

//! Padded scratch copy of a sequence/quality so that inactive lanes can read past the ends safely
struct Padded {
    vector<char> buf;
    const char * p;   // points at element 0
    void assign(const char * src, int len){
	buf.assign(size_t(len) + 2*PAD, 0);
	if(len > 0)
	    memcpy(&buf[PAD], src, len);
	p = &buf[PAD];
    }
};

//! prefix[n] = 0.0 + r + r + ... (n times, one addition at a time, as in measureOverlap())
static const double * randPrefix(double r, size_t n){
    static thread_local vector<double> prefix;
    static thread_local double         cachedR = 0.0;
    if(prefix.empty() || cachedR != r){
	prefix.assign(1, 0.0);
	cachedR = r;
    }
    while(prefix.size() <= n)
	prefix.push_back(prefix.back() + r);
    return &prefix[0];
}


//! What the quick mode needs to decide about a candidate: prefix sums of the quality characters of both reads
struct QuickCtx {
    const KernelTables * kt;
    const int *          P1;      // P1[x] = sum of the first x quality characters of read 1
    const int *          P2;      // same for the (reverse complemented) read 2, in the order of read2rev
    int                  len1, len2, N;
    uint32_t             Minv[QUICK_GROUP];   // Minv[d] = ceil(2^32/(N+d)), to divide the sums by N+d with a multiplication
};

static void quickInit(QuickCtx & c, const KernelTables & kt, int len1, int len2, const vector<int> & P1, const vector<int> & P2){
    c.kt = &kt; c.len1 = len1; c.len2 = len2; c.N = kt.quickBases;
    c.P1 = &P1[0]; c.P2 = &P2[0];
    for(int d=0; d<QUICK_GROUP; d++)
	c.Minv[d] = uint32_t( ((uint64_t(1) << 32) + (c.N + d) - 1) / (c.N + d) );
}

//! Fills the prefix sums
static void quickPrefix(const char * qual, int len, vector<int> & P){
    P.assign(size_t(len) + 1, 0);
    for(int i=0;i<len;i++)
	P[i+1] = P[i] + int((unsigned char) qual[i]);
}

//! Step after which the candidates of the group of k0 (k0 is a multiple of QUICK_GROUP) are looked at
static inline int quickCheckpoint(const QuickCtx & c, int k0){
    return max(0, k0 + QUICK_GROUP - 1 - c.len2) + c.N;
}

//! True if the candidate k has valid bases left after the checkpoint cp
static inline bool quickApplies(const QuickCtx & c, int k, int cp){
    return min(k, c.len1) > cp;
}

/*!
  Called once the candidate k has compared its first bases up to the step cp, s0 is the value before any of them. If it
  is below what random sequences give, s gets the rest of the terms as their expected value and true is returned.
  (sum / n is done with a multiplication, exact because the sums are far below 2^32/n)
*/
static inline bool quickCheck(const QuickCtx & c, int k, int cp, double s0, double & s){
    const int    A   = max(0, k - c.len2);
    const int    Ev  = min(k, c.len1);
    const int    nv  = Ev - A;
    const int    n   = cp - A;                         // compared bases, between N and N+QUICK_GROUP-1
    const int    i2s = c.len2 - k + A;
    const uint64_t M = c.Minv[n - c.N];
    const int    qa1 = int( (uint64_t(c.P1[cp]      - c.P1[A])   * M) >> 32 );
    const int    qa2 = int( (uint64_t(c.P2[i2s + n] - c.P2[i2s]) * M) >> 32 );
    const double E   = quickExpected(c.kt->overlapTab, qa1, qa2);
    const double r   = c.kt->randMatch;
    if( (s - s0) < quickThreshold(c.kt->overlapTab, qa1, qa2, n, r, c.kt->quickFraction) ){
	const double nvl = double(nv - n);
	s += nvl * E;
	s += nvl * r;
	s += double(k - Ev) * r;
	return true;
    }
    return false;
}


//! Same value as MergeTrimReads::detectAdapter() for a single offset
static inline double adapterOne(const KernelTables & kt,
				const char * read, const char * qual, int len,
				const char * adapter, int adapterLen,
				int off){
    const double * tabMis = kt.adapterTab;
    const double * tabMat = kt.adapterTab + 128;

    double   s        = 0.0;
    double   qualSum  = 0.0;
    unsigned maxIter  = (unsigned) min(min(size_t(len-off), size_t(adapterLen)), kt.maxAdapterComp);
    unsigned i        = 0;
    unsigned indexRead= off;

    while(i < maxIter){
	indexRead = i + off;
	const unsigned q = qidx(qual[indexRead]);
	if(read[indexRead] == adapter[i] || adapter[i] == 'I'){
	    s += tabMat[q];
	}else{
	    s += tabMis[q];
	}
	qualSum += qual[indexRead];
	i++;

	if(i == 10){
	    if( (s - (i*kt.randMatch)) < kt.fastModeThreshold ){
		double iterationLeft = maxIter - i;
		int qavg = int(qualSum/10);
		s += iterationLeft*0.75*tabMis[qavg & 127] + iterationLeft*0.25*tabMat[qavg & 127];
		indexRead = (maxIter-1) + off;
		break;
	    }
	}
    }

    int extraBases = max(0, len - int(indexRead) - 1);
    s += double(extraBases) * kt.randMatch;
    return s;
}


/////////////////////////////////////////////////////////////////////
//                    portable versions                            //
/////////////////////////////////////////////////////////////////////

static void adapterScanGeneric(const KernelTables & kt,
			       const char * read,
			       const char * qual,
			       int          len,
			       const char * adapter,
			       int          adapterLen,
			       int          count,
			       double *     ll){

    if(count <= 0)
	return;

    const int M0 = (int) min(size_t(adapterLen), kt.maxAdapterComp);

    if(M0 < 1){   //nothing to compare, the generic code path does not apply
	for(int off=0; off<count; off++)
	    ll[off] += adapterOne(kt, read, qual, len, adapter, adapterLen, off);
	return;
    }

    //padded copies, the vector loop may look beyond the read for inactive lanes
    static thread_local Padded pr, pq;
    pr.assign(read, len);
    pq.assign(qual, len);
    const char * rd = pr.p;
    const char * ql = pq.p;

    const double * tabMis = kt.adapterTab;
    const double * tabMat = kt.adapterTab + 128;
    const double   r      = kt.randMatch;

    int off = 0;
    for( ; off + G <= count ; off += G){

	int    M[G];
	for(int l=0;l<G;l++)
	    M[l] = min(len - (off+l), M0);
	const int jmax = M[0];       //the maximum, M is not increasing with the lane

	double s[G], qs[G];
	bool   hop[G];
	for(int l=0;l<G;l++){ s[l]=0.0; qs[l]=0.0; hop[l]=false; }

	const int jfast = min(jmax, 10);

	for(int j=0;j<jfast;j++){
	    const char a = adapter[j];
	    for(int l=0;l<G;l++){
		const int  pos = off + l + j;
		const unsigned q = qidx(ql[pos]);
		const bool eq  = (rd[pos] == a) | (a == 'I');
		const bool act = (j < M[l]);
		const double t = tabMis[q + (eq ? 128 : 0)];
		s[l]   += act ? t : 0.0;
		qs[l]  += act ? double(ql[pos]) : 0.0;
	    }
	}

	if(jmax >= 10){
	    const double thr = kt.fastModeThreshold;
	    const double ref = 10*r;
	    bool anyAlive = false;
	    for(int l=0;l<G;l++){
		if(M[l] >= 10){
		    if( (s[l] - ref) < thr ){
			double iterationLeft = M[l] - 10;
			int qavg = int(qs[l]/10);
			s[l] += iterationLeft*0.75*tabMis[qavg & 127] + iterationLeft*0.25*tabMat[qavg & 127];
			hop[l] = true;
		    }else{
			anyAlive = true;
		    }
		}
	    }

	    if(anyAlive){
		for(int j=10;j<jmax;j++){
		    const char a = adapter[j];
		    for(int l=0;l<G;l++){
			const int  pos = off + l + j;
			const unsigned q = qidx(ql[pos]);
			const bool eq  = (rd[pos] == a) | (a == 'I');
			const bool act = (j < M[l]) & !hop[l];
			const double t = tabMis[q + (eq ? 128 : 0)];
			s[l]   += act ? t : 0.0;
		    }
		}
	    }
	}

	for(int l=0;l<G;l++){
	    const int extraBases = len - (off+l) - M[l];     // = max(0, len - indexRead - 1)
	    s[l] += double(extraBases) * r;
	    ll[off+l] += s[l];
	}
    }

    //remaining offsets
    for( ; off < count ; off++)
	ll[off] += adapterOne(kt, read, qual, len, adapter, adapterLen, off);
}



static void overlapScanGeneric(const KernelTables & kt,
			       const char * read1,
			       const char * qual1,
			       int          len1,
			       const char * read2rev,
			       const char * qual2rev,
			       int          len2,
			       int          count,
			       double *     ll){

    if(count <= 0)
	return;

    static thread_local Padded pb1, pq1, pb2, pq2;
    pb1.assign(read1,    len1);
    pq1.assign(qual1,    len1);
    pb2.assign(read2rev, len2);
    pq2.assign(qual2rev, len2);
    const char * b1 = pb1.p;
    const char * q1 = pq1.p;
    const char * b2 = pb2.p;
    const char * q2 = pq2.p;

    const double * tab = kt.overlapTab;
    const double   r   = kt.randMatch;
    const double * rps = randPrefix(r, size_t(count) + G);

    //quick mode
    const int N = kt.quickBases;
    static thread_local vector<int> vP1, vP2;
    QuickCtx qc;
    if(N > 0){
	quickPrefix(qual1,    len1, vP1);
	quickPrefix(qual2rev, len2, vP2);
	quickInit(qc, kt, len1, len2, vP1, vP2);
    }

    int k0 = 0;
    for( ; k0 < count ; k0 += G){

	int    K[G], A[G];
	for(int l=0;l<G;l++){
	    K[l] = k0 + l;
	    A[l] = max(0, K[l] - len2);    // positions where the second read has not started yet
	}
	const int kmax = k0 + G - 1;       // largest K in the group
	const int Amin = A[0];

	double s[G], s0[G];
	bool   hop[G], appl[G];
	const int cp = (N > 0) ? quickCheckpoint(qc, k0) : 0;
	bool   anyAppl = false;
	for(int l=0;l<G;l++){
	    //the positions before the start of read2 only add likeRandomMatch, one at a time starting from 0.0
	    s[l]   = rps[A[l]];
	    s0[l]  = s[l];
	    hop[l] = false;
	    appl[l]= (N > 0) && quickApplies(qc, K[l], cp);
	    anyAppl= anyAppl || appl[l];
	}

	//both reads can have a base: the full computation
	const int jend1 = min(kmax, len1);
	for(int j=Amin; j<jend1; j++){
	    const unsigned qa = qidx(q1[j]);
	    const char     ba = b1[j];
	    const double * row = tab + size_t(qa)*256;
	    for(int l=0;l<G;l++){
		const int      i2  = len2 - K[l] + j;
		const bool     act = (j >= A[l]) & (j < K[l]) & !hop[l];
		const bool     eq  = (ba == b2[i2]);
		const double   t   = row[(eq ? 128 : 0) + qidx(q2[i2])];
		s[l]   += act ? t : 0.0;
		s[l]   += act ? r : 0.0;
	    }

	    if(anyAppl && j + 1 == cp){
		bool allDone = true;
		for(int l=0;l<G;l++){
		    if(appl[l] && quickCheck(qc, K[l], cp, s0[l], s[l]))
			hop[l] = true;
		    if(!hop[l] && (k0 + l) < count)
			allDone = false;
		}
		if(allDone)
		    break;
	    }
	}

	//read1 is exhausted, only likeRandomMatch is added
	for(int j=max(Amin,len1); j<kmax; j++){
	    for(int l=0;l<G;l++){
		const bool act = (j >= A[l]) & (j < K[l]) & !hop[l];
		s[l] += act ? r : 0.0;
	    }
	}

	const int nl = min(G, count - k0);
	for(int l=0;l<nl;l++)
	    ll[k0+l] += s[l];
    }
}


void ChimeraLanes::build(const vector<string> & chimeras){
    nlanes = int(2*chimeras.size());
    stride = (nlanes + 7) & ~7;
    steps  = 0;
    for(size_t c=0;c<chimeras.size();c++)
	steps = max(steps, int(chimeras[c].size()));

    chars.assign(size_t(steps)*stride + 8, 0);
    length.assign(stride + 8, 0);

    for(size_t c=0;c<chimeras.size();c++){
	for(int o=0;o<2;o++){
	    const int lane = int(2*c) + o;
	    const int n    = int(chimeras[c].size());
	    length[lane]   = (n >= o) ? (n - o) : 0;
	    for(int i=0;i+o<n;i++)
		chars[size_t(i)*stride + lane] = (unsigned char) chimeras[c][i+o];
	}
    }
}


static void chimeraScanGeneric(const KernelTables & kt,
			       const ChimeraLanes & cl,
			       const char * read,
			       const char * qual,
			       int          len,
			       double *     out){
    const double * tabMis = kt.adapterTab;
    const double * tabMat = kt.adapterTab + 128;

    for(int lane=0; lane<cl.nlanes; lane++){
	const int maxidx = min(len, cl.length[lane]);
	if(maxidx <= 0){
	    out[lane] = -DBL_MAX;     // as detectChimera()
	    continue;
	}
	double s = 0.0;
	for(int i=0;i<maxidx;i++){
	    const unsigned char c = cl.chars[size_t(i)*cl.stride + lane];
	    const unsigned q = qidx(qual[i]);
	    if( (unsigned char)read[i] == c || c == 'I' ){
		s += tabMat[q];
	    }else{
		s += tabMis[q];
	    }
	}
	s += double(len - maxidx) * kt.randMatch;
	out[lane] = s;
    }
}


/////////////////////////////////////////////////////////////////////
//                           AVX2                                  //
/////////////////////////////////////////////////////////////////////

#if defined(__x86_64__) && defined(__GNUC__)
#define LEEHOM_HAVE_AVX2 1
}
#include <immintrin.h>
namespace leehom {

static inline int load32(const char * p){ int v; memcpy(&v, p, 4); return v; }

//! 4 lanes of 32 bits all ones/zeros -> 4 lanes of 64 bits all ones/zeros, as a double vector usable as a mask
#define MASK64(x) _mm256_castsi256_pd(_mm256_cvtepi32_epi64(x))

/*
  The adapter scan: candidates (offsets) in the lanes, the steps (positions of the adapter) in the loop.

  Everything that depends only on a position of the read is computed once per read: the log-likelihood
  of a match and of a mismatch of the base at that position (tm, tx) and the running sum of the qualities.
  The loop then only needs contiguous loads.
*/
static const int NVA  = 4;           // vectors of 4 candidates processed together
static const int GGA  = 4*NVA;

__attribute__((target("avx2")))
static void adapterScanAVX2(const KernelTables & kt,
			    const char * read,
			    const char * qual,
			    int          len,
			    const char * adapter,
			    int          adapterLen,
			    int          count,
			    double *     ll){
    if(count <= 0)
	return;

    const int M0 = (int) min(size_t(adapterLen), kt.maxAdapterComp);
    if(M0 < 1){
	adapterScanGeneric(kt, read, qual, len, adapter, adapterLen, count, ll);
	return;
    }

    const double * tabMis = kt.adapterTab;
    const double * tabMat = kt.adapterTab + 128;
    const double   r      = kt.randMatch;

    static thread_local Padded pr;
    static thread_local vector<double> vtm, vtx, vqp;
    pr.assign(read, len);
    vtm.assign(size_t(len) + 2*PAD, 0.0);
    vtx.assign(size_t(len) + 2*PAD, 0.0);
    vqp.assign(size_t(len) + 1, 0.0);
    for(int p=0;p<len;p++){
	const unsigned q = qidx(qual[p]);
	vtm[PAD+p] = tabMis[q];
	vtx[PAD+p] = tabMat[q];
	vqp[p+1]   = vqp[p] + double(qual[p]);       //sums of integers: exact
    }
    const char   * rd = pr.p;
    const double * tm = &vtm[PAD];
    const double * tx = &vtx[PAD];
    const double * qp = &vqp[0];

    int off = 0;
    for( ; off + GGA <= count ; off += GGA){

	alignas(16) int M[GGA];
	for(int l=0;l<GGA;l++)
	    M[l] = min(len - (off+l), M0);
	const int jmax = M[0];       //the maximum, M is not increasing with the lane

	__m256d s[NVA];
	__m128i Mv[NVA];
	for(int v=0;v<NVA;v++){
	    s[v]  = _mm256_setzero_pd();
	    Mv[v] = _mm_load_si128((const __m128i*)(M+4*v));
	}

	const int jfast = min(jmax, 10);
	for(int j=0;j<jfast;j++){
	    const char    a   = adapter[j];
	    const __m128i av  = _mm_set1_epi32((unsigned char)a);
	    const __m128i isI = _mm_set1_epi32(a == 'I' ? -1 : 0);
	    const __m128i jv  = _mm_set1_epi32(j);
	    for(int v=0;v<NVA;v++){
		const int pos = off + 4*v + j;
		const __m128i rb  = _mm_cvtepu8_epi32(_mm_cvtsi32_si128(load32(rd+pos)));
		const __m128i eq  = _mm_or_si128(_mm_cmpeq_epi32(rb, av), isI);
		const __m256d t   = _mm256_blendv_pd(_mm256_loadu_pd(tm+pos), _mm256_loadu_pd(tx+pos), MASK64(eq));
		const __m128i act = _mm_cmpgt_epi32(Mv[v], jv);
		s[v] = _mm256_add_pd(s[v], _mm256_and_pd(MASK64(act), t));
	    }
	}

	alignas(32) double sArr[GGA];
	for(int v=0;v<NVA;v++)
	    _mm256_store_pd(sArr+4*v, s[v]);

	if(jmax >= 10){
	    alignas(16) int hopArr[GGA];
	    const double thr = kt.fastModeThreshold;
	    const double ref = 10*r;
	    bool anyAlive = false;
	    for(int l=0;l<GGA;l++){
		hopArr[l] = 0;
		if(M[l] >= 10){
		    if( (sArr[l] - ref) < thr ){
			double iterationLeft = M[l] - 10;
			int qavg = int((qp[off+l+10] - qp[off+l])/10);
			sArr[l] += iterationLeft*0.75*tabMis[qavg & 127] + iterationLeft*0.25*tabMat[qavg & 127];
			hopArr[l] = -1;
		    }else{
			anyAlive = true;
		    }
		}
	    }

	    for(int v=0;v<NVA;v++)
		s[v] = _mm256_load_pd(sArr+4*v);

	    if(anyAlive){
		__m128i hopv[NVA];
		for(int v=0;v<NVA;v++)
		    hopv[v] = _mm_load_si128((const __m128i*)(hopArr+4*v));

		for(int j=10;j<jmax;j++){
		    const char    a   = adapter[j];
		    const __m128i av  = _mm_set1_epi32((unsigned char)a);
		    const __m128i isI = _mm_set1_epi32(a == 'I' ? -1 : 0);
		    const __m128i jv  = _mm_set1_epi32(j);
		    for(int v=0;v<NVA;v++){
			const int pos = off + 4*v + j;
			const __m128i rb  = _mm_cvtepu8_epi32(_mm_cvtsi32_si128(load32(rd+pos)));
			const __m128i eq  = _mm_or_si128(_mm_cmpeq_epi32(rb, av), isI);
			const __m256d t   = _mm256_blendv_pd(_mm256_loadu_pd(tm+pos), _mm256_loadu_pd(tx+pos), MASK64(eq));
			const __m128i act = _mm_andnot_si128(hopv[v], _mm_cmpgt_epi32(Mv[v], jv));
			s[v] = _mm256_add_pd(s[v], _mm256_and_pd(MASK64(act), t));
		    }
		}
	    }
	}

	for(int v=0;v<NVA;v++)
	    _mm256_store_pd(sArr+4*v, s[v]);
	for(int l=0;l<GGA;l++){
	    const int extraBases = len - (off+l) - M[l];
	    sArr[l] += double(extraBases) * r;
	    ll[off+l] += sArr[l];
	}
    }

    //remaining offsets
    for( ; off < count ; off++)
	ll[off] += adapterOne(kt, read, qual, len, adapter, adapterLen, off);
}


/*
  The overlap scan: candidates (positions of the adapter) in the lanes, the positions of read 1 in the loop.
  Read 2 is stored in the order in which it was sequenced so that the bases that the candidates of a group
  compare with a given base of read 1 are consecutive in memory.
*/
static const int NV  = 2;            // vectors of 4 candidates processed together
static const int GG  = 4*NV;

//! One step of the loop of overlapScanAVX2 for the vectors of a group
template<bool MASKED>
__attribute__((target("avx2"), always_inline))
static inline void overlapStep(int j, int k0,
			       const char * b1, const char * q1, const char * rb2, const char * rq2,
			       const double * tab, __m256d rv,
			       const __m128i * Av, const __m128i * Kv, __m256d * s){
    const double * row = tab + size_t(qidx(q1[j]))*256;
    const __m128i  bav = _mm_set1_epi32((unsigned char)b1[j]);
    const __m128i  c127= _mm_set1_epi32(127);
    const __m128i  c128= _mm_set1_epi32(128);
    const __m128i  jv  = _mm_set1_epi32(j);
    for(int v=0;v<NV;v++){
	const int m0 = k0 - 1 - j + 4*v;
	const __m128i b2v = _mm_cvtepu8_epi32(_mm_cvtsi32_si128(load32(rb2+m0)));
	const __m128i q2v = _mm_and_si128(_mm_cvtepu8_epi32(_mm_cvtsi32_si128(load32(rq2+m0))), c127);
	const __m128i eq  = _mm_cmpeq_epi32(b2v, bav);
	const __m128i idx = _mm_add_epi32(q2v, _mm_and_si128(eq, c128));
	__m256d t = _mm256_i32gather_pd(row, idx, 8);
	if(MASKED){
	    const __m128i act = _mm_andnot_si128(_mm_cmpgt_epi32(Av[v], jv), _mm_cmpgt_epi32(Kv[v], jv));
	    const __m256d am  = MASK64(act);
	    t = _mm256_and_pd(am, t);
	    s[v] = _mm256_add_pd(s[v], t);
	    s[v] = _mm256_add_pd(s[v], _mm256_and_pd(am, rv));
	}else{
	    s[v] = _mm256_add_pd(s[v], t);
	    s[v] = _mm256_add_pd(s[v], rv);
	}
    }
}


//! The same step with the candidates that were given up (Hv) masked out
__attribute__((target("avx2"), always_inline))
static inline void overlapStepQ(int j, int k0,
				const char * b1, const char * q1, const char * rb2, const char * rq2,
				const double * tab, __m256d rv,
				const __m128i * Av, const __m128i * Kv, const __m128i * Hv, __m256d * s){
    const double * row = tab + size_t(qidx(q1[j]))*256;
    const __m128i  bav = _mm_set1_epi32((unsigned char)b1[j]);
    const __m128i  c127= _mm_set1_epi32(127);
    const __m128i  c128= _mm_set1_epi32(128);
    const __m128i  jv  = _mm_set1_epi32(j);
    for(int v=0;v<NV;v++){
	const int m0 = k0 - 1 - j + 4*v;
	const __m128i b2v = _mm_cvtepu8_epi32(_mm_cvtsi32_si128(load32(rb2+m0)));
	const __m128i q2v = _mm_and_si128(_mm_cvtepu8_epi32(_mm_cvtsi32_si128(load32(rq2+m0))), c127);
	const __m128i eq  = _mm_cmpeq_epi32(b2v, bav);
	const __m128i idx = _mm_add_epi32(q2v, _mm_and_si128(eq, c128));
	__m256d t = _mm256_i32gather_pd(row, idx, 8);
	const __m128i act = _mm_andnot_si128(Hv[v], _mm_andnot_si128(_mm_cmpgt_epi32(Av[v], jv), _mm_cmpgt_epi32(Kv[v], jv)));
	const __m256d am  = MASK64(act);
	t = _mm256_and_pd(am, t);
	s[v] = _mm256_add_pd(s[v], t);
	s[v] = _mm256_add_pd(s[v], _mm256_and_pd(am, rv));
    }
}


__attribute__((target("avx2")))
static void overlapScanAVX2(const KernelTables & kt,
			    const char * read1,
			    const char * qual1,
			    int          len1,
			    const char * read2rev,
			    const char * qual2rev,
			    int          len2,
			    int          count,
			    double *     ll){
    if(count <= 0)
	return;

    static thread_local Padded pb1, pq1, prb2, prq2;
    static thread_local vector<char> tmp;
    pb1.assign(read1, len1);
    pq1.assign(qual1, len1);
    //read 2 is stored reversed (i.e. as it was sequenced) so that consecutive candidates are consecutive in memory
    tmp.resize(len2 > 0 ? len2 : 1);
    for(int m=0;m<len2;m++) tmp[m] = read2rev[len2-1-m];
    prb2.assign(&tmp[0], len2);
    for(int m=0;m<len2;m++) tmp[m] = qual2rev[len2-1-m];
    prq2.assign(&tmp[0], len2);

    const char * b1  = pb1.p;
    const char * q1  = pq1.p;
    const char * rb2 = prb2.p;
    const char * rq2 = prq2.p;

    const double * tab = kt.overlapTab;
    const double   r   = kt.randMatch;
    const double * rps = randPrefix(r, size_t(count) + GG);
    const __m256d  rv  = _mm256_set1_pd(r);

    //quick mode
    const int N = kt.quickBases;
    static thread_local vector<int> vP1, vP2;
    QuickCtx qc;
    if(N > 0){
	quickPrefix(qual1,    len1, vP1);
	quickPrefix(qual2rev, len2, vP2);
	quickInit(qc, kt, len1, len2, vP1, vP2);
    }

    for(int k0 = 0; k0 < count; k0 += GG){

	alignas(32) int    A[GG], K[GG];
	alignas(32) double sArr[GG];
	for(int l=0;l<GG;l++){
	    K[l]    = k0 + l;
	    A[l]    = max(0, K[l] - len2);
	    sArr[l] = rps[A[l]];
	}
	const int kmax = k0 + GG - 1;
	const int Amin = A[0];

	__m256d s[NV];
	__m128i Av[NV], Kv[NV];
	for(int v=0;v<NV;v++){
	    s[v]   = _mm256_load_pd(sArr+4*v);
	    Av[v]  = _mm_load_si128((const __m128i*)(A+4*v));
	    Kv[v]  = _mm_load_si128((const __m128i*)(K+4*v));
	}

	const int jend1 = min(kmax, len1);

	if(N <= 0){
	    //the candidates are all active for the steps in [A[GG-1], K[0]), the masks are only needed outside
	    const int a     = min(max(A[GG-1], Amin), jend1);
	    const int b     = min(max(K[0],    a),    jend1);

	    for(int j=Amin; j<a; j++)    overlapStep<true >(j, k0, b1, q1, rb2, rq2, tab, rv, Av, Kv, s);
	    for(int j=a;    j<b; j++)    overlapStep<false>(j, k0, b1, q1, rb2, rq2, tab, rv, Av, Kv, s);
	    for(int j=b; j<jend1; j++)   overlapStep<true >(j, k0, b1, q1, rb2, rq2, tab, rv, Av, Kv, s);

	    //read1 is exhausted, only likeRandomMatch is added
	    for(int j=max(Amin,len1); j<kmax; j++){
		const __m128i jv = _mm_set1_epi32(j);
		for(int v=0;v<NV;v++){
		    const __m128i act = _mm_andnot_si128(_mm_cmpgt_epi32(Av[v], jv), _mm_cmpgt_epi32(Kv[v], jv));
		    s[v] = _mm256_add_pd(s[v], _mm256_and_pd(MASK64(act), rv));
		}
	    }
	}else{
	    //quick mode: all the candidates of the group are looked at together, once they have compared their first bases
	    alignas(32) int    Hm[GG];
	    double s0[GG];
	    bool   appl[GG];
	    const int cp = quickCheckpoint(qc, k0);
	    bool   anyAppl = false;
	    for(int l=0;l<GG;l++){
		Hm[l]   = 0;
		s0[l]   = sArr[l];
		appl[l] = quickApplies(qc, K[l], cp);
		anyAppl = anyAppl || appl[l];
	    }
	    __m128i Hv[NV];
	    for(int v=0;v<NV;v++) Hv[v] = _mm_setzero_si128();

	    bool allDone = false;

	    //The steps in [jFrom,jTo) of the group: the cheap ones while all the candidates are active, which is
	    //the case in [aa,bb), and the ones that mask the candidates that were given up once there are some
	    const int aa = min(max(A[GG-1], Amin), jend1);
	    const int bb = min(max(K[0],    aa),   jend1);
	    #define QUICK_RANGE(jFrom, jTo, WITHHOP)                                                               \
		{                                                                                                   \
		    const int jf_ = (jFrom), jt_ = (jTo);                                                           \
		    if(WITHHOP){                                                                                    \
			for(int j=jf_; j<jt_; j++) overlapStepQ(j, k0, b1, q1, rb2, rq2, tab, rv, Av, Kv, Hv, s);   \
		    }else{                                                                                          \
			for(int j=jf_; j<min(jt_,aa); j++)                                                          \
			    overlapStep<true >(j, k0, b1, q1, rb2, rq2, tab, rv, Av, Kv, s);                        \
			for(int j=max(jf_,aa); j<min(jt_,bb); j++)                                                  \
			    overlapStep<false>(j, k0, b1, q1, rb2, rq2, tab, rv, Av, Kv, s);                        \
			for(int j=max(jf_,bb); j<jt_; j++)                                                          \
			    overlapStep<true >(j, k0, b1, q1, rb2, rq2, tab, rv, Av, Kv, s);                        \
		    }                                                                                               \
		}

	    bool anyHop = false;
	    if(anyAppl){
		QUICK_RANGE(Amin, cp, false);

		for(int v=0;v<NV;v++) _mm256_store_pd(sArr+4*v, s[v]);
		for(int l=0;l<GG;l++)
		    if(appl[l] && quickCheck(qc, K[l], cp, s0[l], sArr[l])){
			Hm[l] = -1;
			anyHop = true;
		    }
		for(int v=0;v<NV;v++){
		    s[v]  = _mm256_load_pd(sArr+4*v);
		    Hv[v] = _mm_load_si128((const __m128i*)(Hm+4*v));
		}
		allDone = true;
		for(int l=0;l<GG;l++)
		    if(!Hm[l] && (k0 + l) < count) allDone = false;

		if(!allDone)
		    QUICK_RANGE(cp, jend1, anyHop);
	    }else{
		QUICK_RANGE(Amin, jend1, false);
	    }
	    #undef QUICK_RANGE

	    //read1 is exhausted, only likeRandomMatch is added
	    if(!allDone){
		for(int j=max(Amin,len1); j<kmax; j++){
		    const __m128i jv = _mm_set1_epi32(j);
		    for(int v=0;v<NV;v++){
			const __m128i act = _mm_andnot_si128(Hv[v], _mm_andnot_si128(_mm_cmpgt_epi32(Av[v], jv), _mm_cmpgt_epi32(Kv[v], jv)));
			s[v] = _mm256_add_pd(s[v], _mm256_and_pd(MASK64(act), rv));
		    }
		}
	    }
	}

	for(int v=0;v<NV;v++)
	    _mm256_store_pd(sArr+4*v, s[v]);
	const int nl = min(GG, count - k0);
	for(int l=0;l<nl;l++)
	    ll[k0+l] += sArr[l];
    }
}


/*
  The chimera scan: the sequences (with the 2 offsets) in the lanes, the positions in the loop.
  At a given position all the lanes compare the same base of the read to their own base, so there
  is no table lookup, only 2 values to choose from.
*/
__attribute__((target("avx2")))
static void chimeraScanAVX2(const KernelTables & kt,
			    const ChimeraLanes & cl,
			    const char * read,
			    const char * qual,
			    int          len,
			    double *     out){
    const int ng = cl.stride / 8;
    if(cl.nlanes == 0)
	return;
    if(ng > 4){
	chimeraScanGeneric(kt, cl, read, qual, len, out);
	return;
    }

    const double * tabMis = kt.adapterTab;
    const double * tabMat = kt.adapterTab + 128;

    alignas(32) int maxidx[32];
    int imax = 0;
    for(int l=0;l<32;l++){
	maxidx[l] = (l < cl.stride) ? min(len, cl.length[l]) : 0;
	imax = max(imax, maxidx[l]);
    }

    __m256d s[4][2];
    for(int g=0;g<4;g++){ s[g][0] = _mm256_setzero_pd(); s[g][1] = _mm256_setzero_pd(); }

    const __m128i cI = _mm_set1_epi8('I');

    for(int i=0;i<imax;i++){
	const unsigned q  = qidx(qual[i]);
	const __m256d  M  = _mm256_set1_pd(tabMat[q]);
	const __m256d  X  = _mm256_set1_pd(tabMis[q]);
	const __m128i  rv = _mm_set1_epi8(read[i]);
	const __m256i  iv = _mm256_set1_epi32(i);
	const unsigned char * ch = &cl.chars[size_t(i)*cl.stride];

	for(int g=0;g<ng;g++){
	    const __m128i c8   = _mm_loadl_epi64((const __m128i*)(ch + 8*g));
	    const __m128i eq8  = _mm_or_si128(_mm_cmpeq_epi8(c8, rv), _mm_cmpeq_epi8(c8, cI));
	    const __m256i eq32 = _mm256_cvtepi8_epi32(eq8);
	    const __m256i act  = _mm256_cmpgt_epi32(_mm256_load_si256((const __m256i*)(maxidx+8*g)), iv);
	    const __m256d eqlo = MASK64(_mm256_castsi256_si128(eq32));
	    const __m256d eqhi = MASK64(_mm256_extracti128_si256(eq32, 1));
	    const __m256d aclo = MASK64(_mm256_castsi256_si128(act));
	    const __m256d achi = MASK64(_mm256_extracti128_si256(act, 1));
	    s[g][0] = _mm256_add_pd(s[g][0], _mm256_and_pd(aclo, _mm256_blendv_pd(X, M, eqlo)));
	    s[g][1] = _mm256_add_pd(s[g][1], _mm256_and_pd(achi, _mm256_blendv_pd(X, M, eqhi)));
	}
    }

    alignas(32) double sArr[32];
    for(int g=0;g<ng;g++){
	_mm256_store_pd(sArr+8*g,   s[g][0]);
	_mm256_store_pd(sArr+8*g+4, s[g][1]);
    }
    for(int lane=0; lane<cl.nlanes; lane++){
	if(maxidx[lane] <= 0){
	    out[lane] = -DBL_MAX;     // as detectChimera()
	}else{
	    out[lane] = sArr[lane] + double(len - maxidx[lane]) * kt.randMatch;
	}
    }
}
#endif


static bool useAVX2(){
#ifdef LEEHOM_HAVE_AVX2
    static const bool b = __builtin_cpu_supports("avx2") && getenv("LEEHOM_NOAVX2") == NULL;
    return b;
#else
    return false;
#endif
}


void adapterScan(const KernelTables & kt, const char * read, const char * qual, int len,
		 const char * adapter, int adapterLen, int count, double * ll){
#ifdef LEEHOM_HAVE_AVX2
    if(useAVX2()){
	adapterScanAVX2(kt, read, qual, len, adapter, adapterLen, count, ll);
	return;
    }
#endif
    adapterScanGeneric(kt, read, qual, len, adapter, adapterLen, count, ll);
}

void overlapScan(const KernelTables & kt, const char * read1, const char * qual1, int len1,
		 const char * read2rev, const char * qual2rev, int len2, int count, double * ll){
#ifdef LEEHOM_HAVE_AVX2
    if(useAVX2()){
	overlapScanAVX2(kt, read1, qual1, len1, read2rev, qual2rev, len2, count, ll);
	return;
    }
#endif
    overlapScanGeneric(kt, read1, qual1, len1, read2rev, qual2rev, len2, count, ll);
}

void chimeraScan(const KernelTables & kt, const ChimeraLanes & cl, const char * read, const char * qual, int len, double * out){
#ifdef LEEHOM_HAVE_AVX2
    if(useAVX2()){
	chimeraScanAVX2(kt, cl, read, qual, len, out);
	return;
    }
#endif
    chimeraScanGeneric(kt, cl, read, qual, len, out);
}

}
