#include <cstring>
#include <cstdlib>
#include <iostream>

#include "BamRecord.h"

using namespace std;

BamRecord::BamRecord():
    AlignmentFlag(BAM_FUNMAP),
    RefID(-1),
    Position(-1),
    MateRefID(-1),
    MatePosition(-1),
    InsertSize(0),
    MapQuality(0){
}


void BamRecord::fromBam(const bam1_t * b){
    const int32_t lseq = b->core.l_qseq;

    Name.assign(bam_get_qname(b));

    const uint8_t * s = bam_get_seq(b);
    QueryBases.resize(lseq);
    for(int32_t i=0;i<lseq;i++)
	QueryBases[i] = seq_nt16_str[bam_seqi(s,i)];

    const uint8_t * q = bam_get_qual(b);
    if(lseq > 0 && q[0] == 0xFF){       //qualities are not stored
	Qualities.assign(lseq, char(0xFF));
    }else{
	Qualities.resize(lseq);
	for(int32_t i=0;i<lseq;i++)
	    Qualities[i] = char(q[i] + 33);
    }

    AlignmentFlag = b->core.flag;
    RefID         = b->core.tid;
    Position      = int32_t(b->core.pos);
    MateRefID     = b->core.mtid;
    MatePosition  = int32_t(b->core.mpos);
    InsertSize    = int32_t(b->core.isize);
    MapQuality    = b->core.qual;

    const uint32_t * c = bam_get_cigar(b);
    Cigar.assign(c, c + b->core.n_cigar);

    TagData.assign((const char *)bam_get_aux(b), bam_get_l_aux(b));
}


void BamRecord::toBam(bam1_t * b) const {
    vector<uint8_t> rawQual;
    const uint8_t * qptr = NULL;           //NULL means not stored
    if(!Qualities.empty() && uint8_t(Qualities[0]) != 0xFF){
	rawQual.resize(Qualities.size());
	for(size_t i=0;i<Qualities.size();i++)
	    rawQual[i] = uint8_t(Qualities[i] - 33);
	qptr = &rawQual[0];
    }

    if(bam_set1(b,
		Name.size(), Name.c_str(),
		AlignmentFlag,
		RefID, Position,
		MapQuality,
		Cigar.size(), Cigar.empty() ? NULL : &Cigar[0],
		MateRefID, MatePosition,
		InsertSize,
		QueryBases.size(), QueryBases.c_str(), (const char *)qptr,
		TagData.size()) < 0){
	cerr<<"Unable to build the BAM record "<<Name<<endl;
	exit(1);
    }

    //bam_set1() only reserves the room for the auxiliary data
    if(!TagData.empty()){
	memcpy(b->data + b->l_data, TagData.data(), TagData.size());
	b->l_data += int(TagData.size());
    }
}



//! Pointer past the value of type 'type' that starts at p, NULL if it is truncated or has an unknown type
static const char * skipValue(const char * p, const char * end, char type){
    size_t n = 0;
    switch(type){
    case 'A': case 'c': case 'C':           n = 1; break;
    case 's': case 'S':                     n = 2; break;
    case 'i': case 'I': case 'f':           n = 4; break;
    case 'Z': case 'H': {
	const void * z = memchr(p, 0, end - p);
	if(z == NULL) return NULL;
	return (const char *)z + 1;
    }
    case 'B': {
	if(end - p < 5) return NULL;
	uint32_t count;
	memcpy(&count, p+1, 4);
	size_t elem;
	switch(p[0]){
	case 'c': case 'C':                 elem = 1; break;
	case 's': case 'S':                 elem = 2; break;
	case 'i': case 'I': case 'f':       elem = 4; break;
	default: return NULL;
	}
	n = 5 + size_t(count)*elem;
	break;
    }
    default:
	return NULL;
    }
    if(size_t(end - p) < n) return NULL;
    return p + n;
}


int BamRecord::findTag(const string & tag, char * type) const {
    if(tag.size() != 2)
	return -1;

    const char * start = TagData.data();
    const char * end   = start + TagData.size();
    const char * p     = start;

    while(end - p >= 3){
	const char t = p[2];
	if(p[0] == tag[0] && p[1] == tag[1]){
	    if(type != NULL) *type = t;
	    return int(p - start);
	}
	p = skipValue(p+3, end, t);
	if(p == NULL)
	    return -1;
    }
    return -1;
}


bool BamRecord::HasTag(const string & tag) const {
    return findTag(tag) >= 0;
}


bool BamRecord::GetTagType(const string & tag, char & type) const {
    return findTag(tag, &type) >= 0;
}


bool BamRecord::GetTag(const string & tag, string & destination) const {
    char type;
    const int off = findTag(tag, &type);
    if(off < 0 || (type != 'Z' && type != 'H'))
	return false;
    destination.assign(TagData.c_str() + off + 3);
    return true;
}


bool BamRecord::GetTag(const string & tag, int & destination) const {
    char type;
    const int off = findTag(tag, &type);
    if(off < 0)
	return false;

    const char * v = TagData.data() + off + 3;
    switch(type){
    case 'c': { int8_t   x; memcpy(&x, v, 1); destination = x; return true; }
    case 'C': { uint8_t  x; memcpy(&x, v, 1); destination = x; return true; }
    case 's': { int16_t  x; memcpy(&x, v, 2); destination = x; return true; }
    case 'S': { uint16_t x; memcpy(&x, v, 2); destination = x; return true; }
    case 'i': { int32_t  x; memcpy(&x, v, 4); destination = x; return true; }
    case 'I': { uint32_t x; memcpy(&x, v, 4); destination = int(x); return true; }
    default:
	return false;
    }
}


bool BamRecord::AddTag(const string & tag, const string & type, const string & value){
    if(tag.size() != 2 || type.size() != 1 || (type[0] != 'Z' && type[0] != 'H'))
	return false;
    if(HasTag(tag))
	return false;

    TagData.append(tag);
    TagData.push_back(type[0]);
    TagData.append(value);
    TagData.push_back('\0');
    return true;
}


bool BamRecord::AddTag(const string & tag, const string & type, int32_t value){
    if(tag.size() != 2 || type != "i")
	return false;
    if(HasTag(tag))
	return false;

    TagData.append(tag);
    TagData.push_back('i');
    TagData.append((const char *)&value, 4);
    return true;
}


bool BamRecord::EditTag(const string & tag, const string & type, const string & value){
    if(HasTag(tag))
	RemoveTag(tag);
    return AddTag(tag, type, value);
}


void BamRecord::RemoveTag(const string & tag){
    char type;
    const int off = findTag(tag, &type);
    if(off < 0)
	return;

    const char * start = TagData.data();
    const char * next  = skipValue(start + off + 3, start + TagData.size(), type);
    if(next == NULL)
	return;
    TagData.erase(off, next - (start + off));
}
