#ifndef leeHom_BamRecord_h
#define leeHom_BamRecord_h

#include <string>
#include <vector>
#include <stdint.h>

#include <htslib/sam.h>

/*!
  A BAM record held as plain C++ members, converted from/to htslib's bam1_t at the edges.

  The names of the members and the semantics of the tag functions follow what the rest of
  leeHom expects: the sequence is ASCII, the qualities are PHRED+33 ASCII (or a string of 0xFF
  when they are not stored) and the tags are kept in their raw BAM encoding in TagData.
*/
class BamRecord {
public:
    std::string           Name;
    std::string           QueryBases;
    std::string           Qualities;
    uint32_t              AlignmentFlag;
    int32_t               RefID;
    int32_t               Position;
    int32_t               MateRefID;
    int32_t               MatePosition;
    int32_t               InsertSize;
    uint8_t               MapQuality;
    std::vector<uint32_t> Cigar;
    std::string           TagData;    //raw BAM auxiliary data

    BamRecord();

    //! Fill this record from an htslib record
    void fromBam(const bam1_t * b);
    //! Rebuild the htslib record from this one
    void toBam(bam1_t * b) const;

    bool IsPaired()      const { return (AlignmentFlag & BAM_FPAIRED)   != 0; }
    bool IsMapped()      const { return (AlignmentFlag & BAM_FUNMAP)    == 0; }
    bool IsFirstMate()   const { return (AlignmentFlag & BAM_FREAD1)    != 0; }
    bool IsSecondMate()  const { return (AlignmentFlag & BAM_FREAD2)    != 0; }
    bool IsFailedQC()    const { return (AlignmentFlag & BAM_FQCFAIL)   != 0; }
    bool IsDuplicate()   const { return (AlignmentFlag & BAM_FDUP)      != 0; }

    void SetIsMapped(bool ok)     { setFlag(BAM_FUNMAP,  !ok); }
    void SetIsFailedQC(bool ok)   { setFlag(BAM_FQCFAIL,  ok); }
    void SetIsDuplicate(bool ok)  { setFlag(BAM_FDUP,     ok); }

    bool HasTag(const std::string & tag) const;
    //! Type character of the tag, false if it is not present
    bool GetTagType(const std::string & tag, char & type) const;
    //! Value of a tag of type Z or H
    bool GetTag(const std::string & tag, std::string & destination) const;
    //! Value of a tag of any integer type
    bool GetTag(const std::string & tag, int & destination) const;
    //! Appends a new tag of type Z (string), fails if the tag is already present
    bool AddTag(const std::string & tag, const std::string & type, const std::string & value);
    //! Appends a new tag of type i (32-bit integer), fails if the tag is already present
    bool AddTag(const std::string & tag, const std::string & type, int32_t value);
    //! Replaces (removes then appends) a string tag
    bool EditTag(const std::string & tag, const std::string & type, const std::string & value);
    void RemoveTag(const std::string & tag);

private:
    void setFlag(uint32_t f, bool on){ if(on){ AlignmentFlag |= f; }else{ AlignmentFlag &= ~f; } }
    //! Offset in TagData of the tag, -1 if absent. The type is returned in type
    int  findTag(const std::string & tag, char * type = NULL) const;
};

#endif
