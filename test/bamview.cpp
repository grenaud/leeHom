// Minimal BAM dumper for test/test.sh so that the tests do not need samtools.
// Prints one line per record: name, flag, sequence, qualities (ASCII, phred+33) and the tags (sorted by appearance).
#include <cstdio>
#include <cstdint>
#include <string>

#include "htslib/sam.h"

int main(int argc, char ** argv){
    if(argc != 2){
	fprintf(stderr, "usage: %s file.bam\n", argv[0]);
	return 1;
    }
    samFile *   fp  = sam_open(argv[1], "r");
    if(fp == NULL){ fprintf(stderr, "cannot open %s\n", argv[1]); return 1; }
    sam_hdr_t * hdr = sam_hdr_read(fp);
    if(hdr == NULL){ fprintf(stderr, "cannot read the header of %s\n", argv[1]); return 1; }
    bam1_t *    b   = bam_init1();
    while(sam_read1(fp, hdr, b) >= 0){
	std::string seq, qual;
	const uint8_t * s = bam_get_seq(b);
	const uint8_t * q = bam_get_qual(b);
	for(int i=0;i<b->core.l_qseq;i++){
	    seq  += seq_nt16_str[bam_seqi(s,i)];
	    qual += (q[0] == 0xff) ? '*' : char(q[i]+33);
	}
	kstring_t tags = {0,0,NULL};
	std::string tagstr;
	// sam_format1 gives all the fields, keep only the optional ones (after the 11th)
	if(sam_format1(hdr, b, &tags) >= 0){
	    int tabs=0;
	    for(size_t i=0;i<tags.l;i++){
		if(tags.s[i]=='\t') tabs++;
		if(tabs>=11) tagstr += tags.s[i];
	    }
	}
	free(tags.s);
	printf("%s\t%d\t%s\t%s%s\n", bam_get_qname(b), (int) b->core.flag, seq.c_str(), qual.c_str(), tagstr.c_str());
    }
    bam_destroy1(b);
    sam_hdr_destroy(hdr);
    sam_close(fp);
    return 0;
}
