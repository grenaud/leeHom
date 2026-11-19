#!/bin/bash

RVAL=0
set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

echo -n "Running on BAM format:"
../src/leeHom -f AGATCGGAAGAGCACACGTCTGAACTCCAG -s GGAAGAGCGTCGTGTAGGGAAAGAGTGTAG   -o out.bam ../testData/rawAncientDNA.bam |& md5sum > bam.md5sum
echo -e " ${GREEN}ok${NC}"

echo -n "testing md5sum:"
if diff bam.md5sum bam.md5sum_ > /dev/null
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi


echo -n "Running on paired-end fastq format:"
../src/leeHom  -fqo out -fq1 ../testData/rawAncientDNA.f1.gz -fq2 ../testData/rawAncientDNA.f2.gz |& md5sum > fq.md5sum
echo -e " ${GREEN}ok${NC}"

echo -n "testing md5sum:"
if diff fq.md5sum fq.md5sum_ > /dev/null
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi
zcat out.fq.gz |head -n 60 |md5sum > fq.head60out.md5sum

echo -n "testing output md5sum:"
if diff fq.head60out.md5sum fq.head60out.md5sum_ > /dev/null
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi




echo -n "Running on single-end fastq format:"
../src/leeHom  -fqo out_s -fq1 ../testData/rawAncientDNA.f1.gz |& md5sum > fq_s.md5sum
echo -e " ${GREEN}ok${NC}"

echo -n "testing md5sum:"
if diff fq_s.md5sum fq_s.md5sum_ > /dev/null
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi
zcat out_s.fq.gz |head -n 60 |md5sum > fq_s.head60out.md5sum

echo -n "testing output md5sum:"
if diff fq_s.head60out.md5sum fq_s.head60out.md5sum_ > /dev/null
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi




echo -n "Running ancient DNA mode:"
../src/leeHom --ancientdna -f AGATCGGAAGAGCACACGTCTGAACTCCAG -s GGAAGAGCGTCGTGTAGGGAAAGAGTGTAG   -o out.bam ../testData/rawAncientDNA.bam |& md5sum > bamanc.md5sum
echo -e " ${GREEN}ok${NC}"

echo -n "testing md5sum:"
if diff bamanc.md5sum bamanc.md5sum_ > /dev/null
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi

echo -n "Running discovering adapters mode:"
../src/leeHom --auto --ancientdna -f AGATCGGAAGAGCACACGTCTGAACTCCAG -s GGAAGAGCGTCGTGTAGGGAAAGAGTGTAG   -o out.bam ../testData/rawAncientDNA.bam |& md5sum > bamancadp.md5sum
echo -e " ${GREEN}ok${NC}"

echo -n "testing md5sum:"
if diff bamancadp.md5sum bamancadp.md5sum_ > /dev/null
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi


echo -n "Checking that the vectorized code gives the same likelihoods as the original code:"
fuzzA=$(../src/fuzz_kernels)                   || { echo -e " ${RED}test failed${NC}"; echo "$fuzzA"; exit 1; }
fuzzB=$(LEEHOM_NOAVX2=1 ../src/fuzz_kernels)   || { echo -e " ${RED}test failed${NC}"; echo "$fuzzB"; exit 1; }
fuzzC=$(LEEHOM_REFERENCE=1 ../src/fuzz_kernels) || { echo -e " ${RED}test failed${NC}"; echo "$fuzzC"; exit 1; }
if [ "${fuzzA%%likelihood*}" == "${fuzzB%%likelihood*}" ] && [ "${fuzzA%%likelihood*}" == "${fuzzC%%likelihood*}" ]
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi


# The output must not depend on the number of threads and must be in the order of the input
echo -n "Running paired-end fastq with 1 and with 4 threads:"
../src/leeHom --seed 7 --ancientdna -t 1 -fqo thr1 -fq1 ../testData/rawAncientDNA.f1.gz -fq2 ../testData/rawAncientDNA.f2.gz 2> thr1.log
../src/leeHom --seed 7 --ancientdna -t 4 -z 500 -fqo thr4 -fq1 ../testData/rawAncientDNA.f1.gz -fq2 ../testData/rawAncientDNA.f2.gz 2> thr4.log
echo -e " ${GREEN}ok${NC}"

echo -n "testing that the output is identical (same reads, same order):"
same=1
for suffix in .fq.gz _r1.fq.gz _r2.fq.gz .fail.fq.gz _r1.fail.fq.gz _r2.fail.fq.gz
do
    if [ "$(zcat thr1$suffix | md5sum)" != "$(zcat thr4$suffix | md5sum)" ]; then same=0; fi
done
if [ $same == 1 ] && diff thr1.log thr4.log > /dev/null
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi

echo -n "testing that the fastq files are valid gzip files with the reads in the order of the input:"
zcat ../testData/rawAncientDNA.f1.gz | awk 'NR%4==1 {sub(/ .*/,""); sub(/\/1$/,""); print}' > thr_in_names.txt
# every cluster ends up in exactly one of the files: merged, kept or failed. The names in each file must follow the order of the input
ok=1
for suffix in .fq.gz _r1.fq.gz .fail.fq.gz _r1.fail.fq.gz
do
    zcat thr4$suffix | awk 'NR%4==1 {sub(/ .*/,""); sub(/\/1$/,""); print}' > thr_out_names.txt
    # the names of the output must be a subsequence of the input names, if the order was changed that fails
    if ! awk 'NR==FNR {want[++n]=$0; next} {while (i<n && want[++i]!=$0) {} if (want[i]!=$0) {bad=1; exit}} END {exit bad}' thr_in_names.txt thr_out_names.txt
    then
        ok=0
    fi
done
if [ $ok == 1 ]
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi


# Quick mode: the output must not depend on the number of threads, and must stay close to the exact computation
echo -n "Running the quick mode with 1 and with 4 threads:"
../src/leeHom --seed 7 --ancientdna --quick -t 1 -fqo qk1 -fq1 ../testData/rawAncientDNA.f1.gz -fq2 ../testData/rawAncientDNA.f2.gz 2> qk1.log
../src/leeHom --seed 7 --ancientdna --quick -t 4 -z 500 -fqo qk4 -fq1 ../testData/rawAncientDNA.f1.gz -fq2 ../testData/rawAncientDNA.f2.gz 2> qk4.log
echo -e " ${GREEN}ok${NC}"

echo -n "testing that the output is identical with 1 and 4 threads:"
same=1
for suffix in .fq.gz _r1.fq.gz _r2.fq.gz .fail.fq.gz _r1.fail.fq.gz _r2.fail.fq.gz
do
    if [ "$(zcat qk1$suffix | md5sum)" != "$(zcat qk4$suffix | md5sum)" ]; then same=0; fi
done
if [ $same == 1 ] && diff qk1.log qk4.log > /dev/null
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi

echo -n "testing that the quick mode changes less than 1% of the reads:"
changed=$(diff <(zcat thr1.fq.gz) <(zcat qk1.fq.gz) | grep -c '^<')
total=$(zcat thr1.fq.gz | awk 'NR%4==1' | wc -l)
if [ $((changed * 100)) -lt $total ]
then
    echo -e " ${GREEN}test passed${NC} ($changed of $total merged reads differ)"
else
    echo -e " ${RED}test failed${NC} ($changed of $total merged reads differ)"
    exit 1
fi

# The compression level changes the size of the output, not what it contains
echo -n "testing that --gzlevel does not change the content of the output:"
../src/leeHom --seed 7 --ancientdna --gzlevel 1 -t 2 -fqo gz1 -fq1 ../testData/rawAncientDNA.f1.gz -fq2 ../testData/rawAncientDNA.f2.gz 2> gz1.log
same=1
for suffix in .fq.gz _r1.fq.gz _r2.fq.gz .fail.fq.gz _r1.fail.fq.gz _r2.fail.fq.gz
do
    if [ "$(zcat thr1$suffix | md5sum)" != "$(zcat gz1$suffix | md5sum)" ]; then same=0; fi
done
if [ $same == 1 ] && [ "$(stat -c %s gz1.fq.gz)" -gt "$(stat -c %s thr1.fq.gz)" ]
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi

# Fastq in, BAM out: same reads as the fastq output, in the same order, whatever the number of threads (needs samtools to read the BAM)
if ! command -v samtools > /dev/null
then
    echo "samtools not found, skipping the fastq to BAM tests"
else
echo -n "Running fastq input with BAM output:"
../src/leeHom --seed 7 --ancientdna -fqt fqi.fq.gz -fq1 ../testData/rawAncientDNA.f1.gz -fq2 ../testData/rawAncientDNA.f2.gz 2> fqi.log
../src/leeHom --seed 7 --ancientdna -t 1 -o fq2bam1.bam -fq1 ../testData/rawAncientDNA.f1.gz -fq2 ../testData/rawAncientDNA.f2.gz 2> fq2bam1.log
../src/leeHom --seed 7 --ancientdna -t 4 -z 500 -u -fq1 ../testData/rawAncientDNA.f1.gz -fq2 ../testData/rawAncientDNA.f2.gz > fq2bam4.bam 2> fq2bam4.log
echo -e " ${GREEN}ok${NC}"

echo -n "testing that the BAM has the sequences and qualities of the fastq output, in the same order:"
if [ "$(zcat fqi.fq.gz | awk 'NR%4==2 || NR%4==0' | md5sum)" == "$(samtools view fq2bam1.bam | awk '{print $10; print $11}' | md5sum)" ]
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi

echo -n "testing that the BAM is identical with 1 and 4 threads:"
if [ "$(samtools view fq2bam1.bam | md5sum)" == "$(samtools view fq2bam4.bam | md5sum)" ] && diff fq2bam1.log fq2bam4.log > /dev/null
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi
fi


###############################################################################
# Fastq in, one BAM (or uncompressed BAM) out: --keepOrig and the QC failed flag
# These use src/bamview (built by make test) and not samtools.
###############################################################################
if [ ! -x ../src/bamview ]
then
    echo "src/bamview not found (make -C src bamview), skipping the fastq to BAM, keepOrig and QC fail tests"
else
BV=../src/bamview
ADP="-f AGATCGGAAGAGCACACGTCTGAACTCCAG -s GGAAGAGCGTCGTGTAGGGAAAGAGTGTAG"

# Illumina (Casava >= 1.8) style headers: every 3rd cluster is flagged as filtered (1:Y:), the others pass (1:N:)
for n in 1 2
do
    zcat ../testData/rawAncientDNA.f$n.gz | head -n 4000 | awk -v n=$n 'NR%4==1 {r++; sub(/\/[12]$/,""); $0=$0" "n":"((r%3==0)?"Y":"N")":0:ACGT"} {print}' | gzip > qc.f$n.fq.gz
done
# names of the filtered clusters and the original reads (name sequence quality) of read 1 and read 2
zcat qc.f1.fq.gz | awk 'NR%4==1 {n=substr($1,2); f[n]=($2 ~ /^[12]:Y:/)} END{for(n in f) print n, f[n]}' | sort > qc.filtered.txt
{ zcat qc.f1.fq.gz | awk 'NR%4==1{n=substr($1,2)} NR%4==2{s=$1} NR%4==0{print n"\t"s"\t"$1}';
  zcat qc.f2.fq.gz | awk 'NR%4==1{n=substr($1,2)} NR%4==2{s=$1} NR%4==0{print n"\t"s"\t"$1}'; } | sort > qc.input.txt

echo -n "Running fastq input to BAM and uncompressed BAM, one output file, 1 and 4 threads:"
../src/leeHom --seed 7 --ancientdna $ADP -t 1 -o qc1.bam   -fq1 qc.f1.fq.gz -fq2 qc.f2.fq.gz 2> qc1.log
../src/leeHom --seed 7 --ancientdna $ADP -t 4 -z 100 -u    -fq1 qc.f1.fq.gz -fq2 qc.f2.fq.gz > qc4.bam 2> qc4.log
echo -e " ${GREEN}ok${NC}"

echo -n "testing that the BAM output is identical for 1 and 4 threads, compressed or not:"
if [ "$($BV qc1.bam | md5sum)" == "$($BV qc4.bam | md5sum)" ] && [ "$($BV qc1.bam | wc -l)" -gt 0 ] && diff qc1.log qc4.log > /dev/null
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi

echo -n "testing that -u really gives an uncompressed BAM:"
# stored (level 0) BGZF blocks are larger than the deflated blocks of the same records
if [ "$(stat -c %s qc4.bam)" -gt "$(stat -c %s qc1.bam)" ]
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi

echo -n "testing that the QC failed flag (0x200) of the fastq header (1:Y:) is carried to the BAM, for merged and unmerged reads:"
# every record of a filtered cluster must be QC failed, and the others must not all be (leeHom fails some of them itself)
if $BV qc1.bam | awk 'NR==FNR {filt[$1]=$2; next} {n=$1; sub(/_M$/,"",n); fail=int($2/512)%2; if(filt[n]==1 && !fail) bad++; if(filt[n]==0 && !fail) pass++; if(filt[n]==1 && ($2%2==0)) mergedfilt++} END{exit (bad>0 || pass==0 || mergedfilt==0)}' qc.filtered.txt -
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi

echo -n "Running fastq input to BAM with --keepOrig:"
../src/leeHom --seed 7 --ancientdna $ADP --keepOrig -t 3 -z 100 -o qk.bam -fq1 qc.f1.fq.gz -fq2 qc.f2.fq.gz 2> qk.log
echo -e " ${GREEN}ok${NC}"

echo -n "testing that --keepOrig only adds the original reads (duplicate flag) and the _M suffix to the merged reads:"
# without the originals (0x400) and the suffix, the records are those of the run without --keepOrig
if [ "$($BV qk.bam | awk -F'\t' -v OFS='\t' '{if(int($2/1024)%2==0){sub(/_M$/,"",$1); print}}' | md5sum)" == "$($BV qc1.bam | md5sum)" ] &&
   [ "$($BV qk.bam | awk 'int($2/1024)%2==1' | wc -l)" -gt 0 ]
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi

echo -n "testing that each merged read has both original reads, unchanged, with the QC failed flag of the input:"
ok=1
# 2 originals per merged read (the merged reads are the unpaired records)
merged=$($BV qk.bam | awk '$2%2==0' | wc -l)
origs=$($BV qk.bam | awk 'int($2/1024)%2==1' | wc -l)
[ "$origs" == "$((2*merged))" ] || ok=0
# the originals are the input reads (name, sequence, qualities), mate 1 (0x40) and mate 2 (0x80) are both there
[ "$($BV qk.bam | awk 'int($2/1024)%2==1 && int($2/64)%2==1' | wc -l)" == "$merged" ] || ok=0
[ "$($BV qk.bam | awk 'int($2/1024)%2==1 && int($2/128)%2==1' | wc -l)" == "$merged" ] || ok=0
$BV qk.bam | awk 'int($2/1024)%2==1 {print $1"\t"$3"\t"$4}' | sort > qk.origs.txt
[ "$(comm -23 qk.origs.txt qc.input.txt | wc -l)" == 0 ] || ok=0
# the flags: originals of a filtered cluster are QC failed, the others are not
$BV qk.bam | awk 'NR==FNR {filt[$1]=$2; next} int($2/1024)%2==1 {fail=int($2/512)%2; if(fail != filt[$1]) bad++} END{exit bad>0}' qc.filtered.txt - || ok=0
if [ $ok == 1 ]
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi

echo -n "Running single-end fastq input to BAM and to uncompressed BAM, with and without --keepOrig:"
../src/leeHom --seed 7 --ancientdna $ADP -t 1 -o qs1.bam -fq1 qc.f1.fq.gz 2> qs1.log
../src/leeHom --seed 7 --ancientdna $ADP -t 4 -z 100 -u   -fq1 qc.f1.fq.gz > qs4.bam 2> qs4.log
../src/leeHom --seed 7 --ancientdna $ADP --keepOrig -t 4 -z 100 -o qsk.bam -fq1 qc.f1.fq.gz 2> qsk.log
echo -e " ${GREEN}ok${NC}"

echo -n "testing the single-end BAM: same for 1 and 4 threads, QC failed flag carried, one record per read:"
ok=1
[ "$($BV qs1.bam | md5sum)" == "$($BV qs4.bam | md5sum)" ] || ok=0
[ "$($BV qs1.bam | wc -l)" == "$(zcat qc.f1.fq.gz | awk 'NR%4==1' | wc -l)" ] || ok=0
$BV qs1.bam | awk 'NR==FNR {filt[$1]=$2; next} {fail=int($2/512)%2; if(filt[$1]==1 && !fail) bad++; if(filt[$1]==0 && !fail) pass++} END{exit (bad>0 || pass==0)}' qc.filtered.txt - || ok=0
if [ $ok == 1 ]
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi

echo -n "testing --keepOrig with single-end reads: originals only for the changed reads, which get the _M suffix, QC failed flag carried:"
ok=1
[ "$($BV qsk.bam | awk -F'\t' -v OFS='\t' '{if(int($2/1024)%2==0){sub(/_M$/,"",$1); print}}' | md5sum)" == "$($BV qs1.bam | md5sum)" ] || ok=0
nm=$($BV qsk.bam | awk '$1 ~ /_M$/' | wc -l)
no=$($BV qsk.bam | awk 'int($2/1024)%2==1' | wc -l)
[ "$nm" == "$no" ] && [ "$no" -gt 0 ] || ok=0
# reads that are not marked as changed have no _M, and the originals are the input reads
[ "$($BV qsk.bam | awk '$1 !~ /_M$/ && int($2/1024)%2==0' | wc -l)" == "$(( $($BV qs1.bam | wc -l) - nm ))" ] || ok=0
$BV qsk.bam | awk 'int($2/1024)%2==1 {print $1"\t"$3"\t"$4}' | sort > qsk.origs.txt
[ "$(comm -23 qsk.origs.txt qc.input.txt | wc -l)" == 0 ] || ok=0
$BV qsk.bam | awk 'NR==FNR {filt[$1]=$2; next} int($2/1024)%2==1 {fail=int($2/512)%2; if(fail != filt[$1]) bad++} END{exit bad>0}' qc.filtered.txt - || ok=0
if [ $ok == 1 ]
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi
fi
