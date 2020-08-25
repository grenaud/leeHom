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


echo -n "Running on fastq format:"
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
