==========================================================
  Maximum likelihood reconstruction of ancient DNA fragments
==========================================================

QUESTIONS :
   gabriel [dot] reno [ at sign ] gmail.com


About
----------------------

leeHom is programs for the maximum likelihood reconstruction of ancient DNA (the binary is mergeTrimReadsBAM)


Downloading:
----------------------

Go to https://github.com/grenaud/leeHom and either:

1) Download ZIP 

or

2) Do a "git clone --recursive https://github.com/grenaud/leeHom.git"


Installation:
----------------------

1) Build Bamtools first:

    cd bamtools/   
    mkdir build/   
    cd build/
    cmake ..
    make 
    cd ../..

2) Build the submodules and main code by typing :

    make



running the program:
----------------------

To launch the program simply type:

    src/mergeTrimReadsBAM

Test data:
----------------------

A small test set of 50,000 ancient DNA sequences from the Altai Neandertal can be found here:

    testData/rawAncientDNA.bam

To launch the reconstruction program for this test data set, use the following command :

    src/mergeTrimReadsBAM -f AGATCGGAAGAGCACACGTCTGAACTCCAG -s GGAAGAGCGTCGTGTAGGGAAAGAGTGTAG --ancientdna  -o testData/reconsAncientDNA.bam testData/rawAncientDNA.bam

The reconstructed ancient DNA fragments will be in 

    testData/reconsAncientDNA.bam

for fastq data:

    src/mergeTrimReadsBAM -f AGATCGGAAGAGCACACGTCTGAACTCCAG -s GGAAGAGCGTCGTGTAGGGAAAGAGTGTAG --ancientdna  -fq1 testData/rawAncientDNA.f1.gz -fq2 testData/rawAncientDNA.f2.gz -fqo testData/outfq

The output files will have the  testData/outfq* prefix.
