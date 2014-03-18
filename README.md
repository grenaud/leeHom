==========================================================
  Maximum likelihood reconstruction of ancient DNA fragments
==========================================================

QUESTIONS :
   gabriel [dot] reno [ at sign ] gmail.com


About
----------------------

leeHom is a program for the maximum likelihood reconstruction of ancient DNA (the binary is mergeTrimReadsBAM)


Downloading:
----------------------

Go to https://github.com/grenaud/leeHom and either:

1) Download ZIP 

or

2) Do a "git clone --recursive https://github.com/grenaud/leeHom.git"


Installation:
----------------------

Make sure you have cmake installed, if not install it by 
typing "apt-get install cmake" for Ubuntu for example.

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

    src/leeHom

Test data:
----------------------

A small test set of 50,000 ancient DNA sequences from the Altai Neandertal can be found here:

    testData/rawAncientDNA.bam

To launch the reconstruction program for this test data set, use the following command :

    src/leeHom -f AGATCGGAAGAGCACACGTCTGAACTCCAG -s GGAAGAGCGTCGTGTAGGGAAAGAGTGTAG --ancientdna  -o testData/reconsAncientDNA.bam testData/rawAncientDNA.bam

The reconstructed ancient DNA fragments will be in 

    testData/reconsAncientDNA.bam

for fastq data:

    src/leeHom -f AGATCGGAAGAGCACACGTCTGAACTCCAG -s GGAAGAGCGTCGTGTAGGGAAAGAGTGTAG --ancientdna  -fq1 testData/rawAncientDNA.f1.gz -fq2 testData/rawAncientDNA.f2.gz -fqo testData/outfq

The output files will have the  testData/outfq* prefix.
