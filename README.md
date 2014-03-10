About
----------------------

FlintstoneXL is series of programs for both ancient and modern DNA:

1) Maximum likelihood reconstruction of ancient DNA (mergeTrimReadsBAM)
2) Maximum likelihood demultiplexing of reads (assignRG)
3) Filtering low quality clusters (filterReads)


Downloading:
----------------------

Go to https://github.com/grenaud/FlintstoneXL and either:

1) Download ZIP 

or

2) Do a "git clone --recursive https://github.com/grenaud/FlintstoneXL.git"


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
