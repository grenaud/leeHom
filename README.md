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
