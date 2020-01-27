
#BAMTOOLS= $(realpath ../bamtools/)
#LIBGAB= $(realpath ../libgab/)


CXX      = g++   #-g -pg
CXXFLAGS+= -O3 -Wall -I../lib/bamtools/include/ -I../lib/bamtools/src/ -I../lib/libgab -I../lib/samtools/ -I../lib/htslib/ -I../lib/libgab/gzstream/  
LDFLAGS  = -lgsl -lcurl -lgslcblas -llzma -lbz2 -lm -lz  -lpthread 
LDLIBS+=  ../lib/bamtools/build/src/utils/CMakeFiles/BamTools-utils.dir/*cpp.o ../lib/bamtools/build/src/api/libbamtools.a ../lib/libgab/utils.o 
LDLIBS+= ../lib/samtools/bedidx.o ../lib/htslib/libhts.a ../lib/samtools/libbam.a ../lib/samtools/libst.a 



all:    ../lib/bamtools/build/src/api/libbamtools.a ../lib/htslib/libhts.a ../lib/samtools/bedidx.o ../lib/libgab/utils.o MergeTrimReads.o leeHom leeHomMulti

%.o: %.cpp
	${CXX} -c ${CXXFLAGS} $^ 

leeHom: leeHom.o  MergeTrimReads.o ../lib/libgab/PutProgramInHeader.o ../lib/libgab/FastQObj.o  ../lib/libgab/FastQParser.o  ../lib/libgab/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)


leeHomMulti: leeHomMulti.o  MergeTrimReads.o ../lib/libgab/PutProgramInHeader.o ../lib/libgab/PutProgramInHeaderHTS.o  ../lib/libgab/FastQObj.o  ../lib/libgab/FastQParser.o ../lib/libgab/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

../lib/htslib/libhts.a:  ../lib/htslib/hts_internal.h
	cd ../lib/htslib/ && make && cd ../..

../lib/htslib/hts_internal.h:
	rm -rfv ../lib/htslib/
	git clone --recursive https://github.com/samtools/htslib.git ../lib/htslib/

../lib/samtools/bedidx.o:  ../lib/samtools/sam.h
	cd ../lib/samtools/ && make && cd ../..

../lib/samtools/sam.h:
	rm -rf ../lib/samtools/
	git clone --recursive https://github.com/samtools/samtools.git ../lib/samtools/


../lib/libgab/utils.h:
	rm -rf libgab/
	git clone --recursive https://github.com/grenaud/libgab.git ../lib/libgab/

../lib/libgab/utils.o: ../lib/bamtools/build/src/api/libbamtools.a  ../lib/libgab/utils.h
	make -C ../lib/libgab

../lib/bamtools/src/api/BamAlignment.h:
	rm -rf ../lib/bamtools/
	mkdir ../lib/bamtools/
	git clone --recursive https://github.com/pezmaster31/bamtools.git ../lib/bamtools/ && cd ../lib/bamtools/ && git reset --hard d24d850de17134fe4e7984b26493c5c0a1844b35

../lib/bamtools/build/src/api/libbamtools.a: ../lib/bamtools/src/api/BamAlignment.h
	cd ../lib/bamtools/ &&  mkdir -p build/  && cd build/ && cmake .. && make && cd ../..

clean:
	rm -f leeHom leeHomMulti *.o

