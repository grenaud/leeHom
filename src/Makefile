
#BAMTOOLS= $(realpath ../bamtools/)
#LIBGAB= $(realpath ../libgab/)


CXX      = g++   #-g -pg
CXXFLAGS+= -O3 -Wall -I../bamtools/include/ -I../bamtools/src/ -I../libgab -I../libgab/gzstream/  
LDFLAGS  = -lm -lz  -lpthread 
LDLIBS+=  ../bamtools/build/src/utils/CMakeFiles/BamTools-utils.dir/*cpp.o ../bamtools/build/src/api/libbamtools.a ../libgab/utils.o 



all:    ../bamtools/build/src/api/libbamtools.a ../libgab/utils.o MergeTrimReads.o leeHom leeHomMulti

%.o: %.cpp
	${CXX} -c ${CXXFLAGS} $^ 

leeHom: leeHom.o  MergeTrimReads.o ../libgab/PutProgramInHeader.o ../libgab/FastQObj.o  ../libgab/FastQParser.o  ../libgab/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)


leeHomMulti: leeHomMulti.o  MergeTrimReads.o ../libgab/PutProgramInHeader.o  ../libgab/FastQObj.o  ../libgab/FastQParser.o ../libgab/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

../libgab/utils.h:
	rm -rf libgab/
	git clone --recursive https://github.com/grenaud/libgab.git ../libgab/

../libgab/utils.o: ../bamtools/build/src/api/libbamtools.a  ../libgab/utils.h
	make -C ../libgab

../bamtools/src/api/BamAlignment.h:
	rm -rf ../bamtools/
	mkdir ../bamtools/
	git clone --recursive https://github.com/pezmaster31/bamtools.git ../bamtools/

../bamtools/build/src/api/libbamtools.a: ../bamtools/src/api/BamAlignment.h
	cd ../bamtools/ && mkdir -p build/  && cd build/ && cmake .. && make && cd ../..

clean:
	rm -f leeHom leeHomMulti *.o

