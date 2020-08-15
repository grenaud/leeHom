
#BAMTOOLS= $(realpath ../bamtools/)
#LIBGAB= $(realpath ../libgab/)


CXX      = g++   #-g -pg
CXXFLAGS+= -O3 -Wall -Wno-char-subscripts -I../bamtools/include/ -I../bamtools/src/ -I../libgab -I../libgab/gzstream/  
LDFLAGS  = -lm -lz  -lpthread 
LDLIBS+=  ../bamtools/build/src/utils/CMakeFiles/BamTools-utils.dir/*cpp.o ../bamtools/build/src/api/libbamtools.a ../libgab/utils.o 



all:    ../bamtools/build/src/api/libbamtools.a ../libgab/utils.o MergeTrimReads.o leeHom

%.o: %.cpp
	${CXX} -c ${CXXFLAGS} $^ 

leeHom: leeHom.o  MergeTrimReads.o ../libgab/PutProgramInHeader.o  ../libgab/FastQObj.o  ../libgab/FastQParser.o ../libgab/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

../libgab/utils.h:
	rm -rf libgab/
	git clone --recursive https://github.com/grenaud/libgab.git ../libgab/

../libgab/utils.o: ../bamtools/build/src/api/libbamtools.a  ../libgab/utils.h
	make -C ../libgab

../bamtools/src/api/BamAlignment.h:
	rm -rf ../bamtools/
	mkdir ../bamtools/
	git clone --recursive https://github.com/pezmaster31/bamtools.git ../bamtools/ && cd ../bamtools/ && git reset --hard d24d850de17134fe4e7984b26493c5c0a1844b35

../bamtools/build/src/api/libbamtools.a: ../bamtools/src/api/BamAlignment.h
	cd ../bamtools/ &&  mkdir -p build/  && cd build/ && if cmake ..; then echo ""; else if cmake3 ..; then echo ""; else echo "cmake failed, please install cmake v3"; fi  fi  && make && cd ../..

clean:
	rm -f leeHom *.o

