
BAMTOOLS= $(realpath ../bamtools/)
LIBGAB= $(realpath ../libgab/)


CXX      = g++ # -g -pg
CXXFLAGS+= -O3 -Wall -I${BAMTOOLS}/include/ -I${LIBGAB} -I${LIBGAB}/gzstream/  
LDFLAGS  = -lm -lz  -lpthread 
LDLIBS+=  ${BAMTOOLS}/build/src/utils/CMakeFiles/BamTools-utils.dir/*cpp.o ${BAMTOOLS}/lib/libbamtools.a ${LIBGAB}/utils.o 



all:    leeHom leeHomMulti

%.o: %.cpp
	${CXX} -c ${CXXFLAGS} $^ 

leeHom: leeHom.o  MergeTrimReads.o ${LIBGAB}/PutProgramInHeader.o ${LIBGAB}/FastQObj.o  ${LIBGAB}/FastQParser.o  ${LIBGAB}/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)


leeHomMulti: leeHomMulti.o  MergeTrimReads.o ${LIBGAB}/PutProgramInHeader.o  ${LIBGAB}/FastQObj.o  ${LIBGAB}/FastQParser.o ${LIBGAB}/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

${LIBGAB}/utils.o: ../lib/bamtools/lib/libbamtools.a  ${LIBGAB}/utils.h
	make -C ${LIBGAB}

${BAMTOOLS}/src/api/BamAlignment.h:
	rm -rf ${BAMTOOLS}/
	mkdir ${BAMTOOLS}/
	git clone --recursive https://github.com/pezmaster31/bamtools.git ${BAMTOOLS}/

${BAMTOOLS}/lib/libbamtools.a: ${BAMTOOLS}/src/api/BamAlignment.h
	cd ${BAMTOOLS}/ && mkdir -p build/  && cd build/ && cmake .. && make && cd ../..

clean:
	rm -f leeHom leeHomMulti *.o

