
BAMTOOLS= $(realpath ../bamtools/)
LIBGAB= $(realpath ../libgab/)


<<<<<<< HEAD
CXX      = g++ #-g#-pg# -g
=======
CXX      = g++  -pg
>>>>>>> 544f462ee01ff394a11c36d394d35dc290f83120
CXXFLAGS+= -O3 -Wall -I${BAMTOOLS}/include/ -I${LIBGAB} -I${LIBGAB}/gzstream/  
LDFLAGS  = -lm -lz  -lpthread #-lprofiler
LDLIBS+=  ${BAMTOOLS}/build/src/utils/CMakeFiles/BamTools-utils.dir/*cpp.o ${BAMTOOLS}/lib/libbamtools.a ${LIBGAB}/utils.o 


all:    leeHom leeHomMulti

%.o: %.cpp
	${CXX} -c ${CXXFLAGS} $^ 


leeHom: leeHom.o  MergeTrimReads.o ${LIBGAB}/PutProgramInHeader.o ${LIBGAB}/FastQObj.o  ${LIBGAB}/FastQParser.o  ${LIBGAB}/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)


leeHomMulti: leeHomMulti.o  MergeTrimReads.o ${LIBGAB}/PutProgramInHeader.o  ${LIBGAB}/FastQObj.o  ${LIBGAB}/FastQParser.o ${LIBGAB}/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

clean:
	rm -f leeHom leeHomMulti *.o

