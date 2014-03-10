
BAMTOOLS= $(realpath ../bamtools/)
LIBGAB= $(realpath ../libgab/)



CXXFLAGS+= -O3 -Wall -I${BAMTOOLS}/include/ -I${LIBGAB} -I${LIBGAB}/gzstream/  
LDLIBS+= ${BAMTOOLS}/lib/libbamtools.a ${LIBGAB}/utils.o -lm -lz


all:    mergeTrimReadsBAM

%.o: %.cpp
	${CXX} -c ${CXXFLAGS} $^ 


mergeTrimReadsBAM: mergeTrimReadsBAM.o  MergeTrimReads.o ${LIBGAB}/PutProgramInHeader.o ${LIBGAB}/gzstream/libgzstream.a ${LIBGAB}/FastQObj.o  ${LIBGAB}/FastQParser.o
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)




clean:
	rm -f mergeTrimReadsBAM *.o

