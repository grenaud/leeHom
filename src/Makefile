
BAMTOOLS= $(realpath ../bamtools/)
LIBGAB= $(realpath ../libgab/)



CXXFLAGS+= -O3 -Wall -I${BAMTOOLS}/include/ -I${LIBGAB} -I${LIBGAB}/gzstream/  
LDLIBS+= ${BAMTOOLS}/lib/libbamtools.a ${LIBGAB}/utils.o -lm -lz


all:    leeHom

%.o: %.cpp
	${CXX} -c ${CXXFLAGS} $^ 


leeHom: leeHom.o  MergeTrimReads.o ${LIBGAB}/PutProgramInHeader.o ${LIBGAB}/gzstream/libgzstream.a ${LIBGAB}/FastQObj.o  ${LIBGAB}/FastQParser.o
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)




clean:
	rm -f leeHom *.o

