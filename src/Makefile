
BAMTOOLS= $(realpath ../bamtools/)
LIBGAB= $(realpath ../libgab/)
SIMPLEJSON= $(realpath ../SimpleJSON/)



CXXFLAGS+= -O3 -Wall -I${BAMTOOLS}/include/ -I${LIBGAB} -I${LIBGAB}/gzstream/  -I${SIMPLEJSON}/src/
LDLIBS+= ${BAMTOOLS}/lib/libbamtools.a ${LIBGAB}/utils.o -lm -lz


all:    mergeTrimReadsBAM assignRG filterReads errorRatePerCycle 
	#make -C SimpleJSON
	#make -C libgab

%.o: %.cpp
	${CXX} -c ${CXXFLAGS} $^ 

assignRG: assignRG.o PrefixTree.o RGAssign.o ${LIBGAB}/PutProgramInHeader.o ${SIMPLEJSON}/obj/JSON.o ${SIMPLEJSON}/obj/JSONValue.o
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)

mergeTrimReadsBAM: mergeTrimReadsBAM.o  MergeTrimReads.o ${LIBGAB}/PutProgramInHeader.o ${LIBGAB}/gzstream/libgzstream.a ${LIBGAB}/FastQObj.o  ${LIBGAB}/FastQParser.o
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)

filterReads: filterReads.o FilterBAMal.o ${LIBGAB}/PutProgramInHeader.o 
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)

errorRatePerCycle: errorRatePerCycle.o   ${LIBGAB}/ReconsReferenceBAM.o
	${CXX} -o $@ $^ $(LDLIBS)




clean:
	rm -f mergeTrimReadsBAM  assignRG filterReads errorRatePerCycle *.o

