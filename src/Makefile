#src
BAMTOOLS= $(realpath ../bamtools/)
LIBGAB= $(realpath ../libgab/)


#defined via user
ifdef BAMTOOLSINC
	BAMTOOLSINCAPI=${BAMTOOLSINC}
	BAMTOOLSINCSHARED=${BAMTOOLSINC}
else
	BAMTOOLSINCAPI=${BAMTOOLS}/src/
	BAMTOOLSINCSHARED=${BAMTOOLS}/src/
endif

ifdef BAMTOOLSLIB
	BAMTOOLSLIBOBJ=${BAMTOOLSLIB}/libbamtools.a 
else
	BAMTOOLSLIBOBJ=${BAMTOOLS}/build/src/api/libbamtools.a 
endif

ifdef LIBGABINC
	LIBGABINC=${LIBGABINC}
else
	LIBGABINC=${LIBGAB}
endif

ifdef LIBGABLIB
	LIBGABLIBOBJ=${LIBGABLIB}/libgab.a 
else
	LIBGABLIBOBJ=${LIBGAB}/libgab.a
	LIBGABLIB=${LIBGAB}
endif


ifeq ($(CXX),)
CXX := g++ #-g  -pg 
endif

print-%: ; @echo $*=$($*)


#CXX      = ${CXX}#g++   #-g -pg
CXXFLAGS+= -O3 -Wall -Wno-char-subscripts -I${BAMTOOLSINCAPI} -I${BAMTOOLSINCSHARED}  -I${LIBGABINC} -I${LIBGABINC}/gzstream/  
LDFLAGS  = -lm -lpthread -lz

LDLIBS+=  ${BAMTOOLSLIBOBJ} ${LIBGABLIBOBJ} ${LIBGABLIB}/libgabbamtools.a  ${LIBGABLIB}/gzstream/libgzstream.a

all:    ${BAMTOOLSLIBOBJ} ${LIBGABLIBOBJ} MergeTrimReads.o leeHom

%.o: %.cpp
	${CXX} -c ${CXXFLAGS} $^ 

leeHom: leeHom.o  MergeTrimReads.o ${LIBGABLIBOBJ} 
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

${LIBGABINC}/libgab.h:
	rm -rf ../libgab/
	mkdir ../libgab/
	git clone --recursive https://github.com/grenaud/libgab.git ../libgab/

${LIBGABLIBOBJ}: ${BAMTOOLSLIBOBJ}  ${LIBGABINC}/libgab.h
	make -C ../libgab

${BAMTOOLSINCAPI}/api/BamAlignment.h:
	@echo ${BAMTOOLSINCAPI}/api/BamAlignment.h
	@echo ${BAMTOOLSLIBOBJ}
	exit
	rm -rf ../bamtools/
	mkdir ../bamtools/
	git clone --recursive https://github.com/pezmaster31/bamtools.git ../bamtools/ && cd ../bamtools/ && git reset --hard d24d850de17134fe4e7984b26493c5c0a1844b35

${BAMTOOLSLIBOBJ}: ${BAMTOOLSINCAPI}/api/BamAlignment.h
	cd ../bamtools/ &&  mkdir -p build/  && cd build/ && if cmake ..; then echo ""; else if cmake3 ..; then echo ""; else echo "cmake failed, please install cmake v3"; fi  fi  && make && cd ../..

clean:
	rm -f leeHom *.o

