#src
HTSLIB= $(abspath ../htslib)
LIBGAB= $(realpath ../libgab/)


#defined via user
ifdef HTSLIBINC
	HTSLIBINCDIR=${HTSLIBINC}
else
	HTSLIBINCDIR=${HTSLIB}
endif

ifdef HTSLIBLIB
	HTSLIBLIBOBJ=-L${HTSLIBLIB} -lhts
	HTSLIBDEP=
else
	HTSLIBLIBOBJ=${HTSLIB}/libhts.a
	HTSLIBDEP=${HTSLIB}/libhts.a
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
CXXFLAGS+= -O3 -ffp-contract=off -pthread -Wall -Wno-char-subscripts -I${HTSLIBINCDIR} -I${LIBGABINC} -I${LIBGABINC}/gzstream/  
LDFLAGS+= -lm -lpthread -lz

LDLIBS+=  ${LIBGABLIBOBJ} ${LIBGABLIB}/gzstream/libgzstream.a ${HTSLIBLIBOBJ}

OBJS= leeHom.o MergeTrimReads.o Kernels.o BamRecord.o

all:    ${HTSLIBDEP} ${LIBGABLIBOBJ} leeHom

%.o: %.cpp
	${CXX} -c ${CXXFLAGS} $< 

${OBJS}: MergeTrimReads.h Kernels.h BamRecord.h Pipeline.h
leeHom.o: version.h

version.h:
	./get_version.sh

leeHom: version.h ${OBJS} ${LIBGABLIBOBJ} ${HTSLIBDEP}
	@echo "linking"
	@echo ${LDLIBS}
	@echo ${LDFLAGS}
	${CXX} -o $@ ${OBJS} $(LDLIBS) $(LDFLAGS) 

static: version.h ${OBJS} ${LIBGABLIBOBJ} ${HTSLIBDEP}
	@echo "staic linking"
	@echo ${LDLIBS}
	@echo ${LDFLAGS}
	${CXX} -static -o leeHom ${OBJS} $(LDLIBS) $(LDFLAGS) 


# checks that the vectorized kernels give the same results as the original scalar code, see ../test/test.sh
fuzz_kernels: ../test/fuzz_kernels.cpp MergeTrimReads.o Kernels.o BamRecord.o ${LIBGABLIBOBJ} ${HTSLIBDEP}
	${CXX} ${CXXFLAGS} -I. -o $@ $< MergeTrimReads.o Kernels.o BamRecord.o $(LDLIBS) $(LDFLAGS)

${LIBGABINC}/libgab.h:
	rm -rf ../libgab/
	mkdir ../libgab/
	git clone --recursive https://github.com/grenaud/libgab.git ../libgab/

${LIBGABLIBOBJ}: ${LIBGABINC}/libgab.h
	make -C ../libgab

# htslib is a submodule, it is built as a static library with only what is needed for BAM and gzip (no libcurl, bz2, lzma)
${HTSLIB}/htslib/sam.h:
	cd .. && (git submodule update --init --recursive htslib || git clone --recursive --branch 1.22.1 https://github.com/samtools/htslib.git htslib)

${HTSLIB}/libhts.a: ${HTSLIB}/htslib/sam.h
	cd ../htslib && (test -f config.mk || (autoreconf -i && ./configure --disable-libcurl --disable-plugins --disable-bz2 --disable-lzma --disable-gcs --disable-s3 --without-libdeflate)) && $(MAKE) libhts.a

clean:
	rm -f leeHom fuzz_kernels *.o
