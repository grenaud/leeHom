#main
LIBGAB= $(realpath  libgab/)

PASSVARIABLES=

ifdef CFLAGS
	PASSVARIABLES+= CFLAGS="${CFLAGS}" 
endif

ifdef LDFLAGS
	PASSVARIABLES+= LDFLAGS="${LDFLAGS}" 
endif

ifdef HTSLIBINC
	PASSVARIABLES+= HTSLIBINC="${HTSLIBINC}" 
endif

ifdef HTSLIBLIB
	PASSVARIABLES+= HTSLIBLIB="${HTSLIBLIB}" 
endif

ifdef LIBGABINC
	PASSVARIABLES+= LIBGABINC="${LIBGABINC}" 
endif

ifdef LIBGABLIB
	PASSVARIABLES+= LIBGABLIB="${LIBGABLIB}" 
endif


ifeq ($(CXX),)
CXX := g++ #-g  -pg 
endif

print-%: ; @echo $*=$($*)

#ifeq ($(strip $(LIBGABINC)),)
#	PASSVARIABLES =  "BAMTOOLSINC="${BAMTOOLSINC}" BAMTOOLSLIB="${BAMTOOLSLIB}" LIBGABINC="${LIBGABINC}" LIBGABLIB="${LIBGABLIB}""
#endif

all:
	make  CXX="${CXX}" ${PASSVARIABLES}  -C src

clean:
	make -C src clean

test:   all
	make  CXX="${CXX}" ${PASSVARIABLES}  -C src fuzz_kernels
	cd test/ && bash test.sh && cd ..

.PHONY: all
