#main
BAMTOOLS= $(realpath bamtools/)
LIBGAB= $(realpath  libgab/)

PASSVARIABLES=

ifdef CFLAGS
	PASSVARIABLES+= CFLAGS="${CFLAGS}" 
endif

ifdef LDFLAGS
	PASSVARIABLES+= LDFLAGS="${LDFLAGS}" 
endif

ifdef BAMTOOLSINC
	PASSVARIABLES+= BAMTOOLSINC="${BAMTOOLSINC}" 
endif

ifdef BAMTOOLSLIB
	PASSVARIABLES+= BAMTOOLSLIB="${BAMTOOLSLIB}" 
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
	cd test/ && bash test.sh && cd ..

.PHONY: all
