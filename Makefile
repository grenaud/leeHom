all: 
	make -C libgab
	make -C src


clean:
	make -C libgab clean
	make -C src clean


.PHONY: all
