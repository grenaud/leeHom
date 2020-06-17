all: 
	make -C src

clean:
	make -C src clean

test:   all
	cd test/ && bash test.sh && cd ..

.PHONY: all
