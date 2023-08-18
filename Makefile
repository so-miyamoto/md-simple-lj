

all: clean
	make clean -C src
	make -C src
	mkdir -p bin
	ln -s ../src/build/nve.out bin
	ln -s ../src/build/nvt.out bin

clean:
	rm -rf ./bin/
