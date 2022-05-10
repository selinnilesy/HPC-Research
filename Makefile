all: csr_encoder SSSconflictFree

csr_encoder: CSR_encoder.o
	g++ -o csr_encoder CSR_encoder.o -std=c++17  -lstdc++fs

CSR_encoder.o: CSR_encoder.cpp header.h
	g++ -Wall -g -c CSR_encoder.cpp header.h

SSSconflictFree: SSSconflictFree.cpp
	cd ./xianyi-OpenBLAS-0b678b1 && make &&  make PREFIX=/home/selin/HPC-Research/xianyi-OpenBLAS-0b678b1  install
	g++ SSSconflictFree.cpp -o SSSconflictFree -I/home/selin/HPC-Research/xianyi-OpenBLAS-0b678b1/include/ -L/home/selin/HPC-Research/xianyi-OpenBLAS-0b678b1/lib -Wl,-rpath,/home/selin/HPC-Research/xianyi-OpenBLAS-0b678b1/lib -lopenblas -std=c++17

clean:
	 rm CSR_encoder.o
	 rm csr_encoder
	 rm SSSconflictFree