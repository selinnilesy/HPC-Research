all: csr_encoder SSSconflictFree

csr_encoder: CSR_encoder.o
	g++ -o csr_encoder CSR_encoder.o -std=c++17  -lstdc++fs

CSR_encoder.o: CSR_encoder.cpp header.h
	g++ -Wall -g -c CSR_encoder.cpp header.h

SSSconflictFree: SSSconflictFree.cpp
	g++ SSSconflictFree.cpp -o SSSconflictFree -I/home/selin/xianyi-OpenBLAS-0b678b1/ -L/home/selin/xianyi-OpenBLAS-0b678b1/ -Wl,-rpath,/home/selin/xianyi-OpenBLAS-0b678b1/ -lopenblas  -std=c++17

clean:
	 rm CSR_encoder.o
	 rm csr_encoder
	 rm SSSconflictFree