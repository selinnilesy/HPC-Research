all: csr_encoder SSSconflictFree

csr_encoder: CSR_encoder.o
	g++ -o csr_encoder CSR_encoder.o -std=c++17  -lstdc++fs

CSR_encoder.o: CSR_encoder.cpp header.h
	g++ -Wall -g -c CSR_encoder.cpp header.h

SSSconflictFree: SSSconflictFree.cpp
	mpic++ SSSconflictFree.cpp -o SSSconflictFree -lstdc++fs -std=c++17

clean:
	 rm CSR_encoder.o
	 rm csr_encoder
	 rm SSSconflictFree