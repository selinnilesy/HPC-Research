all: csr_encoder SSSconflictFree

csr_encoder: CSR_encoder.o
	g++ -o csr_encoder CSR_encoder.o -std=c++17  -lstdc++fs

CSR_encoder.o: CSR_encoder.cpp header.h
	g++ -Wall -g -c CSR_encoder.cpp header.h

SSSconflictFree: SSSconflictFree.cpp header.h rcm.f90
	gfortran -c rcm.f90
	mpic++  --showme:compile SSSconflictFree.cpp  header.h -c SSSconflictFree.o -std=c++17 -lstdc++fs
	mpic++ rcm.o SSSconflictFree.o -o SSSconflictFree -lgfortran -lstdc++fs -std=c++17

clean:
	 rm CSR_encoder.o
	 rm csr_encoder
	 rm SSSconflictFree
