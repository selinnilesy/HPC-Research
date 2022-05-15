all: csr_encoder SSSconflictFree

csr_encoder: CSR_encoder.o
	g++ -o csr_encoder CSR_encoder.o -std=c++17  -lstdc++fs

CSR_encoder.o: CSR_encoder.cpp header.h
	g++ -Wall -g -c CSR_encoder.cpp header.h

SSSconflictFree: SSSconflictFree.cpp header.h
	cd /home/selin/SPARSKIT2 && make clean && make all
	g++  ../SPARSKIT2/FORMATS/formats.o  SSSconflictFree.cpp  header.h ../SPARSKIT2/MATGEN/FDIF/functns.o ../SPARSKIT2/libskit.a  -o SSSconflictFree -std=c++17 -lstdc++fs -lgfortran

clean:
	 rm CSR_encoder.o
	 rm csr_encoder
	 rm SSSconflictFree