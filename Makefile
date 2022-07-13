sall: csr_encoder SSSconflictFree

csr_encoder: CSR_encoder.o
	g++ -o csr_encoder CSR_encoder.o -std=c++17  -lstdc++fs

CSR_encoder.o: CSR_encoder.cpp header.h
	g++ -Wall -g -c CSR_encoder.cpp header.h

SSSconflictFree: SSSconflictFree.cpp
	cd /home/selin/SPARSKIT2 && make clean && make all
	mpic++ ../SPARSKIT2/FORMATS/formats.o  SSSconflictFree.cpp  header.h ../SPARSKIT2/MATGEN/FDIF/functns.o ../SPARSKIT2/libskit.a -o SSSconflictFree -lgfortran -lstdc++fs -std=c++17

clean:
	 rm CSR_encoder.o
	 rm csr_encoder
	 rm SSSconflictFree