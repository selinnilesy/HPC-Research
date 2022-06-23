#!/bin/bash

for thread in  5 
do
	for i  in 0 1 2 3 4 5
	do
		./SSSconflictFree $i 0 0 #$((2 ** thread))

	done
	for i  in 0 1 2 3 4 5
	do
		./SSSconflictFree $i 1 0 # $((2 ** thread))
	done
done

