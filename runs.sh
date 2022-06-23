#!/bin/bash

for nthread in  3 4 5 6
do
	for i  in 0 1 2 3 4 5
	do
		./SSSconflictFree $i 0 1 $((2 **nthread))

	done
	for i  in 0 1 2 3 4 5
	do
		./SSSconflictFree $i 1 1 $((2 ** nthread))
	done
done
