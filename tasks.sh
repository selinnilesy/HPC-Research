#!/bin/bash

for i in  0 1 2 3 4 5
do
	for nthread  in 2 3 4 5 6
	do
		 mpiexec -n $((2 **nthread))  ./SSSconflictFree $i $((2 **nthread))
	done
done
