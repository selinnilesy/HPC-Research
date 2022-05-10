#! /bin/sh
mpiexec -np 4  -mca btl ^openib SSSconflictFree 0
mpiexec -np 4  -mca btl ^openib SSSconflictFree 1
mpiexec -np 4  -mca btl ^openib SSSconflictFree 2
mpiexec -np 4  -mca btl ^openib SSSconflictFree 3
mpiexec -np 4  -mca btl ^openib SSSconflictFree 4
mpiexec -np 4  -mca btl ^openib SSSconflictFree 5
