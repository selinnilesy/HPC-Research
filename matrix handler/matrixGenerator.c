/*
*   Matrix Market I/O example program
*
*   Read a real (non-complex) sparse matrix from a Matrix Market (v. 2.0) file.
*   and copies it to stdout.  This porgram does nothing useful, but
*   illustrates common usage of the Matrix Matrix I/O routines.
*   (See http://math.nist.gov/MatrixMarket for details.)
*
*   Usage:  a.out [filename] > output
*
*
*   NOTES:
*
*   1) Matrix Market files are always 1-based, i.e. the index of the first
*      element of a matrix is (1,1), not (0,0) as in C.  ADJUST THESE
*      OFFSETS ACCORDINGLY offsets accordingly when reading and writing
*      to files.
*
*   2) ANSI C requires one to use the "l" format modifier when reading
*      double precision floating point numbers in scanf() and
*      its variants.  For example, use "%lf", "%lg", or "%le"
*      when reading doubles, otherwise errors will occur.
*/

#include <stdio.h>
#include <stdlib.h>
#include "mmio.h"

int main(int argc, char *argv[])
{
    int ret_code;
    MM_typecode matcode;
    FILE *f;
    int M, N, nz;
    int i, *I, *J;
    double *val;

    M=N=914898;
    nz=14589714;

     if ((f = fopen(argv[1], "r")) == NULL)
            exit(1);


    /* reseve memory for matrices */

    I = (int *) malloc(M * sizeof(int));
    J = (int *) malloc(N * sizeof(int));
    val = (double *) malloc(nz * sizeof(double));

    /************************/
    /* now write out matrix */
    /************************/

    for (i=0; i<M; i++){
        for (j=0; i<N; j++) {
                fscanf("%d %d %20.19g\n", I[i], J[i], val[i]);
        }
    }

    if (f !=stdin) fclose(f);

    return 0;
}
