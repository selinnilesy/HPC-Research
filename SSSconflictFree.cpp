#include <cblas.h>
#include <stdio.h>
#include <iostream>
#include <cmath>

using namespace std;


int main()
{
    int i=0;
    int size= 10;
    double** matrix = new double*[size];
    double** matrix2 = new double*[5];

    // generate symmetric banded matrix

    for(int i=0; i<10; i++){
        matrix[i] = new double[10];
        for(int j=0; j<10; j++) {
            if(j < i-2 || j> i+2 ) matrix[i][j] = 0.0;
            else if(i==j){
                matrix[i][j] = 1.0;
                if(i > 0 && j < 9){
                    matrix[i][j-1] = 5.0;
                    matrix[i][j+1] = 5.0;
                }

                if(i > 1 && j < 8){
                    matrix[i][j-2] = 2.0;
                    matrix[i][j+2] = 2.0;
                }

            }
        }
    }
    matrix[9][8] = 5.0;
    matrix[9][7] = 2.0;
    matrix[0][1] = 5.0;
    matrix[0][2] = 2.0;
    matrix[8][6] = 2.0;
    matrix[1][3] = 2.0;


    /*
    for(int i=0; i<5; i++){
        matrix[i] = new float[5];
        for(int j=0; j<5; j++) {
            if(j < i-1 || j> i+1 ) matrix[i][j] = 0.0;
            else if(i==j){
                matrix[i][j] = 1.0;
                if(i > 0 && j < 4){
                    matrix[i][j-1] = 1.0;
                    matrix[i][j+1] = 1.0;
                }
            }
        }
    }
    matrix[0][1] = matrix[4][3]  = 1.0;
    */


    for(int i=0; i<5; i++){
        matrix2[i] = new double[5];
        for(int j=0; j<5; j++) {
            matrix2[i][j] = 1.0;
        }
    }
    matrix2[0][3] = matrix2[0][4] = matrix2[1][4]  = matrix2[3][0] = matrix2[4][0] = matrix2[4][1] =0.0;

    cout << "Matrix" << endl;
    for(int i=0; i<size; i++){
        for(int j=0; j<size; j++) {
            cout << matrix[i][j] << " " ;
        }
        cout << endl;
    }
    cout << endl;
    /*
    cout << "Matrix 2" << endl;
    for(int i=0; i<5; i++){
        for(int j=0; j<5; j++) {
            cout << matrix2[i][j] << " " ;
        }
        cout << endl;
    }
     */

    double* X = new double[size];
    for(int i=0; i<size; i++) X[i] = 1.0;
    double* Y = new double[size];
    for(int i=0; i<size; i++) Y[i] = 0.0;
    double alpha = 1;
    double beta = 0;
    //int k = 2;
    int kl = 2;
    int ku=2;
    //int lda = k+1;
    int lda=kl+ku+1;
    int incx = 1;
    int incy = 1;
    double* A = new double[lda*size];
    double* B = new double[lda*size-1];

    int m;
    // for LOWER only


    /*
    for(int j=1; j<=size; j++) {
        m = 1 - j;
        for(int i=j ; i<=min(size, j+k); i++){
            int ind_x = i-1;
            int ind_y = j-1;
            A[(m+i-1)*size + ind_y] = matrix[ind_x][ind_y];
        }
    }


    // for upper storage of A with size=10


    A[0]=0;
    A[1]=0;
    A[size]=0;
    A[size -1 + size]=1;
    A[2*size + size - 1]=1;
    A[2*size + size - 2]=1;
     */
    // for upper storage of A with size=5 k=1
    /*
    A[0]=0;
    A[2*size -1 ]=1;
     */
    // for lower storage of A with size=10 k=2
    /*
    A[0]=1;
    A[1]=1;
    A[size]=1;
    A[3*size - 1 ]=0;
    A[3*size - 2 ]=0;
    A[2*size - 1 ]=0;
    */


    for(int i=0; i<size; i++) {
        for(int j=0 ; j<lda; j++){
            A[(lda)*i + j]= 1;
        }
    }


   for(int i=0; i<size; i++) {
      A[i*lda + 2]= matrix[i][i];
      if(i >=1) A[(i)*lda + 1]= matrix[i][i-1];
      if(i >=2) A[(i)*lda +0]= matrix[i][i-2];
   }
   A[0] =0;
   A[1]=0;
   A[lda]=0;
    for(int i=0; i<size; i++) {
        A[(i)*lda + (lda-2)]= matrix[i][i+1];
        A[(i)*lda +(lda-1)]= matrix[i][i+2];
    }

    A[lda*(size-1) + lda-1]=0;
    A[lda*(size-1) + lda-1]=0;
    A[(lda)*(size-2) + lda-1]=0;


    /*
    for(int i=0; i<size; i++) {
        A[i*lda + (lda-1)]= matrix[i][i];
        if(i >=1) A[(i)*lda + (lda-2)]= matrix[i][i-1];
        if(i >=2) A[(i)*lda +(lda-3)]= matrix[i][i-2];
    }
    for(int i=0; i<size-1; i++) {
        if(i >=1) B[(i)*lda + 0]= matrix[i][i-2];
    }
    B[(1)*lda + 0] = 2;
     */


    // col major - upper  10 x 3

    //A[0]=0;
    //A[1]=0;
    //A[lda]=0;

    // col major - lower  10 x 3

    //A[size*lda-1]=0;
    //A[size*lda-2]=0;
    //A[(size-1)*lda-1]=0;


    cout << "Formed A: " << endl;

    for(int i=0; i<size; i++) {
        for(int j=0 ; j<lda; j++){
            cout << A[(lda)*i + j] << " " ;
        }
        cout <<  endl;
    }
    /*
     * cout << "Formed B: " << endl;

    for(int i=0; i<size-1; i++) {
        for(int j=0 ; j<lda; j++){
            cout << B[(lda)*i + j] << " " ;
        }
        cout <<  endl;
    }
     */
    cout << "Call cblas_dgbmv. " << endl ;
    // CblasRowMajor = 101
    // CblasLower = 122
    // CBLAS_TRANSPOSE
    //const CBLAS_UPLO Uplo = CblasLower;
    // column major : upper verince kernel lower'a giriyor. lower verince upper'a.
    // column major : CblasLower -> ( ! kernel'de upper) 10x10 1 bandwith icin calisiyor.
    cblas_dgbmv(CblasRowMajor, CblasNoTrans , size, size, kl, ku, alpha, A, lda, X, incx, beta, Y, incy);

    cout << "Output: " << endl;
    for(int j=0; j<size; j++) {
        cout << Y[j] << endl;
    }
}
