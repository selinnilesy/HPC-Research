#include <cblas.h>
#include <stdio.h>
#include <iostream>
#include <cmath>

using namespace std;


int main()
{
    int i=0;
    float** matrix = new float*[10];
    float** matrix2 = new float*[5];
    /*
    double matrix0[5] = {0.0,0.0,0.0,0.0,0.0,};
    double matrix1[5] ={0.0,0.0,0.0,0.0,0.0,};
    double matrix2[5] = {0.0,0.0,0.0,0.0,0.0,};
    double matrix3[5] = {0.0,0.0,0.0,0.0,0.0,};
    double matrix4[5] = {0.0,0.0,0.0,0.0,0.0,};
     */
    // generate symmetric banded matrix.
    for(int i=0; i<10; i++){
        matrix[i] = new float[10];
        for(int j=0; j<10; j++) {
            if(j < i-2 || j> i+2 ) matrix[i][j] = 0.0;
            else if(i==j){
                matrix[i][j] = 1.0;
                if(i > 0 && j < 9){
                    matrix[i][j-1] = 1.0;
                    matrix[i][j+1] = 1.0;
                }
                if(i > 1 && j < 8){
                    matrix[i][j-2] = 1.0;
                    matrix[i][j+2] = 1.0;
                }
            }
        }
    }
    matrix[0][1] = 1.0;
    matrix[1][3] = 1.0;
    matrix[0][2] = 1.0;
    matrix[9][8] = 1.0;
    matrix[8][6] = 1.0;
    matrix[9][7] = 1.0;
    //matrix[2][3] = 1.0;
    //matrix[3][4] = 1.0;
    //matrix[4][3] = 1.0;

    for(int i=0; i<5; i++){
        matrix2[i] = new float[5];
        for(int j=0; j<5; j++) {
            matrix2[i][j] = 1.0;
        }
    }
    matrix2[0][3] = matrix2[0][4] = matrix2[1][4]  = matrix2[3][0] = matrix2[4][0] = matrix2[4][1] =0.0;

    for(int i=0; i<10; i++){
        for(int j=0; j<10; j++) {
            cout << matrix[i][j] << " " ;
        }
        cout << endl;
    }
    cout << endl;
    for(int i=0; i<5; i++){
        for(int j=0; j<5; j++) {
            cout << matrix2[i][j] << " " ;
        }
        cout << endl;
    }

    float* X = new float[10];
    for(int i=0; i<10; i++) X[i] = 1.0;
    float* Y = new float[10];
    for(int i=0; i<10; i++) Y[i] = 0.0;
    float alpha = 1;
    float beta = 0;
    int size= 10;
    int k = 2;
    int lda = 3;
    int incx = 1;
    int incy = 1;
    float* A = new float[lda*size];

    int m;
    // for LOWER only
    for(int j=1; j<=size; j++) {
        m = 1 - j;
        for(int i=j ; i<=min(size, j+k); i++){
            int ind_x = i-1;
            int ind_y = j-1;
            A[(m+i-1)*size + ind_y] = matrix[ind_x][ind_y];
        }
    }
    /* for upper storage of A */
    /*
    A[0]=0;
    A[1]=0;
    A[size]=0;
    A[size -1 + size]=1;
    A[2*size + size - 1]=1;
    A[2*size + size - 2]=1;
     */

    cout << "Formed A: " << endl;

    for(int i=0; i<lda; i++) {
        for(int j=0 ; j<size; j++){
            cout << A[(i)*size + j] << " " ;
        }
        cout <<  endl;
    }
    cout << "Call cblas_ssbmv. " << endl ;
    // CblasRowMajor = 101
    // CblasLower = 122
    //const CBLAS_UPLO Uplo = CblasLower;
    // upper verince kernel lower'a giriyo. lower verince upper'a.
    cblas_ssbmv(CblasColMajor, CblasLower, size, k, alpha, A, lda, X, incx, beta, Y, incy);


    // MY C++ IMPLEMENTATION FOR SSBMV.F
    /*
    int temp1, temp2, l;
    for(int j=1; j<=size; j++){
        temp1 = X[j-1];
        temp2 = 0;
        //cout << "Y[j-1] with j-1: " << j-1 << endl;
        //cout << "1st body - add " << temp1*A[j-1] << endl;
        Y[j-1] = Y[j-1] + temp1*A[j-1];
        l = 1 - j;
        for(int i=j+1 ; i<=min(size, j+k); i++) {
            Y[i-1] = Y[i-1] + temp1 * A[(l+i -1) * size + j-1];
            //cout << "2nd body - add Y[i-1] with i: " << i << " " << temp1 * A[(l+i -1) * size + j-1] << endl;
            temp2 = temp2 - A[(l+i -1) * size + j-1] * X[i-1];
        }
        //cout << "3rd body - add " <<temp2 << endl;
        Y[j-1] = Y[j-1] + temp2;
    }
     */


    cout << "Output: " << endl;
    for(int j=0; j<10; j++) {
        cout << Y[j] << endl;
    }

}