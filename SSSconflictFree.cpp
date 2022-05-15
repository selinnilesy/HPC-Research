#include <cblas.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <cmath>

using namespace std;
vector<int> colVec;
vector<int> rowVec;
vector<double> valVec;

int readCooFormat(int z, double ratio) {
    cout <<  " start reading coo files..." << endl;
    double doubleVal;
    int intVal;
    const fs::path matrixFolder{"/home/selin/Split-Data/" + matrix_names[z]};
    for(auto const& dir_entry: fs::directory_iterator{matrixFolder}) {
        if (dir_entry.path().stem == "inner") {
            string rowFileName = c + "/coordinate-" + ratio + "-row.txt";
            string colFileName = dir_entry.path() + "/coordinate-" + ratio + "-col.txt";
            string valFileName = dir_entry.path() + "/coordinate-" + ratio + "-val.txt";

            std::fstream myfile(rowFileName, std::ios_base::in);
            // else, start reading doubles.
            while (myfile >> intVal) {
                rowVec.push_back(intVal);
            }
            myfile.close();
            cout << rowFileName << " has been read." << endl;

            myfile.open(colFileName, std::ios_base::in);
            while (myfile >> intVal) {
                colVec.push_back(intVal);
            }
            myfile.close();
            cout << colFileName << " has been read." << endl;

            myfile.open(valFileName, std::ios_base::in);
            while (myfile >> doubleVal) {
                valVec.push_back(doubleVal);
            }
            myfile.close();
            cout << valFileName << " has been read." << endl;
        }
    }
    else {
        cout << "not you wanted: " << dir_entry.path().stem ;
    }
    return 0;
}

int main(int argc, char **argv)
{

    int n = matrixSize[atoi(argv[1])];
    int inputType = atoi(argv[1]);
    double inputRatio = atof(argv[2]);
    cout << "input ratio: " << inputRatio << endl;
    readCooFormat(inputType, inputRatio);

    int i=0;
    int size= 10;

    float* X = new float[size];
    for(int i=0; i<size; i++) X[i] = 1.0;
    float* Y = new float[size];
    for(int i=0; i<size; i++) Y[i] = 0.0;
    float alpha = 1;
    float beta = 0;
    int k = 2;
    int lda = k+1;
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
    // for upper storage of A with size=10

    /*
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
    // col major - upper  10 x 3

    A[0]=0;
    A[1]=0;
    A[lda]=0;

    // col major - lower  10 x 3
    /*
    A[size*lda-1]=0;
    A[size*lda-2]=0;
    A[(size-1)*lda-1]=0;
     */



    cout << "Formed A: " << endl;

    for(int i=0; i<size; i++) {
        for(int j=0 ; j<lda; j++){
            cout << A[(lda)*i + j] << " " ;
        }
        cout <<  endl;
    }
    cout << "Call cblas_ssbmv. " << endl ;


    // column major : upper verince kernel lower'a giriyor. lower verince upper'a.
    // column major : CblasLower -> ( ! kernel'de upper) 10x10 1 bandwith icin calisiyor.
    cblas_ssbmv(CblasColMajor, CblasUpper, size, k, alpha, A, lda, X, incx, beta, Y, incy);

    cout << "Output: " << endl;
    for(int j=0; j<size; j++) {
        cout << Y[j] << endl;
    }
}
