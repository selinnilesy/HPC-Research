#include <cblas.h>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <string>
#include <cmath>
#include "header.h"

using namespace std;

vector<double> inner_banded;
vector<double> middle_banded;

int readBandedStorage(int z, double ratio, double middleRatio, bool inner, double *ptr) {
    cout <<  " start reading banded file..." << endl;
    double doubleVal;
    int intVal;

    string fileName;
    int counter=0;
    if(inner) fileName = "/home/selin/Split-Data/" + matrix_names[z]  + "/inner-outer-equal";
    else fileName = "/home/selin/Split-Data/" + matrix_names[z];

        if (inner) {
            fileName += "/inner-banded-A" + to_string(ratio) + ".txt";

            std::fstream myfile(fileName, std::ios_base::in);
            // else, start reading doubles.
            while (myfile >> doubleVal) {
                inner_banded.push_back(doubleVal);
                ptr[counter++] = doubleVal;
            }
            myfile.close();
            cout << fileName << " has been read with size: " << inner_banded.size() << endl;
        }
        // two parameters used. we read middle for this.
        else if (!inner) {
            fileName +=  + "/middle-banded-A" + to_string(ratio)+ "-" + to_string(middleRatio) + ".txt";

            std::fstream myfile(fileName, std::ios_base::in);
            // else, start reading doubles.
            while (myfile >> doubleVal) {
                middle_banded.push_back(doubleVal);
                ptr[counter++] = doubleVal;
            }
            myfile.close();
            cout << fileName << " has been read with size: " << middle_banded.size() << endl;
        }
    return 0;
}


int main(int argc, char **argv)
{
    int n = matrixSize[atoi(argv[1])];
    int inputType = atoi(argv[1]);
    double inputRatio = atof(argv[2]);
    double middleRatio = atof(argv[3]);
    bool inner = atof(argv[4]);
    cout << "input ratio: " << inputRatio << endl;
    cout << "middle ratio: " << middleRatio << endl;
    cout << "inner bool: " << inner << endl;
    // inner read = 1 , middle read = 0 !!!;

    int innerBandwith,middleBandwith;
    innerBandwith = (int) (nnz_n_Ratios[inputType]*bandwithProportions[inputType] * inputRatio);
    //middleBandwith = (int) ((bandwithSize[inputType] - innerBandwith)*middleRatio);

    cout << "inner bandwith: " << innerBandwith << endl;
    cout << "middle bandwith: " << middleBandwith << endl;
    //cout << "total bandwith: " << bandwithSize[inputType] << endl;

    int i,j,size,k,lda,size_1,size_2;
    double** A;

    if(inner) {
        k = innerBandwith;
        lda = k+1;
        size= n;
        size_1 = size;
        size_2=lda;
        A = new double*[size];
        for( i=0; i<size; i++) {
            A[i]  = new double[lda];
        }
        for( i=0; i<size; i++) {
            for( j=0; j<lda; j++) {
                A[i][j]  = 0.0;
            }
        }
    }
    else if(!inner){
        k = middleBandwith;
        lda = k+1;
        size_1 = size-innerBandwith-1;
        size_2=middleBandwith;
        A = new double*[size_1];
        for( i=0; i<size_1; i++) {
            A[i]  = new double[size_2];
        }

        for( i=0; i<size_1; i++) {
            for( j=0; j<size_2; j++) {
                A[i][j]  = 0.0;
            }
        }
    }
    cout << "A initialized." <<  endl;

    double* X = new double[size_1];
    for(int i=0; i<size_1; i++) X[i] = 1.0;
    double* Y = new double[size_1];
    for(int i=0; i<size_1; i++) Y[i] = 0.0;
    double alpha = 1;
    double beta = 0;
    int incx = 1;
    int incy = 1;

    double *B = *A;
    readBandedStorage(inputType, inputRatio, middleRatio, inner, B);
    cout << "banded file read onto B." <<  endl;

    cout << "Call cblas_dsbmv. " << endl ;
    // column major : upper verince kernel lower'a giriyor. lower verince upper'a.
    // column major : CblasLower -> ( ! kernel'de upper) 10x10 1 bandwith icin calisiyor.
    cblas_dsbmv(CblasColMajor, CblasUpper, size_1, k, alpha, B, lda, X, incx, beta, Y, incy);

    ofstream myfile;
    string output;
    if(inner) output = "/home/selin/Outputs/" + matrix_names[inputType] + "/inner-"  + to_string(inputRatio) + ".txt";
    if(!inner) output =  "/home/selin/Outputs/" + matrix_names[inputType] + "/middle-"  + to_string(inputRatio) + "-" + to_string(middleRatio) + ".txt";
    myfile.open(output, ios::out | ios::trunc);

    cout << "Writing Y: " << endl;
    for( i=0; i<size_1; i++) {
        myfile << Y[i] << " " ;
    }
    myfile.close();
    cout << "Output completed." << endl;

    delete [] X;
    delete [] Y;
    delete [] A;

}