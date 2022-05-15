#include <cblas.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <cmath>

using namespace std;
vector<int> colVec;
vector<int> rowVec;
vector<double> valVec;
vector<double> diagVec;

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
        else {
            cout << "not you wanted: " << dir_entry.path().stem ;
        }
    }
    return 0;
}
int readDiag(int z) {
    cout <<  " start reading diag file..." << endl;
    double doubleVal;
    string diagFile = "/home/selin/SSS-Data/" + matrix_names[z] + "/diag.txt" ;

    std::fstream myfile(diagFile, std::ios_base::in);
    // else, start reading doubles.
    while (myfile >> doubleVal) {
        diagVec.push_back(doubleVal);
    }
    myfile.close();
    cout << diagFile << " has been read." << endl;
    return 0;
}

int main(int argc, char **argv)
{

    int n = matrixSize[atoi(argv[1])];
    int inputType = atoi(argv[1]);
    double inputRatio = atof(argv[2]);
    cout << "input ratio: " << inputRatio << endl;
    readCooFormat(inputType, inputRatio);
    readDiag(inputType);

    if(rowVec.size() != colVec.size()){
        cout << "not equal row and col (NNZ) !!! " << endl;
        return -1;
    }

    double innerBandwith, middleBandwith;
    innerBandwith = nnz_n_Ratios[inputType]*bandwithProportions[inputType] * inputRatio;
    middleBandwith = bandwithSize[inputType] - 2*innerBandwith;
    cout << "inner bandwith: " << innerBandwith << endl;
    cout << "middle bandwith: " << middleBandwith << endl;
    cout << "total bandwith: " << bandwithSize[inputType] << endl;

    int i,j;
    int size= n;
    // this is except diagonal.
    int k = innerBandwith;
    int lda = k+1;

    float** A = new float*[size];
    for( i=0; i<size; i++) {
        A[i]  = new float[lda];
    }
    memset(A, 0, lda*size*sizeof(float));


    int row, col, val, counter=0;
    // i keeps track of whole element count.
    for( i=0; i<rowVec.size(); i++) {
        row = rowVec[i] - 1;
        col = colVec[i] - 1;
        val = valVec[i];
        if(counter==innerBandwith) counter=0;

        if(row < innerBandwith){
            // insert first element
            A[row-2][innerBandwith-1 - i] =  val;
            for(int k=0; k<row-2; k++){
                val = valVec[i + (k+1)];
                A[row-2][innerBandwith-1 - i + (k+1)] =  val;
            }
            i+=k;
        }
        // soldan saga doldurmaya baslayabilirsin artik.
        else{
            A[row-2][counter] =  val;
            counter++;
        }
        cout <<  endl;
    }
    ofstream myfile;
    string output =  "/home/selin/HPC-Research/banded-A.txt";
    myfile.open(output, ios::out | ios::trunc);


    cout << "Formed A: " << endl;
    for( i=0; i<size; i++) {
        for( j=0 ; j<lda; j++){
            myfile << A[i][j] << " " ;
        }
        myfile <<  endl;
        myfile <<  endl;
    }

    /* needed data for openblas ssbmv
   float* X = new float[size];
   for(int i=0; i<size; i++) X[i] = 1.0;
   float* Y = new float[size];
   for(int i=0; i<size; i++) Y[i] = 0.0;
   float alpha = 1;
   float beta = 0;
   int incx = 1;
   int incy = 1;
    */

    //cout << "Call cblas_ssbmv. " << endl ;

    // column major : upper verince kernel lower'a giriyor. lower verince upper'a.
    // column major : CblasLower -> ( ! kernel'de upper) 10x10 1 bandwith icin calisiyor.
    //cblas_ssbmv(CblasColMajor, CblasUpper, size, k, alpha, A, lda, X, incx, beta, Y, incy);

    //cout << "Output: " << endl;
    //for(int j=0; j<size; j++) {
    //    cout << Y[j] << endl;
    //}
}
