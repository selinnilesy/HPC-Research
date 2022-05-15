#include <cblas.h>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <string>
#include <cmath>
#include "header.h"

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
        if (dir_entry.path().stem() == "inner") {
            string rowFileName = dir_entry.path().string() + "/coordinate-" + to_string(ratio) + "-row.txt";
            string colFileName = dir_entry.path().string()+ "/coordinate-" + to_string(ratio) + "-col.txt";
            string valFileName = dir_entry.path().string() + "/coordinate-" + to_string(ratio) + "-val.txt";

            std::fstream myfile(rowFileName, std::ios_base::in);
            // else, start reading doubles.
            while (myfile >> intVal) {
                rowVec.push_back(intVal);
            }
            myfile.close();
            cout << rowFileName << " has been read with size: " << rowVec.size() << endl;

            myfile.open(colFileName, std::ios_base::in);
            while (myfile >> intVal) {
                colVec.push_back(intVal);
            }
            myfile.close();
            cout << colFileName << " has been read with size: " << colVec.size() << endl;

            myfile.open(valFileName, std::ios_base::in);
            while (myfile >> doubleVal) {
                valVec.push_back(doubleVal);
            }
            myfile.close();
            cout << valFileName << " has been read with size: " << valVec.size() << endl;
        }
        else {
            cout << "not you wanted: " << dir_entry.path().stem() ;
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
    /*
     * [1 0 0 0 0    ]
     * [8 1 0 0 0    ]
     * [9 5 1 0 0    ]
     * [0 7 6 1 0    ]
     * [0 0 4 3 1    ]

    rowVec.push_back(2);
    rowVec.push_back(3);
    rowVec.push_back(3);
    rowVec.push_back(4);
    rowVec.push_back(4);
    rowVec.push_back(5);
    rowVec.push_back(5);

    valVec.push_back(8);
    valVec.push_back(9);
    valVec.push_back(5);
    valVec.push_back(7);
    valVec.push_back(6);
    valVec.push_back(4);
    valVec.push_back(3);
     */


    if(rowVec.size() != valVec.size()){
        cout << "not equal row and valVec (NNZ) !!! " << endl;
        return -1;
    }


    double innerBandwith, middleBandwith;
    innerBandwith = nnz_n_Ratios[inputType]*bandwithProportions[inputType] * inputRatio;
    middleBandwith = bandwithSize[inputType] - 2*innerBandwith;
    //innerBandwith=2;
    //middleBandwith=0;
    cout << "inner bandwith: " << innerBandwith << endl;
    cout << "middle bandwith: " << middleBandwith << endl;
    //cout << "total bandwith: " << bandwithSize[inputType] << endl;

    int i,j;
    int size= n;
    //int size=5;
    // this is except diagonal.
    int k = innerBandwith;
    int lda = k+1;

    float** A = new float*[size];
    for( i=0; i<size; i++) {
        A[i]  = new float[lda];
    }

    for( i=0; i<size; i++) {
        for( j=0; j<lda; j++) {
            A[i][j]  = 0.0;
        }
    }
    //memset(A, 0, k*k);

    cout << "A initialized to 0." <<  endl;
    int row, col, val, counter=0, x;
    // i keeps track of whole element count.
    // upper bound is not n, but instead n-1 !
    for( i=0; i<rowVec.size(); i++) {
        row = rowVec[i] - 1;
        //col = colVec[i] - 1;
        val = valVec[i];
        //cout << "row: " << row << " col: " << col <<  " val: " << val << endl;
        if(counter==innerBandwith) counter=0;

        if(row <= innerBandwith){
            // insert first element
            A[row-1][((int)innerBandwith)-1 - i] =  val;
            //cout << "inserted: " << val << endl;
            for(x=0; x<row-1; x++){
                val = valVec[i + (x+1)];
                A[row-1][((int)innerBandwith)-1 - i + (x+1)] =  val;
                //cout << "before bw wrote " << val <<endl;
            }
            i+=x;
        }
        // soldan saga doldurmaya baslayabilirsin artik.
        else{
            //cout << "finished first part" << " " ;
            A[row-1][counter] =  val;
            //cout << "wrote " << val <<endl;
            counter++;
        }
       // cout << i << " ";
    }
    cout << "algo finished." << endl;

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
