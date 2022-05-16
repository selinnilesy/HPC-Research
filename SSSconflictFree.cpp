#include <cblas.h>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <string>
#include <cmath>
#include "header.h"

using namespace std;
vector<int> inner_colVec;
vector<int> inner_rowVec;
vector<double> inner_valVec;
vector<double> inner_diagVec;

vector<int> middle_colVec;
vector<int> middle_rowVec;
vector<double> middle_valVec;
vector<double> middle_diagVec;

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
                inner_rowVec.push_back(intVal);
            }
            myfile.close();
            cout << rowFileName << " has been read with size: " << inner_rowVec.size() << endl;

            myfile.open(colFileName, std::ios_base::in);
            while (myfile >> intVal) {
                inner_colVec.push_back(intVal);
            }
            myfile.close();
            cout << colFileName << " has been read with size: " << inner_colVec.size() << endl;

            myfile.open(valFileName, std::ios_base::in);
            while (myfile >> doubleVal) {
                inner_valVec.push_back(doubleVal);
            }
            myfile.close();
            cout << valFileName << " has been read with size: " << inner_valVec.size() << endl;
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
        inner_diagVec.push_back(doubleVal);
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
     * [1 0 0 0 0  0 0 0    ]
     * [8 1 0 0 0  0 0 0    ]
     * [9 5 1 0 0  0 0 0    ]
     * [5 7 6 1 0  0 0 0   ]
     * [3 5 4 3 1  0 0 0   ]
     * [0 3 5 3 4  1 0 0   ]
     * [0 0 3 5 9  9 1 0   ]
     * [0 0 0 3 5  8 2 1   ]
     *

    rowVec.push_back(2);
    rowVec.push_back(3);

    rowVec.push_back(3);
    rowVec.push_back(4);
    rowVec.push_back(4);
    rowVec.push_back(4);

    rowVec.push_back(5);
    rowVec.push_back(5);
    rowVec.push_back(5);
    rowVec.push_back(5);

    rowVec.push_back(6);
    rowVec.push_back(6);
    rowVec.push_back(6);
    rowVec.push_back(6);

    rowVec.push_back(7);
    rowVec.push_back(7);
    rowVec.push_back(7);
    rowVec.push_back(7);

    rowVec.push_back(8);
    rowVec.push_back(8);
    rowVec.push_back(8);
    rowVec.push_back(8);

    // ---
    valVec.push_back(8);
    valVec.push_back(9);
    valVec.push_back(5);
    valVec.push_back(5);

    valVec.push_back(7);
    valVec.push_back(6);
    valVec.push_back(3);
    valVec.push_back(5);

    valVec.push_back(4);
    valVec.push_back(3);
    valVec.push_back(3);
    valVec.push_back(5);

    valVec.push_back(3);
    valVec.push_back(4);
    valVec.push_back(3);
    valVec.push_back(5);

    valVec.push_back(9);
    valVec.push_back(9);
    valVec.push_back(3);
    valVec.push_back(5);

    valVec.push_back(8);
    valVec.push_back(2);
 */


    if(inner_rowVec.size() != inner_valVec.size()){
        cout << "not equal size - rowVec and valVec (NNZ) !!! " << endl;
        return -1;
    }
    if(n != inner_diagVec.size()){
        cout << "not equal size - diagVec (n) !!! " << endl;
        return -1;
    }


    double innerBandwith, middleBandwith;
    innerBandwith = (int) (nnz_n_Ratios[inputType]*bandwithProportions[inputType] * inputRatio);
    middleBandwith = bandwithSize[inputType] - 2*innerBandwith;
    //innerBandwith=4;
    //middleBandwith=0;
    cout << "inner bandwith: " << innerBandwith << endl;
    cout << "middle bandwith: " << middleBandwith << endl;
    //cout << "total bandwith: " << bandwithSize[inputType] << endl;

    int i,j;
    int size= n;
    //int size=8;
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

    float** A_middle = new float*[size];
    for( i=0; i<size; i++) {
        A_middle[i]  = new float[lda];
    }

    for( i=0; i<size; i++) {
        for( j=0; j<lda; j++) {
            A_middle[i][j]  = 0.0;
        }
    }
    //memset(A, 0, k*k);


    cout << "A initialized to 0." <<  endl;
    int row, col, val, counter=0, x;
    // i keeps track of whole element count.
    // upper bound is not n, but instead n-1 !
    cout << "start generating banded for inner" <<  endl;
    for( i=0; i<inner_rowVec.size(); i++) {
        row = inner_rowVec[i] - 1;
        //col = colVec[i] - 1;
        val = inner_valVec[i];
        //cout << "row: " << row << " col: " << col <<  " val: " << val << endl;
        if(counter==innerBandwith) counter=0;

        if(row <= innerBandwith){
            // insert first element
            A[row][((int)innerBandwith) - row] =  val;
            //cout << "inserted: " << val << endl;
            for(x=0; x<row-1; x++){
                val = inner_valVec[i + (x+1)];
                A[row][((int)innerBandwith) - row + (x+1)] =  val;
                //cout << "before bw wrote " << val <<endl;
            }
            i+=x;
        }
        // soldan saga doldurmaya baslayabilirsin artik.
        else{
            //cout << "finished first part" << " " ;
            A[row][counter] =  val;
            //cout << "wrote " << val <<endl;
            counter++;
        }
       // cout << i << " ";
    }
    cout << "start generating banded for middle" <<  endl;
    x=counter=0;
    for( i=0; i<middle_rowVec.size(); i++) {
        row = middle_rowVec[i] - innerBandwith - 1;
        //col = colVec[i] - 1;
        val = middle_valVec[i];
        //cout << "row: " << row << " col: " << col <<  " val: " << val << endl;
        if(counter==middleBandwith) counter=0;

        if(row+innerBandwith <= middleBandwith){
            // insert first element
            A[row-1][((int)middleBandwith) - row] =  val;
            //cout << "inserted: " << val << endl;
            for(x=0; x<row-1; x++){
                val = middle_valVec[i + (x+1)];
                A[row-1][((int)middleBandwith) - row + (x+1)] =  val;
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

    for(i=0; i<inner_diagVec.size(); i++){
        A[i][(int) innerBandwith] = inner_diagVec[i];
    }
    cout << "written diag onto A." << endl;

    ofstream myfile;
    string output =  "/home/selin/Split-Data/"+ matrix_names[inputType]+ "/inner-banded-A" + to_string(inputRatio)+ ".txt";
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
