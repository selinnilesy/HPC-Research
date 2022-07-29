//
// Created by Selin Yıldırım on 7.01.2022.
//

#include <iostream>
#include <omp.h>
#include  "header.h"
#include <limits>
#include <iomanip>
//#include <mpi.h>
//#include "rcmtest.cpp"
//#include "geeks.cpp"

using namespace std;
#define MATRIX_COUNT 6

vector<int> colind_up, rowptr_up;
vector<double> offdiag_up;

vector<int> outer_col, outer_row;
vector<double> outer_val;

typedef std::numeric_limits< double > dbl;

int readSSSFormat(int z) {
    fs::path matrixFolder;
     matrixFolder = "/home/selin/SSS-Data/" + matrix_names[z];
    for(auto const& dir_entry: fs::directory_iterator{matrixFolder}){
        std::fstream myfile(dir_entry.path(), std::ios_base::in);
        if(dir_entry.path().stem() == "diag"){
            double tempVal;
            vector<double> tempVec;
            while (myfile >> tempVal) {
                tempVec.push_back(tempVal);
            }
            dvaluesPtrs.push_back(new double[tempVec.size()]);
            double *temp = dvaluesPtrs[0];
            for(int i=0; i<tempVec.size(); i++) temp[i]=tempVec[i];
            dvaluesSize.push_back(tempVec.size());
            cout << dir_entry.path() << " has been read with size: " <<  tempVec.size() << endl;
            myfile.close();
        }
        //else cout << "unused file name: " << dir_entry.path() << endl;
    }
    return 0;
}
int readCSRFormat(int z) {
    fs::path matrixFolder;
    matrixFolder = "/home/selin/Split-Data/" + matrix_names[z] + "/middle/CSR-Data";
    for(auto const& dir_entry: fs::directory_iterator{matrixFolder}){
        std::fstream myfile(dir_entry.path(), std::ios_base::in);
        if(dir_entry.path().stem() == ("1.000000-" + to_string((double) bandwithSize[z]-3) + "-row")) {
            int tempValInt;
            vector<int> tempVecInt;
            while (myfile >> tempValInt) {
                tempVecInt.push_back(tempValInt);
            }
            rowptrPtrs.push_back(new int[tempVecInt.size()]);
            int *temp = rowptrPtrs[0];
            for(int i=0; i<tempVecInt.size(); i++) temp[i]=tempVecInt[i];
            rowptrSize.push_back(tempVecInt.size());
            cout << dir_entry.path() << " has been read with size: " <<tempVecInt.size() << endl;
            myfile.close();
        }
        else if(dir_entry.path().stem() == ("1.000000-" + to_string( (double) bandwithSize[z]-3) + "-col")) {
            int tempValInt;
            vector<int> tempVecInt;
            while (myfile >> tempValInt) {
                tempVecInt.push_back(tempValInt);
            }
            colindPtrs.push_back(new int[tempVecInt.size()]);
            int *temp = colindPtrs[0];
            for(int i=0; i<tempVecInt.size(); i++) temp[i]=tempVecInt[i];
            colindSize.push_back(tempVecInt.size());
            cout << dir_entry.path() << " has been read with size: " <<  tempVecInt.size() << endl;
            myfile.close();
        }
        else if(dir_entry.path().stem() == ("1.000000-" + to_string((double) bandwithSize[z]-3) + "-val")){
            double tempVal;
            vector<double> tempVec;
            while (myfile >> tempVal) {
                tempVec.push_back(tempVal);
            }
            valuesPtrs.push_back(new double[tempVec.size()]);
            double *temp = valuesPtrs[0];
            for(int i=0; i<tempVec.size(); i++) temp[i]=tempVec[i];
            valuesSize.push_back(tempVec.size());
            cout << dir_entry.path() << " has been read with size: " <<  tempVec.size() << endl;
            myfile.close();
        }
        //else cout << "unexpected file name: " << dir_entry.path() << endl;
    }

        //else cout << "unexpected file name: " << dir_entry.path() << endl;
        matrixFolder = "/home/selin/Split-Data/" + matrix_names[z] + "/outer/";
        for(auto const& dir_entry: fs::directory_iterator{matrixFolder}) {
            std::fstream myfile(dir_entry.path(), std::ios_base::in);
            if (dir_entry.path().stem() == ("1.000000-" + to_string((double) bandwithSize[z] - 3) + "-row")) {
                int tempValInt;
                while (myfile >> tempValInt) {
                    outer_row.push_back(tempValInt);
                }
                cout << dir_entry.path() << " has been read with size: " << outer_row.size() << endl;
                myfile.close();
            } else if (dir_entry.path().stem() == ("1.000000-" + to_string((double) bandwithSize[z] - 3) + "-col")) {
                int tempValInt;
                while (myfile >> tempValInt) {
                    outer_col.push_back(tempValInt);
                }
                cout << dir_entry.path() << " has been read with size: " << outer_col.size() << endl;
                myfile.close();
            } else if (dir_entry.path().stem() == ("1.000000-" + to_string((double) bandwithSize[z] - 3) + "val")) {
                double tempVal;
                while (myfile >> tempVal) {
                    outer_val.push_back(tempVal);
                }
                cout << dir_entry.path() << " has been read with size: " << outer_val.size() << endl;
                myfile.close();
            }
        }
    return 0;
}
/*
    A RNxN : matrix in SSS format
    x RN : input vector
    y RN : output vector

 */

int main(int argc, char **argv){
    int n;
    //init();
    if(!argv[1]){
        cout << "please provide input matrix index (int): boneS10, Emilia_923, ldoor, af_5_k101, Serena, audikw_1" << endl;
        return -1;
    }

    readSSSFormat(atoi(argv[1]));
    readCSRFormat(atoi(argv[1]));

    n = matrixSize[atoi(argv[1])];
    int inputType = atoi(argv[1]);


    double *matrixOffDiagonal = valuesPtrs[0];
    double *matrixDiagonal = dvaluesPtrs[0];
    int *matrixColind = colindPtrs[0];
    int *matrixRowptr= rowptrPtrs[0];
    double *x = new double[n];
    double *y = new double[n];
    for(int i=0; i<n; i++) x[i] = 1.0;
    for(int i=0; i<n; i++) y[i] = 0.0;
    ofstream myfile1;

    int colInd;
    cout << "start computing serial SSS mv..." << endl;
    double row_i,row_e, val;
    double t = omp_get_wtime();
        // middle -  lower
        for (int i = 0; i < n; i++) {
            val = matrixDiagonal[i] * x[i];
            row_i = matrixRowptr[i] - 1;
            row_e = matrixRowptr[i + 1] - 1;
            for (int j = row_i; j < row_e; j++) {
                colInd = matrixColind[j] - 1;
                // skew-symm
                val += matrixOffDiagonal[j] * x[colInd];
                //if(i==686172) cout << "accumulating on 238050 - offdiag: " << matrixOffDiagonal[j] << " x: " << x[colInd] << endl;
                // middle -  upper
                y[colInd] -= matrixOffDiagonal[j] * x[i];
                //if(colInd==686172)  cout << "adding colInd to 686171: " <<  -matrixOffDiagonal[j] * x[i] << endl;
            }
            y[i] += val;
        }
        for (int i = 0; i < outer_row.size(); i++) {
            // outer - lower
            y[outer_row[i] - 1] += outer_val[i] * x[outer_col[i] - 1];
            // outer - upper
            y[outer_col[i] - 1] -= outer_val[i] * x[outer_row[i] - 1];
        }
    t = omp_get_wtime() - t;
    printf ("It took me %f seconds for 1000-times serial run.\n", t);
    myfile1.open ("/home/selin/3way-Seq-Results/" + matrix_names[inputType] + "/result.txt", ios::out | ios::trunc);


    cout << "Writing to output... " << endl;
    for (int i=0; i<n; i++) {
        myfile1 << std::fixed << std::setprecision(dbl::max_digits10) << y[i] << '\t';
    }
    myfile1.close();
    cout << "Completed output... " << endl;


    delete [] x;
    delete [] y;
    delete [] matrixRowptr;
    delete [] matrixColind;
    delete [] matrixOffDiagonal;
    delete [] matrixDiagonal;
    return 0;
}