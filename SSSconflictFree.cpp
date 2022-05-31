//
// Created by Selin Yıldırım on 7.01.2022.
//

#include <iostream>
#include  "header.h"
//#include <mpi.h>
//#include "rcmtest.cpp"
//#include "geeks.cpp"

using namespace std;
#define MATRIX_COUNT 6

int readSSSFormat(int z, bool banded) {
    fs::path matrixFolder;
    if(!banded) matrixFolder = "/home/selin/SSS-Data/" + matrix_names[z] + "/unbanded" ;
    else matrixFolder = "/home/selin/SSS-Data/" + matrix_names[z];
    for(auto const& dir_entry: fs::directory_iterator{matrixFolder}){
        std::fstream myfile(dir_entry.path(), std::ios_base::in);
        if(dir_entry.path().stem() == "rowptr") {
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
        else if(dir_entry.path().stem() == "colind") {
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
        else if(dir_entry.path().stem() == "diag"){
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
        else if(dir_entry.path().stem() == "vals"){
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
        else cout << "unexpected file name: " << dir_entry.path() << endl;
    }
    return 0;
}
/*
    A RNxN : matrix in SSS format
    x RN : input vector
    y RN : output vector

 */

int main(int argc, char **argv){
    int n, rowLimit;
    cout << "i call read CSR Format. " << endl;
    //init();
    if(!argv[1]){
        cout << "please provide input matrix index (int): boneS10, Emilia_923, ldoor, af_5_k101, Serena, audikw_1" << endl;
        return -1;
    }
    if(!argv[2]){
        cout << "please provide bool for banded" << endl;
        return -1;
    }
    bool banded = atoi(argv[2]);
    readSSSFormat(atoi(argv[1]) , banded);

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
    //memset(y, 0, (size_t) (n*sizeof(double)));
    cout << x[0] << " " << y[0];

    int colInd;
    cout << "start computing sequential ssbmv..." << endl;
    for(int i=0; i<n; i++){
        y[i] = matrixDiagonal[i] * x[i];
        for(int j=matrixRowptr[i]-1; j<matrixRowptr[i+1]-1; j++){
            colInd = matrixColind[j] - 1;
            y[i] += matrixOffDiagonal[j] * x[colInd];
            y[colInd] += matrixOffDiagonal[j]*x[i];
        }
    }
    cout << "finished computing sequential ssbmv." << endl;


    ofstream myfile1;
    if(!banded) myfile1.open ("/home/selin/Seq-Results/" + matrix_names[inputType] + "/unbanded/result.txt", ios::out | ios::trunc);
    else myfile1.open ("/home/selin/Seq-Results/" + matrix_names[inputType] + "/banded/result.txt", ios::out | ios::trunc);

    cout << "Writing to output. " << endl;
    for (int i=0; i<n; i++) {
        myfile1 << y[i] << '\t';
    }
    myfile1.close();

    return 0;
}
