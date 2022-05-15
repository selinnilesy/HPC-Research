//
// Created by Selin Yıldırım on 7.01.2022.
//

#include <iostream>
#include <string>
#include  "header.h"
//#include <mpi.h>
//#include "rcmtest.cpp"
//#include "geeks.cpp"

using namespace std;
#define MATRIX_COUNT 6
vector<int> col_inner, row_inner;
vector<int> col_outer, row_outer;
vector<int> col_middle, row_middle;
vector<double> vals_inner, vals_outer, vals_middle;

int p;
vector<int> x;
vector<int> y;
vector<int> elm_row; // use this to find out row_ptr values
int *coord_row, *coord_col;
double *coord_val;

void init(){
    double *values, *dvalues;
    int *colind, *rowptr;
    valuesPtrs.push_back(values);
    dvaluesPtrs.push_back(dvalues);
    colindPtrs.push_back(colind);
    rowptrPtrs.push_back(rowptr);
}

int readSSSFormat(int z) {
    double tempVal;
    vector<double> tempVec;

    /*
    fstream myfile("/home/selin/CSR-Data/" + matrix_names[z] + "/banded/CSRout_col.txt", std::ios_base::in);
    int x;
    for(int i=0; i<73; i++){
        myfile >> x;
        cout << x << '\t' ;
    }
    myfile.close();
     */

    const fs::path matrixFolder{"/home/selin/SSS-Data/" + matrix_names[z]};
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

            cout << dir_entry.path() << " has been read." << endl;
            myfile.close();
            continue;
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
            cout << dir_entry.path() << " has been read." << endl;
            myfile.close();
            continue;
        }
        // else, start reading doubles.
        while (myfile >> tempVal) {
            tempVec.push_back(tempVal);
        }

        if(dir_entry.path().stem() == "diag"){
            dvaluesPtrs.push_back(new double[tempVec.size()]);
            double *temp = dvaluesPtrs[0];
            for(int i=0; i<tempVec.size(); i++) temp[i]=tempVec[i];
            dvaluesSize.push_back(tempVec.size());
        }
        else if(dir_entry.path().stem() == "vals"){
            valuesPtrs.push_back(new double[tempVec.size()]);
            double *temp = valuesPtrs[0];
            for(int i=0; i<tempVec.size(); i++) temp[i]=tempVec[i];
            valuesSize.push_back(tempVec.size());
        }
        else cout << "unexpected file name: " << dir_entry.path() << endl;
        cout << dir_entry.path() << " has been read." << endl;

        tempVec.clear();
        myfile.close();
    }
    return 0;
}
int writeCooFormat(int z, double inputRatio) {
    string dirpath="/home/selin/Split-Data/" + matrix_names[z];
    // row index
    ofstream myfile(dirpath + "/inner/coordinate-"+ to_string(inputRatio)+ "-row.txt", std::fstream::out);
    for(int i=0; i<row_inner.size(); i++) myfile << row_inner[i] << '\t' ;
    cout << "inner/coordinate-row.txt" << " has been written." << endl;
    myfile.close();
    // col index
    myfile.open(dirpath + "/inner/coordinate-"+ to_string(inputRatio) + "-col.txt", std::fstream::out);
    for(int i=0; i<col_inner.size() ;i++) myfile << col_inner[i] << '\t' ;
    cout << "inner/coordinate-col.txt" << " has been written." << endl;
    myfile.close();
    // vals
    myfile.open(dirpath + "/inner/coordinate-"+ to_string(inputRatio) + "-val.txt", std::fstream::out);
    for(int i=0; i<vals_inner.size(); i++) myfile << vals_inner[i] << '\t' ;
    cout << "inner/coordinate-val.txt" << " has been written." << endl;
    myfile.close();
    // -----------
    // row index
    myfile.open(dirpath + "/middle/coordinate-"+ to_string(inputRatio) + "-row.txt", std::fstream::out);
    for(int i=0; i<row_middle.size() ;i++) myfile << row_middle[i] << '\t' ;
    cout << "middle/coordinate-row.txt" << " has been written." << endl;
    myfile.close();
    // col index
    myfile.open(dirpath + "/middle/coordinate-"+ to_string(inputRatio) + "-col.txt", std::fstream::out);
    for(int i=0; i<col_middle.size(); i++) myfile << col_middle[i] << '\t' ;
    cout << "middle/coordinate-col.txt"<< " has been written." << endl;
    myfile.close();
    // vals
    myfile.open(dirpath + "/middle/coordinate-"+ to_string(inputRatio) + "-val.txt", std::fstream::out);
    for(int i=0; i<vals_middle.size(); i++) myfile << vals_middle[i] << '\t' ;
    cout << "middle/coordinate-val.txt" << " has been written." << endl;
    myfile.close();
    // -----------
    // row index
    myfile.open(dirpath + "/outer/coordinate-"+ to_string(inputRatio) + "-row.txt", std::fstream::out);
    for(int i=0; i<row_outer.size(); i++) myfile << row_outer[i] << '\t' ;
    cout << "outer/coordinate-row.txt" << " has been written." << endl;
    myfile.close();
    // col index
    myfile.open(dirpath + "/outer/coordinate-"+ to_string(inputRatio) + "-col.txt", std::fstream::out);
    for(int i=0; i<col_outer.size(); i++) myfile << col_outer[i] << '\t' ;
    cout << "outer/coordinate-col.txt" << " has been written." << endl;
    myfile.close();
    // vals
    myfile.open(dirpath + "/outer/coordinate-"+ to_string(inputRatio) + "val.txt", std::fstream::out);
    for(int i=0; i<vals_outer.size(); i++) myfile << vals_outer[i] << '\t' ;
    cout << "outer/coordinate-val.txt" << " has been written." << endl;
    myfile.close();

    return 0;
}
int main(int argc, char **argv){
    int n, rowLimit;
    cout << "i call readSSSFormat. " << endl;
    //init();
    if(!argv[1]){
        cout << "please provide input matrix index (int): boneS10, Emilia_923, ldoor, af_5_k101, Serena, audikw_1" << endl;
        return -1;
    }
    if(!argv[2]){
        cout << "please provide a ratio for bandwiths" << endl;
        return -1;
    }
    readSSSFormat(atoi(argv[1]));

    n = matrixSize[atoi(argv[1])];
    int inputType = atoi(argv[1]);
    double inputRatio = atof(argv[2]);
    cout << "input ratio: " << inputRatio << endl;

    double *matrixOffDiagonal = valuesPtrs[0];
    int *matrixColind = colindPtrs[0];
    int *matrixRowptr= rowptrPtrs[0];

    int elmCountPerRow, colInd, rowBegin;
    double innerBandwith, middleBandwith;
    innerBandwith = (int) (nnz_n_Ratios[inputType]*bandwithProportions[inputType] * inputRatio);
    middleBandwith = bandwithSize[inputType] - 2*innerBandwith;
    cout << "inner bandwith: " << innerBandwith << endl;
    cout << "middle bandwith: " << middleBandwith << endl;
    cout << "total bandwith: " << bandwithSize[inputType] << endl;
    int maxJ=-1;
    double val;

    if(rowptrSize[0] != (matrixSize[inputType]+1)) {
        cout << "Corrupt rowptr. read size does not match original matrix rowptr size." << endl;
        return -1;
    }

    // two new SSS storage for inner and outer regions

    int counter_inner, counter_middle, counter_outer;

    for (int i = 0; i < rowptrSize[0] - 1; i++) {
        // row ptrs start from 1 !!!
        rowBegin = matrixRowptr[i] - 1;
        elmCountPerRow = matrixRowptr[i + 1] - matrixRowptr[i];
        maxJ=-1;
        counter_inner=counter_middle=counter_outer = 0;
        for (int j = 0; j < elmCountPerRow; j++) {
            // i = row indexi
            // colind = column indexi
            // val = degeri
            colInd = matrixColind[rowBegin + j] - 1 ;
            val = matrixOffDiagonal[rowBegin + j];

            if(colInd > maxJ) maxJ=colInd;
            else cout << "ON THE SAME ROW, COL INDEX HAS BEEN SMALLED. NOT IN ASCENDING ORDER: maxj, colInd " << maxJ <<" " << colInd << endl;
            // inner Dense Region
            if(colInd >= i - innerBandwith){
                row_inner.push_back(i+1);
                col_inner.push_back(colInd+1);
                vals_inner.push_back(val);
            }
            // middle Region
            else if(colInd >= i-innerBandwith-middleBandwith){
                row_middle.push_back(i+1);
                col_middle.push_back(colInd+1);
                vals_middle.push_back(val);
            }
            // outer Dense Region
            else{
                row_outer.push_back(i+1);
                col_outer.push_back(colInd+1);
                vals_outer.push_back(val);
            }
        }
    }
    cout << "read SSS format and grouped 3way bandwith" << endl;
    writeCooFormat(inputType, inputRatio);



    delete [] matrixOffDiagonal;
    delete [] matrixColind;
    delete [] matrixRowptr;
}