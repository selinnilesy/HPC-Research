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

int p;
vector<int> x;
vector<int> y;
vector<int> elm_row; // use this to find out row_ptr values
int *coord_row, *coord_col;
double *coord_val;
int nonzeroSize_row, nonzeroSize_col, nonzeroSize_val;



void init(){
    double *values, *dvalues;
    int *colind, *rowptr;
    valuesPtrs.push_back(values);
    dvaluesPtrs.push_back(dvalues);
    colindPtrs.push_back(colind);
    rowptrPtrs.push_back(rowptr);
}

/*extern "C" {
    extern void i4_swap_ ( int *i, int *j );
    extern void  reverse_ ( int *n, int *a );
    extern void degree_  ( int root, int adj_num, int* adj_row, int* adj, int mask, int* deg, int *iccsze, int* ls, int *node_num );
    extern void rcm_  ( int root, int adj_num, int* adj_row, int* adj, int mask, int* perm, int *iccsze, int *node_num );
}*/

extern "C" {
extern void coocsr_(int *nrow, int *nnz, double* a, int* ir,int* jc, double* ao, int* jao, int* iao);
}

int readCooFormatEqual(int z, double inputRatio, double middleRatio, bool reversed) {
    double tempVal;
    int tempValInt;
    vector<double> tempVec;
    fstream myfile;
    const fs::path matrixFolder{"/home/selin/Split-Data/" + matrix_names[z] + "/inner"};
    for(auto const& dir_entry: fs::directory_iterator{matrixFolder}){
        if(dir_entry.path().stem() == ("coordinate-"+ to_string(inputRatio)+"-" + to_string(middleRatio) + "-row")) {
            myfile.open(dir_entry.path(), std::ios_base::in);
            vector<int> tempVecInt;
            while (myfile >> tempValInt) {
                tempVecInt.push_back(tempValInt);
            }
            coord_row = new int[tempVecInt.size()];
            for(int i=0; i<tempVecInt.size(); i++) coord_row[i]=tempVecInt[i];
            nonzeroSize_row=tempVecInt.size();
            cout << dir_entry.path() << " has been read." << endl;
            myfile.close();
        }
        else if(dir_entry.path().stem() == ("coordinate-"+ to_string(inputRatio)+"-" + to_string(middleRatio) + "-col")) {
            myfile.open(dir_entry.path(), std::ios_base::in);
            vector<int> tempVecInt;
            while (myfile >> tempValInt) {
                tempVecInt.push_back(tempValInt);
            }
            coord_col = new int[tempVecInt.size()];
            for(int i=0; i<tempVecInt.size(); i++) coord_col[i]=tempVecInt[i];
            nonzeroSize_col=tempVecInt.size();
            cout << dir_entry.path() << " has been read." << endl;
            myfile.close();
        }
        else if(dir_entry.path().stem() == ("coordinate-"+ to_string(inputRatio)+ "-" + to_string(middleRatio) + "-val")){
            myfile.open(dir_entry.path(), std::ios_base::in);
            while (myfile >> tempVal) {
                tempVec.push_back(tempVal);
            }
            coord_val = new double[tempVec.size()];
            // !reversed means lower and reverse means upper.
            if(!reversed){
                for(int i=0; i<tempVec.size(); i++) coord_val[i] = -tempVec[i];
            }
            else if(reversed){
                // use neg for skew-symmetric part !!!
                for(int i=0; i<tempVec.size(); i++) coord_val[i] = tempVec[i];
            }
            nonzeroSize_val=tempVec.size();
            cout << dir_entry.path() << " has been read." << endl;
            myfile.close();
        }
    }
    return 0;
}
int readCooFormatNotEqual(int z, double inputRatio,double middleRatio,  bool reversed) {
    double tempVal;
    int tempValInt;
    vector<double> tempVec;
    fstream myfile;
    const fs::path matrixFolder{"/home/selin/Split-Data/" + matrix_names[z] + "/middle/"};
    for(auto const& dir_entry: fs::directory_iterator{matrixFolder}){
        if(dir_entry.path().stem() == ("coordinate-"+ to_string(inputRatio)+ "-" + to_string(middleRatio)+ "-row")) {
            myfile.open(dir_entry.path(), std::ios_base::in);
            vector<int> tempVecInt;
            while (myfile >> tempValInt) {
                tempVecInt.push_back(tempValInt);
            }
            coord_row = new int[tempVecInt.size()];
            for(int i=0; i<tempVecInt.size(); i++) coord_row[i]=tempVecInt[i];
            nonzeroSize_row=tempVecInt.size();
            cout << dir_entry.path() << " has been read." << endl;
            myfile.close();
        }
        else if(dir_entry.path().stem() == ("coordinate-"+ to_string(inputRatio) + "-" + to_string(middleRatio) + "-col")) {
            myfile.open(dir_entry.path(), std::ios_base::in);
            vector<int> tempVecInt;
            while (myfile >> tempValInt) {
                tempVecInt.push_back(tempValInt);
            }
            coord_col = new int[tempVecInt.size()];
            for(int i=0; i<tempVecInt.size(); i++) coord_col[i]=tempVecInt[i];
            nonzeroSize_col=tempVecInt.size();
            cout << dir_entry.path() << " has been read." << endl;
            myfile.close();
        }
        else if(dir_entry.path().stem() == ("coordinate-"+ to_string(inputRatio) + "-" + to_string(middleRatio) + "-val")){
            myfile.open(dir_entry.path(), std::ios_base::in);
            while (myfile >> tempVal) {
                tempVec.push_back(tempVal);
            }
            coord_val = new double[tempVec.size()];
            // !reversed means lower and reverse means upper.
            if(!reversed){
                for(int i=0; i<tempVec.size(); i++) coord_val[i] = -tempVec[i];
            }
            else if(reversed){
                // use neg for skew-symmetric part !!!
                for(int i=0; i<tempVec.size(); i++) coord_val[i] = tempVec[i];
            }
            nonzeroSize_val=tempVec.size();
            cout << dir_entry.path() << " has been read." << endl;
            myfile.close();
        }
    }
    return 0;
}

int main(int argc, char **argv){

    int n, rowLimit;
    //init();
    if(!argv[1]){
        cout << "please provide input matrix index (int): boneS10, Emilia_923, ldoor, af_5_k101, Serena, audikw_1" << endl;
        return -1;
    }
    if(!argv[2]){
        cout << "please provide input ratio for inner bandwith" << endl;
        return -1;
    }
    if(!argv[3]){
        cout << "please provide middle bandwith" << endl;
        return -1;
    }
    if(!argv[4]){
        cout << "please provide a boolean for reversed (process only upper corresponding elements)" << endl;
        return -1;
    }
    if(!argv[5]){
        cout << "please provide a boolean for inner-equal (1) or middle (0) directories" << endl;
        return -1;
    }
    bool reversed = atoi(argv[4]);
    bool inner_equal_middle = atoi(argv[5]);
    double middleRatio = atof(argv[3]);
    n = matrixSize[atoi(argv[1])];
    int inputType = atoi(argv[1]);
    double inputRatio =  atof(argv[2]);

    inputRatio = 1;
    middleRatio = bandwithSize[inputType]-20;

    std::cout  <<  "reversed ?: " << reversed << endl;
    std::cout  <<  "inner_equal_middle ?: " << inner_equal_middle << endl;
    if(!inner_equal_middle) readCooFormatNotEqual(inputType,   inputRatio, middleRatio ,reversed);
    else if(inner_equal_middle)  readCooFormatEqual(inputType,  inputRatio,middleRatio,  reversed);

    if(nonzeroSize_row != nonzeroSize_col) std::cout  <<  " nonzeroSize_row and nonzeroSize_col not equal"  << endl;
    if(nonzeroSize_row != nonzeroSize_val) std::cout  <<  " nonzeroSize_row and nonzeroSize_val not equal"  << endl;
    if(nonzeroSize_col != nonzeroSize_val) std::cout  <<  " nonzeroSize_col and nonzeroSize_val not equal"  << endl;

    int *banded_coordRow = coord_row;
    int *banded_coordCol = coord_col;
    double *banded_coordval = coord_val;

    int nrow=matrixSize[inputType];
    int nnz = nonzeroSize_row ;

    int *banded_csrRow = new int[matrixSize[inputType]+1];
    int *banded_csrCol = new int[nnz];
    double *banded_csrval = new double[nnz];


    std::cout  <<  " starts computing coocsr... " << endl;
    // to be able to compare results with cblas_ssbmv, lower part values are negative.
    if(!reversed) coocsr_(&nrow,  &nnz, banded_coordval, banded_coordRow, banded_coordCol, banded_csrval, banded_csrCol, banded_csrRow);
    else if(reversed)  coocsr_(&nrow,  &nnz, banded_coordval, banded_coordCol, banded_coordRow, banded_csrval, banded_csrCol, banded_csrRow);
    std::cout  <<  " finished computing coocsr... :" << banded_csrval[10] << " " << banded_csrCol[10] << " " << banded_csrRow[10] << endl;

    ofstream myfile1, myfile2, myfile3;
    if(inner_equal_middle) {
        if (reversed) {
            myfile1.open(
                    "/home/selin/Split-Data/" + matrix_names[inputType] + "/inner/CSR-Data/upper/" +
                    to_string(inputRatio) + "-row.txt", ios::out | ios::trunc);
            myfile2.open(
                    "/home/selin/Split-Data/" + matrix_names[inputType] + "/inner/CSR-Data/upper/" +
                    to_string(inputRatio) + "-col.txt", ios::out | ios::trunc);
            myfile3.open(
                    "/home/selin/Split-Data/" + matrix_names[inputType] + "/inner/CSR-Data/upper/" +
                    to_string(inputRatio) + "-val.txt", ios::out | ios::trunc);
        } else if (!reversed) {
            myfile1.open("/home/selin/Split-Data/" + matrix_names[inputType] + "/inner/CSR-Data/" +
                         to_string(inputRatio) + "-row.txt", ios::out | ios::trunc);
            myfile2.open("/home/selin/Split-Data/" + matrix_names[inputType] + "/inner/CSR-Data/" +
                         to_string(inputRatio) + "-col.txt", ios::out | ios::trunc);
            myfile3.open("/home/selin/Split-Data/" + matrix_names[inputType] + "/inner/CSR-Data/" +
                         to_string(inputRatio) + "-val.txt", ios::out | ios::trunc);
        }
    }
    else{
        if (reversed) {
            myfile1.open(
                    "/home/selin/Split-Data/" + matrix_names[inputType] + "/middle/CSR-Data/upper/" +
                    to_string(inputRatio)+ "-" + to_string(middleRatio) + "-row.txt", ios::out | ios::trunc);
            myfile2.open(
                    "/home/selin/Split-Data/" + matrix_names[inputType] + "/middle/CSR-Data/upper/" +
                    to_string(inputRatio)+ "-" + to_string(middleRatio) + "-col.txt", ios::out | ios::trunc);
            myfile3.open(
                    "/home/selin/Split-Data/" + matrix_names[inputType] + "/middle/CSR-Data/upper/" +
                    to_string(inputRatio)+ "-" + to_string(middleRatio) + "-val.txt", ios::out | ios::trunc);
        } else if (!reversed) {
            myfile1.open("/home/selin/Split-Data/" + matrix_names[inputType] + "/middle/CSR-Data/" +
                         to_string(inputRatio)+ "-" + to_string(middleRatio) + "-row.txt", ios::out | ios::trunc);
            myfile2.open("/home/selin/Split-Data/" + matrix_names[inputType] + "/middle/CSR-Data/" +
                         to_string(inputRatio)+ "-" + to_string(middleRatio) + "-col.txt", ios::out | ios::trunc);
            myfile3.open("/home/selin/Split-Data/" + matrix_names[inputType] + "/middle/CSR-Data/" +
                         to_string(inputRatio)+ "-" + to_string(middleRatio) + "-val.txt", ios::out | ios::trunc);
        }
    }
    cout << "Writing to " << "-CSR row.txt"  << endl;
    for (int i=0; i<nrow+1; i++) {
        myfile1 << banded_csrRow[i] << '\t';
    }
    cout << "Writing to " << "-CSR col.txt"  << endl;
    for (int i=0; i<nnz; i++) {
        myfile2 << banded_csrCol[i] << '\t';
    }
    cout << "Writing to " << "-CSR val.txt"  << endl;
    for (int i=0; i<nnz; i++){
        myfile3 << banded_csrval[i] << '\t';
    }
    myfile1.close();
    myfile2.close();
    myfile3.close();
}