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

/*extern "C" {
    extern void i4_swap_ ( int *i, int *j );
    extern void  reverse_ ( int *n, int *a );
    extern void degree_  ( int root, int adj_num, int* adj_row, int* adj, int mask, int* deg, int *iccsze, int* ls, int *node_num );
    extern void rcm_  ( int root, int adj_num, int* adj_row, int* adj, int mask, int* perm, int *iccsze, int *node_num );
}*/

extern "C" {
extern void coocsr_(int *nrow, int *nnz, double* a, int* ir,int* jc, double* ao, int* jao, int* iao);
}

int readCooFormat(int z) {
    double tempVal;
    vector<double> tempVec;
    const fs::path matrixFolder{"/home/selin/Coo-Data/unbanded/" + matrix_names[z]};
    for(auto const& dir_entry: fs::directory_iterator{matrixFolder}){
        std::fstream myfile(dir_entry.path(), std::ios_base::in);
        if(dir_entry.path().stem() == "coordinate_row") {
            int tempValInt;
            vector<int> tempVecInt;
            while (myfile >> tempValInt) {
                tempVecInt.push_back(tempValInt);
            }
            if(tempVecInt.size() != nonzerosSize[z]) cout << "Row count not equal." << endl;

            coord_row = new int[tempVecInt.size()];
            for(int i=0; i<tempVecInt.size(); i++) coord_row[i]=tempVecInt[i];

            cout << dir_entry.path() << " has been read: " <<  tempVecInt.size() << endl;
            myfile.close();
        }
        else if(dir_entry.path().stem() == "coordinate_col") {
            int tempValInt;
            vector<int> tempVecInt;
            while (myfile >> tempValInt) {
                tempVecInt.push_back(tempValInt);
            }
            if(tempVecInt.size() != nonzerosSize[z]) cout << "Col count not equal." << endl;
            coord_col = new int[tempVecInt.size()];
            for(int i=0; i<tempVecInt.size(); i++) coord_col[i]=tempVecInt[i];

            cout << dir_entry.path() << " has been read: " << tempVecInt.size() << endl;
            myfile.close();
        }
            // else, start reading doubles.
        else if(dir_entry.path().stem() == "coordinate_val"){
            while (myfile >> tempVal) {
                tempVec.push_back(tempVal);
            }
            coord_val = new double[tempVec.size()];
            for(int i=0; i<tempVec.size(); i++) coord_val[i]=tempVec[i];
            if(tempVec.size() != nonzerosSize[z]) cout << "Vals count not equal." << endl;
            myfile.close();
            cout << dir_entry.path() << " has been read: " << tempVec.size() << endl;
        }
        else cout << "unexpected file name: " << dir_entry.path() << endl;
    }
    return 0;
}

int main(int argc, char **argv){
    int n, rowLimit;
    cout << "calling read Coo Format... " << endl;
    //init();
    if(!argv[1]){
        cout << "please provide input matrix index (int): boneS10, Emilia_923, ldoor, af_5_k101, Serena, audikw_1" << endl;
        return -1;
    }
    readCooFormat(atoi(argv[1]));

    n = matrixSize[atoi(argv[1])];
    int inputType = atoi(argv[1]);


    int *coordRow = coord_row;
    int *coordCol = coord_col;
    double *coordval = coord_val;

    int *csrRow = new int[matrixSize[inputType]+1];
    int *csrCol = new int[nonzerosSize[inputType]];
    double *csrval = new double[nonzerosSize[inputType]];

    int nnz = nonzerosSize[inputType];
    int nrow=matrixSize[inputType];

    std::cout  <<  " starts computing coocsr... " << endl;
    coocsr_(&nrow,  &nnz, coordval, coordRow, coordCol, csrval, csrCol, csrRow);
    std::cout  <<  " FINISHED computing coocsr... " << csrval[10] << " " << csrCol[10] << " " << csrRow[10] << endl;

    ofstream myfile1, myfile2, myfile3;
    myfile1.open ("/home/selin/CSR-Data/" + matrix_names[inputType] + "/unbanded/CSRout_row.txt", ios::out | ios::trunc);
    myfile2.open ("/home/selin/CSR-Data/" + matrix_names[inputType] + "/unbanded/CSRout_col.txt", ios::out | ios::trunc);
    myfile3.open ("/home/selin/CSR-Data/" + matrix_names[inputType] + "/unbanded/CSRout_val.txt", ios::out | ios::trunc);

    cout << "Writing to " << "-CSRout_row.txt"  << endl;
    for (int i=0; i<matrixSize[inputType]+1; i++) {
        myfile1 <<  csrRow[i] << '\t';
    }
    cout << "Writing to " << "-CSRout_col.txt"  << endl;
    for (int i=0; i<nonzerosSize[inputType]; i++) {
        myfile2 <<  csrCol[i] << '\t';
    }
    cout << "Writing to " << "-CSRout_val.txt"  << endl;
    for (int i=0; i<nonzerosSize[inputType]; i++){
        myfile3 <<  csrval[i] << '\t';
    }
    myfile1.close();
    myfile2.close();
    myfile3.close();
}
