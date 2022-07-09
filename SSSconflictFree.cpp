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
int *csr_row, *csr_col;
double *csr_val;

int offd_count_innerbandw;


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
//extern void coocsr_(int *nrow, int *nnz, double* a, int* ir,int* jc, double* ao, int* jao, int* iao);
extern void csrsss_(int *nrow, int *nnz, double* a, int *ja, int *ia, bool* sorted, double *diag, double* al, int *jal, int *ial, double* au);
}

int readCSRFormat(int z, double ratio) {
    double tempVal;
    vector<double> tempVec;
    const fs::path matrixFolder{"/home/selin/Split-Data/" + matrix_names[z] +"/middle/CSR-Data"};
    for(auto const& dir_entry: fs::directory_iterator{matrixFolder}){
        std::fstream myfile(dir_entry.path(), std::ios_base::in);
        if(dir_entry.path().stem() == "1.000000-" + to_string(ratio) + "-row") {
            int tempValInt;
            vector<int> tempVecInt;
            while (myfile >> tempValInt) {
                tempVecInt.push_back(tempValInt);
            }
            csr_row = new int[tempVecInt.size()];
            for(int i=0; i<tempVecInt.size(); i++) csr_row[i]=tempVecInt[i];

            cout << dir_entry.path() << " has been read." << endl;
        }
        else if(dir_entry.path().stem() == "1.000000-" + to_string(ratio) +"-col") {
            int tempValInt;
            vector<int> tempVecInt;
            while (myfile >> tempValInt) {
                tempVecInt.push_back(tempValInt);
            }
            csr_col = new int[tempVecInt.size()];
            for(int i=0; i<tempVecInt.size(); i++) csr_col[i]=tempVecInt[i];
            offd_count_innerbandw = tempVecInt.size();
            cout << dir_entry.path() << " has been read." << endl;
        }
        else if(dir_entry.path().stem() == "1.000000-" + to_string(ratio) +"-val") {
            double tempVall;
            vector<double> tempVec;
            while (myfile >> tempVall) {
                tempVec.push_back(tempVall);
            }
            csr_val = new double[tempVec.size()];
            for(int i=0; i<tempVec.size(); i++) csr_val[i]=tempVec[i];

            cout << dir_entry.path() << " has been read." << endl;
        }
        myfile.close();
    }
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
    int inputType = atoi(argv[1]);
    readCSRFormat(inputType, bandwithSize[inputType]-3);
    n = matrixSize[atoi(argv[1])];

    std::cout  <<  " starts computing csrsss... " << endl;

    int nnz = offd_count_innerbandw;
    int nrow=matrixSize[inputType];

    double *diag = new double[nrow];
    int *rowptr = new int[nrow+1];
    // strict lower part
    int *colinds_lower = new int[nnz];
    double *vals_lower = new double[nnz];
    double *vals_upper = new double[nnz];


    bool sorted = 0;
    /*
    vector<double> test;
    cout << "test diag" << endl;
    int z_count=0;
    for(int i=0; i<923136+1; i++){
        int rowptr = csr_row[i]-1;
        int n_rowptr = csr_row[i+1]-1;
        for(int x=rowptr; x<n_rowptr; x++){
            int col = csr_col[rowptr+x] - 1;
            if(col==i){
                if(csr_val[rowptr+x]==0) z_count++;
                test.push_back(csr_val[rowptr+x]);
            }
        }
    }
    cout << "found: " << z_count ;
    ofstream test_f;
    test_f.open ("/home/selin/Split-Data/" + matrix_names[inputType] + "/inner/CSR-Data/test_diag.txt", ios::out | ios::trunc);
    for (int i=0; i<test.size(); i++) {
        test_f << test[i] << '\t';
    }
    test_f.close();
     */


    cout << "starts computing csrsss_..." << endl;
    csrsss_(&nrow,&nnz, csr_val, csr_col, csr_row, &sorted, diag, vals_lower, colinds_lower, rowptr, vals_upper);
    std::cout  <<  " finished computing csrsss... " << diag[10] << " " << vals_lower[10] << " " << colinds_lower[10]<< " " << rowptr[10] << endl;

    ofstream myfile1, myfile2, myfile3,myfile4;
    myfile1.open ("/home/selin/Split-Data/" + matrix_names[inputType] + "/middle/CSR-Data/rowptr.txt", ios::out | ios::trunc);
    myfile2.open ("/home/selin/Split-Data/" + matrix_names[inputType] + "/middle/CSR-Data/colind.txt", ios::out | ios::trunc);
    myfile3.open ("/home/selin/Split-Data/" + matrix_names[inputType] + "/middle/CSR-Data/vals.txt", ios::out | ios::trunc);
    myfile4.open ("/home/selin/Split-Data/" + matrix_names[inputType] + "/middle/CSR-Data/diag.txt", ios::out | ios::trunc);

    cout << "Writing to " << "-SSSout_rowptr.txt"  << endl;
    for (int i=0; i<nrow+1; i++) {
        myfile1 << rowptr[i] << '\t';
    }
    cout << "Writing to " << "-SSSout_colind.txt"  << endl;
    for (int i=0; i<nnz; i++) {
        myfile2 << colinds_lower[i] << '\t';
    }
    cout << "Writing to " << "-SSSout_vals.txt"  << endl;
    for (int i=0; i<nnz; i++){
        myfile3 << vals_lower[i] << '\t';
    }
    /*
    cout << "Writing to " << "-SSSout_diag.txt"  << endl;
    for (int i=0; i<nrow; i++){
        myfile4 << diag[i] << '\t';
    }
     */
    myfile1.close();
    myfile2.close();
    myfile3.close();
    myfile4.close();

}