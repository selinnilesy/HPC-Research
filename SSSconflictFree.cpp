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

int readSSSFormat(int z) {
    double tempVal;
    vector<double> tempVec;
        const fs::path matrixFolder{"/home/selin/CSR-Data/" + matrix_names[z]};
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
            else if(dir_entry.path().stem() == "coordinate-row") {
                int tempValInt;
                vector<int> tempVecInt;
                while (myfile >> tempValInt) {
                    tempVecInt.push_back(tempValInt);
                }
                coord_row = new int[tempVecInt.size()];
                for(int i=0; i<tempVecInt.size(); i++) coord_row[i]=tempVecInt[i];

                cout << dir_entry.path() << " has been read." << endl;
                myfile.close();
                continue;
            }
            else if(dir_entry.path().stem() == "col") {
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
            else if(dir_entry.path().stem() == "coordinate-col") {
                int tempValInt;
                vector<int> tempVecInt;
                while (myfile >> tempValInt) {
                    tempVecInt.push_back(tempValInt);
                }
                coord_col = new int[tempVecInt.size()];
                for(int i=0; i<tempVecInt.size(); i++) coord_col[i]=tempVecInt[i];

                cout << dir_entry.path() << " has been read." << endl;
                myfile.close();
                continue;
            }
            // else, start reading doubles.
            while (myfile >> tempVal) {
                tempVec.push_back(tempVal);
            }

            if(dir_entry.path().stem() == "dvalues"){
                dvaluesPtrs.push_back(new double[tempVec.size()]);
                double *temp = dvaluesPtrs[0];
                for(int i=0; i<tempVec.size(); i++) temp[i]=tempVec[i];
                dvaluesSize.push_back(tempVec.size());
            }
            else if(dir_entry.path().stem() == "values"){
                valuesPtrs.push_back(new double[tempVec.size()]);
                double *temp = valuesPtrs[0];
                for(int i=0; i<tempVec.size(); i++) temp[i]=tempVec[i];
                valuesSize.push_back(tempVec.size());
            }
            else if(dir_entry.path().stem() == "coordinate-val"){
                coord_val = new double[tempVec.size()];
                for(int i=0; i<tempVec.size(); i++) coord_val[i]=tempVec[i];
            }
            else cout << "unexpected file name: " << dir_entry.path() << endl;
            cout << dir_entry.path() << " has been read." << endl;

            tempVec.clear();
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
    readSSSFormat(atoi(argv[1]));

    n = matrixSize[atoi(argv[1])];
    int inputType = atoi(argv[1]);

    double *matrixOffDiagonal = valuesPtrs[0];
    int *matrixColind = colindPtrs[0];
    int *matrixRowptr= rowptrPtrs[0];

    delete [] matrixOffDiagonal;
    delete [] matrixColind;
    delete [] matrixRowptr;

    int *banded_coordRow = coord_row;
    int *banded_coordCol = coord_col;
    double *banded_coordval = coord_val;

    std::cout  <<  " starts computing coocsr... " << endl;
    int *banded_csrRow = new int[matrixSize[inputType]+1];
    int *banded_csrCol = new int[nonzerosSize[inputType]];
    double *banded_csrval = new double[nonzerosSize[inputType]];
    std::cout  <<  "iao address: " << endl;
    std::cout  << banded_csrRow << endl;
    std::cout  <<  "wrtie on iao : " << endl;
    banded_csrRow[0] = 0;
    std::cout  <<  "read iao: " << endl;
    std::cout  << banded_csrRow[0] << endl;
    int nnz = nonzerosSize[inputType];
    int nrow=matrixSize[inputType];
    coocsr_(&nrow,  &nnz, banded_coordval, banded_coordRow, banded_coordCol, banded_csrval, banded_csrCol, banded_csrRow);
    std::cout  <<  " FINISHED computing coocsr... " << endl;
}
