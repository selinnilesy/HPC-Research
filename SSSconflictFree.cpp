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

int readSSSFormat(int z) {
    double tempVal;
    vector<double> tempVec;
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

    int elmCountPerRow, colInd, rowBegin, val;
    int innerBandwith, middleBandwith ,nnz_n_Ratio;
    innerBandwith = (nnz_n_Ratios[inputType]*bandwithProportions[inputType]/4);
    int maxJ=-1;

    static_assert(rowptrSize[0] == matrixSize[inputType]+1);

    // two new SSS storage for inner and outer regions
    vector<int> colind_inner, rowptr_inner;
    vector<int> colind_outer, rowptr_outer;
    vector<int> colind_middle, rowptr_middle;
    vector<double> vals_inner, vals_outer, vals_middle;
    int counter_inner, counter_middle, counter_outer;
    // initialize row ptr vectors
    for (int i = 0; i < rowptrSize[0] ; i++) {
        rowptr_inner.push_back(0);
        rowptr_middle.push_back(0);
        rowptr_outer.push_back(0);
    }

    for (int i = 0; i < rowptrSize[0] - 1; i++) {
        // row ptrs start from 1 !!!
        rowBegin = matrixRowptr[i] - 1;
        elmCountPerRow = matrixRowptr[i + 1] - (matrixRowptr[i] ;
        maxJ=1;
        counter_inner=counter_middle=counter_outer = 0;
        for (int j = 0; j < elmCountPerRow; j++) {
            // i = row indexi
            // colind = column indexi
            // val = degeri
            colInd = matrixColind[rowBegin + j];
            val = matrixOffDiagonal[rowBegin + j];

            if(colInd > maxJ) maxJ=colInd;
            else cout << "ON THE SAME ROW, COL INDEX HAS BEEN SMALLED. NOT IN ASCENDING ORDER: maxj, colInd " << maxJ <<" " << colInd << endl;
            // inner Dense Region
            if(colInd >= i - innerBandwith){
                for (int k = i+1; k < rowptrSize[0] ; k++) {
                    rowptr_inner[i+1]++;
                }
                colind_inner[rowptr_inner[i] + counter_inner] = colInd;
                vals_inner[rowptr_inner[i] + counter_inner] = val;
                counter_inner++;
            }
            // middle Region
            else if(j > innerBandwith){
                for (int k = i+1; k < rowptrSize[0] ; k++) {
                    rowptr_middle[i+1]++;
                }
                colind_middle[rowptr_middle[i] + counter_middle] = colInd;
                vals_middle[rowptr_middle[i] + counter_middle] = val;
                counter_middle++;
            }
            // outer Dense Region
            else{
                for (int k = i+1; k < rowptrSize[0] ; k++) {
                    rowptr_outer[i+1]++;
                }
                colind_outer[rowptr_outer[i] + counter_outer] = colInd;
                vals_outer[rowptr_outer[i] + counter_outer] = val;
                counter_outer++;
            }
        }
    }


    delete [] matrixOffDiagonal;
    delete [] matrixColind;
    delete [] matrixRowptr;
}