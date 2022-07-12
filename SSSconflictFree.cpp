//
// Created by Selin Yıldırım on 7.01.2022.
//

#include <iostream>
//#include <omp.h>
#include  "header.h"
#include <limits>
#include <assert.h>
#include <iomanip>
#include <mpi.h>
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

int main(int argc, char **argv) {
    int my_rank,  world_size;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    int n, inputType, pieceSize, local_pieceSize, nnz;
    cout <<"World size: "  << world_size  << endl;
    double *matrixOffDiagonal, *matrixDiagonal;
    double *x, *y;
    int *matrixColind, *matrixRowDiff, *matrixRowptr;
    int *pieceSizeArr, *NNZs;
    int myPieceSize, myNNZ;
    int *myRowDiff, *myColInd;
    double *myX, *myOffDiags, *myDiags;
    int *displs;
    if(my_rank==0) {
        if (!argv[1]) {
            cout << "please provide input matrix index (int): boneS10, Emilia_923, ldoor, af_5_k101, Serena, audikw_1"
                 << endl;
            return -1;
        }
        readSSSFormat(atoi(argv[1]));
        readCSRFormat(atoi(argv[1]));

        n = matrixSize[atoi(argv[1])];
        inputType = atoi(argv[1]);
        pieceSize = n / world_size;
        nnz = colindSize[inputType];

        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&inputType, 1, MPI_INT, 0, MPI_COMM_WORLD);

        matrixOffDiagonal = valuesPtrs[0];
        matrixDiagonal = dvaluesPtrs[0];
        matrixColind = colindPtrs[0];
        matrixRowptr = rowptrPtrs[0];
        x = new double[n];
        for (int i = 0; i < n; i++) x[i] = 1.0;

        pieceSizeArr = new int[world_size];
        for (int i = 0; i < world_size - 1; i++) pieceSizeArr[i] = pieceSize;
        pieceSizeArr[world_size - 1] = n - (world_size - 1) * pieceSize;

        matrixRowDiff = new int[n];
        for (int i = 0; i < n; i++) matrixRowDiff[i] = matrixRowptr[i+1] - matrixRowptr[i];

        displs = pieceSizeArr;
    }

    // scatter pieceSizeArr.
    MPI_Scatter(pieceSizeArr, 1, MPI_INT, &myPieceSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    cout << my_rank << " : myPieceSize " << myPieceSize << endl;
    // broadcast expected global piece size.
    if(my_rank==0) {
        MPI_Bcast(&pieceSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
        cout << my_rank << " : local_pieceSize " << local_pieceSize << endl;
    }
    else{
        MPI_Bcast(&local_pieceSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
        cout << my_rank << " : local_pieceSize " << local_pieceSize << endl;
    }

    assert (myPieceSize!=0);
    myX = new double[myPieceSize];
    myRowDiff = new int[myPieceSize];
    myDiags = new double[myPieceSize];
    // scatter x.
    MPI_Scatterv(x, pieceSizeArr, displs, MPI_DOUBLE, myX, myPieceSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // scatter matrixRowptr.
    MPI_Scatterv(matrixRowDiff, pieceSizeArr, displs, MPI_INT, myRowDiff, myPieceSize, MPI_INT, 0, MPI_COMM_WORLD);
    // scatter matrixDiagonal.
    MPI_Scatterv(matrixDiagonal, pieceSizeArr, displs, MPI_DOUBLE, myDiags, myPieceSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if(my_rank==0) {
        NNZs = new int[world_size];
        for (int i = 0; i < world_size; i++) {
            if(i==world_size-1) NNZs[i] = matrixRowptr[n] - matrixRowptr[(world_size-1) * local_pieceSize];
            else NNZs[i] = matrixRowptr[(i+1) * local_pieceSize];
        }
        displs = NNZs;
    }

    MPI_Scatter(NNZs, 1, MPI_INT, &myNNZ, 1, MPI_INT, 0, MPI_COMM_WORLD);
    assert (myNNZ!=0);
    myColInd = new int[myNNZ];
    myOffDiags = new double[myNNZ];
    // scatter matrixColind.
    MPI_Scatterv(matrixColind, NNZs, displs, MPI_INT, myColInd, myNNZ, MPI_INT, 0, MPI_COMM_WORLD);
    // scatter matrixOffDiagonal.
    MPI_Scatterv(matrixOffDiagonal, NNZs, displs, MPI_DOUBLE, myOffDiags, myNNZ, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    y = new double[myPieceSize];
    for (int i = 0; i < myPieceSize; i++) y[i] = 0.0;

    cout << "Conflict: " << my_rank << " my nnzs : " << myNNZ << " " << myColInd[10] << myColInd[50] << endl;
    int accum_colInd=0, colInd;
    for (int i = 0; i < myPieceSize; i++) {
        for (int j =accum_colInd ; j< accum_colInd + myRowDiff[i]; j++) {
            colInd = myColInd[j] - 1;
            //if(colInd < my_rank*pieceSize) cout << "Conflict: " << my_rank << " in square: " << colInd/pieceSize << endl;
        }
        accum_colInd+=myRowDiff[i];
    }

    if(my_rank==0) {
        delete [] x;
        delete [] matrixRowptr;
        delete [] matrixRowDiff;
        delete [] matrixColind;
        delete [] matrixOffDiagonal;
        delete [] matrixDiagonal;
        delete [] NNZs;
        delete [] pieceSizeArr;
        cout << my_rank << ": Bitti." << endl;
    }
    else{
        cout << my_rank << ": Bitti." << endl;
    }

    // Finalize MPI
    // This must always be called after all other MPI functions
    delete [] myX;
    delete [] myColInd;
    delete [] myOffDiags;
    delete [] myDiags;
    delete [] myRowDiff;
    MPI_Finalize();

    return 0;
}