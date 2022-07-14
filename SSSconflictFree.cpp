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
#include <algorithm>
//#include "rcmtest.cpp"
//#include "geeks.cpp"

using namespace std;
#define MATRIX_COUNT 6

vector<int> colind_up, rowptr_up;
vector<double> offdiag_up;

vector<int> outer_col, outer_row;
vector<double> outer_val;

typedef std::numeric_limits< double > dbl;

extern "C" {
    extern void submat_(int* job, int*i1, int* i2, int* j1,int* j2, double* a,int* ja,int* ia,int* nr,int* nc,double* ao,int* jao,int* iao);
}

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

    int n, inputType, pieceSize, nnz;
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
        //cout << my_rank << " : nnz " << nnz << endl;
    }

        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&inputType, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if(my_rank==0) {
        matrixOffDiagonal = valuesPtrs[0];
        matrixDiagonal = dvaluesPtrs[0];
        matrixColind = colindPtrs[0];
        matrixRowptr = rowptrPtrs[0];
        x = new double[n];
        for (int i = 0; i < n; i++) x[i] = 1.0;

        pieceSizeArr = new int[world_size];
        for (int i = 0; i < world_size - 1; i++) pieceSizeArr[i] = pieceSize;
        pieceSizeArr[world_size - 1] = n - (world_size - 1) * pieceSize;
        displs = pieceSizeArr;

        matrixRowDiff = new int[n];
        for (int i = 0; i < n; i++) {
            matrixRowDiff[i] = matrixRowptr[i+1] - matrixRowptr[i];
        }
    }

    // scatter pieceSizeArr.
    MPI_Scatter(pieceSizeArr, 1, MPI_INT, &myPieceSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    //cout << my_rank << " : myPieceSize " << myPieceSize << endl;
    // broadcast expected global piece size.
    MPI_Bcast(&pieceSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    //cout << my_rank << " : pieceSize " << pieceSize << endl;

    assert (myPieceSize!=0);
    myX = new double[myPieceSize];
    myRowDiff = new int[myPieceSize];
    myDiags = new double[myPieceSize];
    // scatter x.
    MPI_Scatterv(x, pieceSizeArr, displs, MPI_DOUBLE, myX, myPieceSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //cout << my_rank << " : X " << myX[10] <<  " " << myX[50] << endl;
    // scatter matrixRowptr.
    MPI_Scatterv(matrixRowDiff, pieceSizeArr, displs, MPI_INT, myRowDiff, myPieceSize, MPI_INT, 0, MPI_COMM_WORLD);
    // scatter matrixDiagonal.
    MPI_Scatterv(matrixDiagonal, pieceSizeArr, displs, MPI_DOUBLE, myDiags, myPieceSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    if(my_rank==0) NNZs = new int[world_size];

    myNNZ=0;
    for (int i = 0; i < myPieceSize; i++) {
        myNNZ += myRowDiff[i];
    }
    MPI_Gather(&myNNZ, 1, MPI_INT, NNZs, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if(my_rank==0)  displs = NNZs;

    //cout << my_rank << " : myNNZ " << myNNZ << endl;
    assert (myNNZ!=0);
    myColInd = new int[myNNZ];
    myOffDiags = new double[myNNZ];
    // scatter matrixColind.
    MPI_Scatterv(matrixColind, NNZs, displs, MPI_INT, myColInd, myNNZ, MPI_INT, 0, MPI_COMM_WORLD);
    // scatter matrixOffDiagonal.
    MPI_Scatterv(matrixOffDiagonal, NNZs, displs, MPI_DOUBLE, myOffDiags, myNNZ, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //cout << my_rank << " : myColInd " << myColInd[10] << " " << myColInd[50] << endl;

    y = new double[myPieceSize];
    for (int i = 0; i < myPieceSize; i++) y[i] = 0.0;

    //int *confSquares = new int[world_size];
    vector<int> confSquares;
    int accum_colInd=0, colInd, lastConfSquare=-1;
    for (int i = 0; i < myPieceSize; i++) {
        for (int j =accum_colInd ; j< accum_colInd + myRowDiff[i]; j++) {
            colInd = myColInd[j] - 1;
            if(colInd < my_rank*pieceSize) {
                if(lastConfSquare!=colInd/pieceSize && find(confSquares.begin(), confSquares.end(), colInd/pieceSize) == confSquares.end() ) {
                    lastConfSquare = colInd/pieceSize;
                    confSquares.push_back(lastConfSquare);
                }
            }
        }
        accum_colInd+=myRowDiff[i];
    }
    for (int i = 0; i < confSquares.size(); i++) cout << "Rank: " << my_rank << " confs in squares: " << confSquares[i] << endl;

    int neighbourSize, *neighbourX;
    MPI_Status status;
    if(my_rank==0){
        MPI_Send(&myPieceSize, 1, MPI_INT, my_rank+1, 0, MPI_COMM_WORLD);
        MPI_Send(myX, myPieceSize, MPI_INT, my_rank+1, 0, MPI_COMM_WORLD);
        cout << "Rank: " << my_rank << " sent X[0]" << endl;
    }
    else if(my_rank && my_rank < world_size-1){
        MPI_Recv(&neighbourSize, 1, MPI_INT, my_rank-1, 0, MPI_COMM_WORLD, &status);
        neighbourX = new int[neighbourSize];
        MPI_Recv(neighbourX, neighbourSize, MPI_INT, my_rank-1, 0, MPI_COMM_WORLD, &status);
        cout << "Rank: " << my_rank << " received X[" << my_rank-1 << "]" << endl;
        MPI_Send(&myPieceSize, 1, MPI_INT, my_rank+1, 0, MPI_COMM_WORLD);
        MPI_Send(myX, myPieceSize, MPI_INT, my_rank+1, 0, MPI_COMM_WORLD);
        cout << "Rank: " << my_rank << " sent X[" << my_rank << "]" << endl;
    }
    else if(my_rank==world_size-1){
        MPI_Recv(&neighbourSize, 1, MPI_INT, my_rank-1, 0, MPI_COMM_WORLD, &status);
        neighbourX = new int[neighbourSize];
        MPI_Recv(neighbourX, neighbourSize, MPI_INT, my_rank-1, 0, MPI_COMM_WORLD, &status);
        cout << "Rank: " << my_rank << " received X[" << my_rank-1 << "]" << endl;
    }
    // send info of remaining confs to root.
    int tobesent, temp_tobesent;
    // store first process id then needed square id.
    vector<int> send_Xids;
    if(my_rank){
        tobesent = confSquares.size() -1 ;
        MPI_Send(&tobesent, 1, MPI_INT, 0, my_rank, MPI_COMM_WORLD);
        if(tobesent) {
            for (int i = 0; i < confSquares.size(); i++) {
                if (confSquares[i] == my_rank - 1) continue;
                else {
                    MPI_Send(&(confSquares[i]), 1, MPI_INT, 0, my_rank, MPI_COMM_WORLD);
                }
            }
        }
    }
    else{
        for(int i=1; i<world_size; i++){
            MPI_Recv(&tobesent, 1, MPI_INT, i, i, MPI_COMM_WORLD, &status);
            while(tobesent--){
                send_Xids.push_back(i);
                MPI_Recv(&temp_tobesent, 1, MPI_INT, i, i, MPI_COMM_WORLD, &status);
                send_Xids.push_back(temp_tobesent);
            }
        }
        for (int i = 0; i < send_Xids.size(); i++) cout << "Rank: " << my_rank << " confs in root to be sent: " << send_Xids[i] << endl;
    }

    if(my_rank==0) {
        //submat_(1 int* i1, int* i2, int* j1,int* j2, double* a,int* ja,int* ia,int* nr,int* nc,double* ao,int* jao,int* iao);
        delete [] x;
        delete [] matrixRowptr;
        delete [] matrixRowDiff;
        delete [] matrixColind;
        delete [] matrixOffDiagonal;
        delete [] matrixDiagonal;
        delete [] NNZs;
        delete [] pieceSizeArr;

        cout << "Rank: " << my_rank << " Bitti." << endl;
    }
    else{
        cout <<"Rank: " << my_rank << " Bitti." << endl;
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