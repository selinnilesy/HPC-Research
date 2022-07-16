//
// Created by Selin Yıldırım on 7.01.2022.
//

#include <iostream>
#include  "header.h"
#include <limits>
#include <assert.h>
#include <iomanip>
#include <mpi.h>
#include <math.h>
#include <algorithm>

using namespace std;
typedef std::numeric_limits< double > dbl;

vector<int> outer_col, outer_row;
vector<double> outer_val;

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
    }
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

        displs = new int[world_size];
        displs[0] = 0;
        for (int i=1; i<world_size; i++)
            displs[i] = displs[i-1] + pieceSizeArr[i-1];

        matrixRowDiff = new int[n];
        for (int i = 0; i < n; i++) {
            matrixRowDiff[i] = matrixRowptr[i+1] - matrixRowptr[i];
        }
    }

    // scatter pieceSizeArr.
    MPI_Scatter(pieceSizeArr, 1, MPI_INT, &myPieceSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // broadcast expected global piece size.
    MPI_Bcast(&pieceSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

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


    if(my_rank==0) NNZs = new int[world_size];

    myNNZ=0;
    for (int i = 0; i < myPieceSize; i++) {
        myNNZ += myRowDiff[i];
    }
    MPI_Gather(&myNNZ, 1, MPI_INT, NNZs, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if(my_rank==0) {
        for (int i=1; i<world_size; i++)
            displs[i] = displs[i-1] + NNZs[i-1];
    }

    assert (myNNZ!=0);
    myColInd = new int[myNNZ];
    myOffDiags = new double[myNNZ];
    // scatter matrixColind.
    MPI_Scatterv(matrixColind, NNZs, displs, MPI_INT, myColInd, myNNZ, MPI_INT, 0, MPI_COMM_WORLD);
    // scatter matrixOffDiagonal.
    MPI_Scatterv(matrixOffDiagonal, NNZs, displs, MPI_DOUBLE, myOffDiags, myNNZ, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    y = new double[myPieceSize];
    for (int i = 0; i < myPieceSize; i++) y[i] = 0.0;

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

    int neighbourSize;
    double *neighbourX;
    MPI_Status status;
    // send X pieces between neighbours.
    if(my_rank==0){
        MPI_Send(&myPieceSize, 1, MPI_INT, my_rank+1, 0, MPI_COMM_WORLD);
        MPI_Send(myX, myPieceSize, MPI_INT, my_rank+1, 0, MPI_COMM_WORLD);
    }
    else if(my_rank && my_rank < world_size-1){
        MPI_Recv(&neighbourSize, 1, MPI_INT, my_rank-1, 0, MPI_COMM_WORLD, &status);
        neighbourX = new double[neighbourSize];
        MPI_Recv(neighbourX, neighbourSize, MPI_INT, my_rank-1, 0, MPI_COMM_WORLD, &status);
        MPI_Send(&myPieceSize, 1, MPI_INT, my_rank+1, 0, MPI_COMM_WORLD);
        MPI_Send(myX, myPieceSize, MPI_INT, my_rank+1, 0, MPI_COMM_WORLD);
    }
    else if(my_rank==world_size-1){
        MPI_Recv(&neighbourSize, 1, MPI_INT, my_rank-1, 0, MPI_COMM_WORLD, &status);
        neighbourX = new double[neighbourSize];
        MPI_Recv(neighbourX, neighbourSize, MPI_INT, my_rank-1, 0, MPI_COMM_WORLD, &status);
    }
    // send info of X pieces of remaining confs to root.
    int tobesent, temp_tobesent;
    vector<double*> Xsquares_in_process;
    for(int i=0; i<world_size; i++){
        double *ptr = nullptr;
        Xsquares_in_process.push_back(ptr);
    }
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
            for (int i = 0; i < confSquares.size(); i++) {
                if (confSquares[i] == my_rank - 1) continue;
                else {
                    double *temp = new double[pieceSize];
                    MPI_Recv(temp, pieceSize, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
                    cout << "Rank: " << my_rank << " received X[" << confSquares[i] << "]" << endl;
                    Xsquares_in_process[confSquares[i]]=temp;
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
        for (int i = 0; i < send_Xids.size(); i+=2){
            MPI_Send(x+send_Xids[i+1]*pieceSize, pieceSize, MPI_DOUBLE, send_Xids[i], 0, MPI_COMM_WORLD);
            cout << "Rank: " << my_rank << " sent X[" << send_Xids[i+1] << "] to process: " << send_Xids[i] << endl;
        }
    }


    double val;
    int accumIndex=0;
    vector<double*> Ysquares_in_process;
    if(my_rank) {
        for (int i = 0; i < world_size; i++) {
            double *ptr = nullptr;
            Ysquares_in_process.push_back(ptr);
        }
        for (int i = 0; i < confSquares.size(); i++) {
            Ysquares_in_process[confSquares[i]] = new double[pieceSize];
            for (int j = 0; j < pieceSize; j++) Ysquares_in_process[confSquares[i]][j] = 0.0;
        }
    }
    int colIndModulo;
    double start_time = MPI_Wtime();
    for (int i = 0; i < myPieceSize; i++) {
        val = myDiags[i] * myX[i];
        for (int j = accumIndex; j < accumIndex+myRowDiff[i]; j++) {
            colInd = myColInd[j] - 1;
            colIndModulo = fmod(colInd,pieceSize);
            if(colInd < my_rank*pieceSize) {
                if(colInd / pieceSize == my_rank-1) val += myOffDiags[j] * neighbourX[colIndModulo];
                else{
                    val += myOffDiags[j] * (Xsquares_in_process[colInd / pieceSize])[colIndModulo];
                }
                (Ysquares_in_process[colInd / pieceSize])[colIndModulo] -= myOffDiags[j] * myX[i];
            }
            else{
                y[colIndModulo] -= myOffDiags[j] * myX[i];
                val += myOffDiags[j] * myX[colIndModulo];
            }
        }
        y[i] += val;
        accumIndex+=myRowDiff[i];
    }

    MPI_Win window;
    MPI_Win_create(y, sizeof(int), sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &window);
    MPI_Win_fence(0, window);
    // process outer region.
    if(!my_rank) {
        int temp;
        for (int i = 0; i < outer_col.size(); i++) {
            temp = -outer_val[i] * x[outer_row[i] -1];
            MPI_Accumulate(&temp, 1, MPI_DOUBLE, (outer_col[i]-1)/pieceSize, fmod((outer_col[i]-1),pieceSize), 1, MPI_DOUBLE, MPI_SUM, window);
            temp = outer_val[i] * x[outer_col[i] -1];
            MPI_Accumulate(&temp, 1, MPI_DOUBLE, (outer_row[i]-1)/pieceSize, fmod((outer_row[i]-1),pieceSize), 1, MPI_DOUBLE, MPI_SUM, window);
        }
    }
    MPI_Win_fence(0, window);
    // accumulate y results.
    MPI_Win_fence(0, window);
    if(my_rank) {
        for (int i = 0; i < confSquares.size(); i++) {
            MPI_Accumulate(Ysquares_in_process[confSquares[i]], pieceSize, MPI_DOUBLE, confSquares[i], 0, pieceSize,
                           MPI_DOUBLE, MPI_SUM, window);
        }
    }
    MPI_Win_fence(0, window);
    double end_time = MPI_Wtime();
    if(!my_rank) printf("It took me %f seconds for parallel run.\n", end_time-start_time);
    // Destroy the window
    MPI_Win_free(&window);


    double *output;
    if(my_rank==0) {
        delete [] x;
        output = new double[n];
        displs[0] = 0;
        for (int i=1; i<world_size; i++)
            displs[i] = displs[i-1] + pieceSizeArr[i-1];
    }

    MPI_Gatherv(y, myPieceSize, MPI_DOUBLE, output, pieceSizeArr, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if(my_rank==0){
        ofstream myfile;
        myfile.open ("/home/selin/3way-Par-Results/" + matrix_names[inputType] + "/result.txt", ios::out | ios::trunc);
        cout << "Writing to output... " << endl;
        for (int i=0; i<n; i++) {
            myfile << std::fixed << std::setprecision(dbl::max_digits10) << output[i] << '\t';
        }
        myfile.close();
        cout << "Completed output... " << endl;
    }


    // termination
    if(my_rank==0) {
        //submat_(1 int* i1, int* i2, int* j1,int* j2, double* a,int* ja,int* ia,int* nr,int* nc,double* ao,int* jao,int* iao);
        //delete[] x;
        delete[] matrixRowptr;
        delete[] matrixRowDiff;
        delete[] matrixColind;
        delete[] matrixOffDiagonal;
        delete[] matrixDiagonal;
        delete[] NNZs;
        delete[] pieceSizeArr;
    }


    // Finalize MPI
    // This must always be called after all other MPI functions
    delete [] myX;
    delete [] myColInd;
    delete [] myOffDiags;
    delete [] myDiags;
    delete [] myRowDiff;

    cout << "Rank: " << my_rank << " Bitti." << endl << flush;

    MPI_Finalize();

    return 0;
}