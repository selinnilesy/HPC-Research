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
int readCSRFormat(int z, int cikartilanBW) {
    fs::path matrixFolder;
    matrixFolder = "/home/selin/Split-Data/" + matrix_names[z] + "/middle/CSR-Data";
    for(auto const& dir_entry: fs::directory_iterator{matrixFolder}){
        std::fstream myfile(dir_entry.path(), std::ios_base::in);
        if(dir_entry.path().stem() == ("1.000000-" + to_string((double) bandwithSize[z]-cikartilanBW) + "-row")) {
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
        else if(dir_entry.path().stem() == ("1.000000-" + to_string( (double) bandwithSize[z]-cikartilanBW) + "-col")) {
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
        else if(dir_entry.path().stem() == ("1.000000-" + to_string((double) bandwithSize[z]-cikartilanBW) + "-val")){
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
        if (dir_entry.path().stem() == ("1.000000-" + to_string((double) bandwithSize[z] - cikartilanBW) + "-row")) {
            int tempValInt;
            while (myfile >> tempValInt) {
                outer_row.push_back(tempValInt);
            }
            cout << dir_entry.path() << " has been read with size: " << outer_row.size() << endl;
            myfile.close();
        } else if (dir_entry.path().stem() == ("1.000000-" + to_string((double) bandwithSize[z] - cikartilanBW) + "-col")) {
            int tempValInt;
            while (myfile >> tempValInt) {
                outer_col.push_back(tempValInt);
            }
            cout << dir_entry.path() << " has been read with size: " << outer_col.size() << endl;
            myfile.close();
        } else if (dir_entry.path().stem() == ("1.000000-" + to_string((double) bandwithSize[z] - cikartilanBW) + "val")) {
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

    // Obtain the group of processes in the world communicator
    MPI_Group world_group, new_group;
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);
    MPI_Comm newworld;

    int n, inputType, pieceSize, nnz;
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
        readCSRFormat(atoi(argv[1]), 20);

        n = matrixSize[atoi(argv[1])];
        inputType = atoi(argv[1]);
        pieceSize = floor(n / world_size);
        nnz = colindSize[inputType];
        //cout << my_rank << " : nnz " << nnz << endl;
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&inputType, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int remainingPieceSize;

    if(my_rank==0) {
        matrixOffDiagonal = valuesPtrs[0];
        matrixDiagonal = dvaluesPtrs[0];
        matrixColind = colindPtrs[0];
        matrixRowptr = rowptrPtrs[0];
        x = new double[n];
        for (int i = 0; i < n; i++) x[i] = 1.0;
    }
    pieceSizeArr = new int[world_size];
    if(my_rank==0)  {
        int balanceFactor = 1;
        pieceSizeArr[0] = balanceFactor * pieceSize + ((n - balanceFactor * pieceSize) % (world_size - 1));
        remainingPieceSize = (n - pieceSizeArr[0]) / (world_size - 1);
        cout << "root okk " << remainingPieceSize << " 0th: " << pieceSizeArr[0] << endl;

        for (int i = 1; i < world_size; i++) pieceSizeArr[i] = remainingPieceSize;
        for (int i = 0; i < world_size; i++) {
            cout << pieceSizeArr[i] << '\t' << flush;
        }
    }
    displs = new int[world_size];
    if(my_rank==0) {
        displs[0] = 0;
        for (int i=1; i<world_size; i++)
            displs[i] = displs[i-1] + pieceSizeArr[i-1];

        matrixRowDiff = new int[n];
        for (int i = 0; i < n; i++) {
            matrixRowDiff[i] = matrixRowptr[i+1] - matrixRowptr[i];
        }
    }

    MPI_Bcast( (void*) displs, world_size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast( (void*) pieceSizeArr, world_size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&remainingPieceSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

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
   cout << "my nnz: " << myNNZ << '\t' << flush;

    MPI_Gather(&myNNZ, 1, MPI_INT, NNZs, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int firstEmptyProcess = world_size;
    if(my_rank==0) {
        for (int i=1; i<world_size; i++) {
            if( NNZs[i] == 0) {
                if(firstEmptyProcess > i) firstEmptyProcess = i;
            }
        }
    }
    MPI_Bcast(&firstEmptyProcess, 1, MPI_INT, 0, MPI_COMM_WORLD);
    y = new double[myPieceSize];
    for (int i = 0; i < myPieceSize; i++) y[i] = 0.0;

    int neighbourSize;
    double* closingYs;
    double *output;
    if(!my_rank){
        output = new double[n];
        for (int i = 0; i < n; i++) output[i]=0.0;
    }
    MPI_Request request;
    MPI_Status status;
    int totalSize=0;
    double time=0.0, end_time, start_time;
    if(firstEmptyProcess < world_size) {
        MPI_Win window2;
        MPI_Win_create(output, myPieceSize*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &window2);
        MPI_Win_fence(0, window2);
        if(my_rank >= firstEmptyProcess){
            // process diagonals before quitting
            start_time = MPI_Wtime();
            for (int i = 0; i < myPieceSize; i++) y[i] = myDiags[i] * myX[i];
            // accumulate diag results onto output vector.
            MPI_Accumulate(y, myPieceSize, MPI_DOUBLE, 0, displs[my_rank],
                           myPieceSize, MPI_DOUBLE, MPI_SUM, window2);
            end_time = MPI_Wtime();
            time+= end_time-start_time;
            cout << my_rank<< ": Communication - empty process " << time << endl;
        }
        MPI_Win_fence(0, window2);
        MPI_Win_free(&window2);
        MPI_Barrier(MPI_COMM_WORLD);
        int ranges[1][3] = {{firstEmptyProcess, world_size - 1, 1}};
        MPI_Group_range_excl(world_group, 1, ranges, &new_group);
        // Create a new communicator
        MPI_Comm_create(MPI_COMM_WORLD, new_group, &newworld);
        if (newworld == MPI_COMM_NULL) {
            delete [] y;
            MPI_Finalize();
            exit(0);
        }
    }
    else{
        MPI_Comm_create(MPI_COMM_WORLD, world_group, &newworld);
    }
    MPI_Comm_rank(newworld, &my_rank);
    MPI_Comm_size(newworld, &world_size);

    myColInd = new int[myNNZ];
    myOffDiags = new double[myNNZ];
    int *NNZs2;

    // scatter matrixColind.
    if(my_rank==0) {
        NNZs2 = new int[firstEmptyProcess];
        for (int i = 0; i < firstEmptyProcess; i++) {
            NNZs2[i] = NNZs[i];
        }
    }
    delete [] displs;
    displs= new int[firstEmptyProcess+1];

    if(my_rank==0){
        displs[0] = 0;
        for (int i=1; i<firstEmptyProcess; i++) {
            displs[i] = displs[i - 1] + NNZs[i - 1];
        }
    }
    MPI_Bcast( (void*) displs, firstEmptyProcess, MPI_INT, 0, newworld);

    MPI_Scatterv(matrixColind, NNZs2, displs, MPI_INT, myColInd, myNNZ, MPI_INT, 0, newworld);
    // scatter matrixOffDiagonal.
    MPI_Scatterv(matrixOffDiagonal, NNZs2, displs, MPI_DOUBLE, myOffDiags, myNNZ, MPI_DOUBLE, 0, newworld);


    displs[0] = 0;
    for (int i=1; i<firstEmptyProcess+1; i++)
        displs[i] = displs[i-1] + pieceSizeArr[i-1];


    vector<int> confSquares;
    int squareIndex;
    int accum_colInd=0, colInd, lastConfSquare=-1;
    for (int i = 0; i < myPieceSize; i++) {
        for (int j =accum_colInd ; j< accum_colInd + myRowDiff[i]; j++) {
            colInd = myColInd[j] - 1;
            if(colInd < displs[1]) {
                squareIndex = 0;
            }
            else{
                squareIndex = 1 + (colInd - displs[1] )/remainingPieceSize;
            }
            if(lastConfSquare!=squareIndex && find(confSquares.begin(), confSquares.end(), squareIndex) == confSquares.end() ) {
                lastConfSquare = squareIndex;
                confSquares.push_back(squareIndex);
            }
        }
        accum_colInd+=myRowDiff[i];
    }
    for (int i = 0; i < confSquares.size(); i++)  cout << "Rank: " << my_rank << " confs in squares: " << confSquares[i] << '\t' << flush;

    double *neighbourX;
    // send X pieces between neighbours.
    if(my_rank==0){
        MPI_Send(&myPieceSize, 1, MPI_INT, my_rank+1, 0, newworld);
        MPI_Send(myX, myPieceSize, MPI_DOUBLE, my_rank+1, 0, newworld);
    }
    else if(my_rank && my_rank < firstEmptyProcess-1){
        MPI_Recv(&neighbourSize, 1, MPI_INT, my_rank-1, 0, newworld, &status);
        neighbourX = new double[neighbourSize];
        MPI_Recv(neighbourX, neighbourSize, MPI_DOUBLE, my_rank-1, 0, newworld, &status);
        MPI_Send(&myPieceSize, 1, MPI_INT, my_rank+1, 0, newworld);
        MPI_Send(myX, myPieceSize, MPI_DOUBLE, my_rank+1, 0, newworld);
    }
    else if(my_rank==firstEmptyProcess-1){
        MPI_Recv(&neighbourSize, 1, MPI_INT, my_rank-1, 0, newworld, &status);
        neighbourX = new double[neighbourSize];
        MPI_Recv(neighbourX, neighbourSize, MPI_DOUBLE, my_rank-1, 0, newworld, &status);
    }

    if(my_rank==firstEmptyProcess-1) cout << my_rank << ": done with sending neigbours " <<  endl << flush;
    // send info of X pieces of remaining confs to root.
    int tobesent, temp_tobesent;
    vector<double*> Xsquares_in_process;
    for(int i=0; i<firstEmptyProcess; i++){
        double *ptr = nullptr;
        Xsquares_in_process.push_back(ptr);
    }
    // store first process id then needed square id.
    vector<int> send_Xids;
    if(my_rank){
        tobesent = confSquares.size() -1 ;
        MPI_Send(&tobesent, 1, MPI_INT, 0, my_rank, newworld);
        if(tobesent) {
            for (int i = 0; i < confSquares.size(); i++) {
                if (confSquares[i] == my_rank - 1) continue;
                else {
                    MPI_Send(&(confSquares[i]), 1, MPI_INT, 0, my_rank, newworld);
                }
            }
            for (int i = 0; i < confSquares.size(); i++) {
                if (confSquares[i] == my_rank - 1) continue;
                else {
                    double *temp = new double[myPieceSize];
                    MPI_Recv(temp, myPieceSize, MPI_DOUBLE, 0, 0, newworld, &status);
                    Xsquares_in_process[confSquares[i]]=temp;
                }
            }
        }
    }
    else{
        for(int i=1; i<firstEmptyProcess; i++){
            int temp_tosend;
            MPI_Recv(&tobesent, 1, MPI_INT, i, i, newworld, &status);
            while(tobesent--){
                MPI_Recv(&temp_tobesent, 1, MPI_INT, i, i, newworld, &status);
                send_Xids.push_back(temp_tobesent);
            }
            for(int j=0; j< send_Xids.size(); j++){
                MPI_Send(x+displs[send_Xids[j]], displs[send_Xids[j]+1]-displs[send_Xids[j]], MPI_DOUBLE, i, 0, newworld);
            }
            send_Xids.clear();
        }
    }
    if(my_rank==firstEmptyProcess-1) cout << my_rank << ": done with sending far neighbours " <<  endl << flush;

    double val;
    vector<double*> Ysquares_in_process;
    if(my_rank) {
        for (int i = 0; i < firstEmptyProcess; i++) {
            double *ptr = nullptr;
            Ysquares_in_process.push_back(ptr);
        }
        for (int i = 0; i < confSquares.size(); i++) {
            Ysquares_in_process[confSquares[i]] = new double[ displs[confSquares[i] + 1] - displs[confSquares[i]] ];
            for (int j = 0; j < displs[confSquares[i] + 1] - displs[confSquares[i]]; j++) Ysquares_in_process[confSquares[i]][j] = 0.0;
        }
    }
    // new row ptr
    int *rowPtr = new int[myPieceSize+1];
    rowPtr[0] = 0;
    for(int i=1; i<myPieceSize+1; i++){
        rowPtr[i] = rowPtr[i-1] + myRowDiff[i-1];
    }

    int colIndModulo, processID, colIndModuloInterP;
    int *colIndOffset, *colIndProcessID, *colIndOffsetInter;
    colIndOffset = new int[myNNZ];
    colIndProcessID = new int[myNNZ];
    colIndOffsetInter = new int[myNNZ];
    for (int i = 0; i < myPieceSize; i++) {
        for (int j = rowPtr[i]; j < rowPtr[i + 1]; j++) {
            colInd = myColInd[j] - 1;
            if(colInd < displs[1]) {
                colIndProcessID[j] = 0;
                colIndOffsetInter[j] = colInd;
                colIndOffset[j] = colInd;
            }
            else {
                colIndProcessID[j] = 1 + (colInd - displs[1]) / remainingPieceSize;
                colIndModulo = (int) fmod(colInd - displs[colIndProcessID[j]],myPieceSize );
                colIndOffsetInter[j] = colIndModulo;
                colIndOffset[j] = colIndModulo;
            }
        }
    }
    cout << my_rank << ": donw with indices for middle " <<  endl << flush;
    // pre-determine outer's offsets
    int *outer_colIndOffset, *outer_colIndProcessID;
    int *outer_rowIndOffset, *outer_rowIndProcessID;
    if(my_rank==0) {
        int rowInd;
        outer_colIndOffset = new int[outer_col.size()];
        outer_colIndProcessID = new int[outer_col.size()];
        outer_rowIndOffset = new int[outer_col.size()];
        outer_rowIndProcessID = new int[outer_col.size()];
        for (int i = 0; i < outer_col.size(); i++) {
            colInd = outer_col[i] - 1;
            rowInd = outer_row[i] - 1;
            if (colInd < displs[1] ) {
                outer_colIndProcessID[i] = 0;
                outer_colIndOffset[i] = colInd;
            } else {
                outer_colIndProcessID[i] = 1 + (colInd - displs[1]) / remainingPieceSize;
                outer_colIndOffset[i] = (int) fmod(colInd - displs[outer_colIndProcessID[i]], remainingPieceSize);
            }
            if (rowInd < displs[1]) {
                outer_rowIndProcessID[i] = 0;
                outer_rowIndOffset[i] = rowInd;
            } else {
                outer_rowIndProcessID[i] = 1 + (rowInd - displs[1]) / remainingPieceSize;
                outer_rowIndOffset[i] = (int) fmod(rowInd - displs[outer_rowIndProcessID[i]], remainingPieceSize);
            }
        }
    }
    cout << my_rank << ": donw with indices for outer " <<  endl << flush;

    double maxTime = 0.0;
    MPI_Barrier(newworld);

    start_time = MPI_Wtime();
    //cout << (bool) MPI_WTIME_IS_GLOBAL << endl;
    for (int i = 0; i < myPieceSize; i++) {
        val = myDiags[i] * myX[i];
        for (int j = rowPtr[i]; j < rowPtr[i+1]; j++) {
            processID = colIndProcessID[j];
            if ((myColInd[j] - 1) <  displs[my_rank]) {
                if (processID == my_rank - 1) {
                    val += myOffDiags[j] * neighbourX[colIndOffsetInter[j]];
                }
                else {
                    val += myOffDiags[j] * (Xsquares_in_process[processID])[colIndOffsetInter[j]];
                }
                (Ysquares_in_process[processID])[colIndOffsetInter[j]] -= myOffDiags[j] * myX[i];
            } else {
                y[colIndOffset[j]] -= myOffDiags[j] * myX[i];
                val += myOffDiags[j] * myX[colIndOffset[j]];
            }
        }
        y[i] += val;
    }
    end_time = MPI_Wtime();
    time += end_time - start_time;
    cout << my_rank << ": computation - multiplication " << end_time - start_time << endl << flush;

    MPI_Win window;
    MPI_Win_create(y, 1*sizeof(double), sizeof(double), MPI_INFO_NULL, newworld, &window);
    MPI_Win_fence(0, window);
    // accumulate y results.
    start_time = MPI_Wtime();
    if (my_rank) {
        for (int i = 0; i < confSquares.size(); i++) {
            {
                if(confSquares[i] == 0 )
                    MPI_Accumulate(Ysquares_in_process[0], displs[1], MPI_DOUBLE, 0, 0,
                                   displs[1], MPI_DOUBLE, MPI_SUM, window);

                else MPI_Accumulate(Ysquares_in_process[confSquares[i]], remainingPieceSize, MPI_DOUBLE, confSquares[i], 0,
                                    remainingPieceSize, MPI_DOUBLE, MPI_SUM, window);
            }
        }
    }
    end_time = MPI_Wtime();
    time += end_time-start_time;
    cout << my_rank << ": communication - accumulation " << end_time - start_time << endl << flush;
    MPI_Win_fence(0, window);
    MPI_Win_free(&window);

    MPI_Reduce(&time, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, newworld);

    MPI_Win window1;
    MPI_Win_create(y, sizeof(double), sizeof(double), MPI_INFO_NULL, newworld, &window1);
    MPI_Win_fence(0, window1);
    if(my_rank == 0) {
        int val1, val2;
        start_time = MPI_Wtime();
        for (int i = 0; i < outer_col.size(); i++) {
            val1 = -outer_val[i] * x[outer_row[i] - 1];  // output[colindoffset]
            val2 = outer_val[i] * x[outer_col[i] - 1]; // output[row]
            MPI_Accumulate(&val1, 1, MPI_DOUBLE, outer_colIndProcessID[i], outer_colIndOffset[i], 1, MPI_DOUBLE,
                           MPI_SUM, window1);
            MPI_Accumulate(&val2, 1, MPI_DOUBLE, outer_rowIndProcessID[i], outer_rowIndOffset[i], 1, MPI_DOUBLE,
                           MPI_SUM, window1);

        }
        end_time = MPI_Wtime();
        maxTime += end_time - start_time;
        cout << my_rank << ": computation - outer " << end_time - start_time << endl << flush;
    }

    MPI_Win_fence(0, window1);
    MPI_Win_free(&window1);
    // gather just to conform code works.
    MPI_Gatherv(y, myPieceSize, MPI_DOUBLE, output, pieceSizeArr, displs, MPI_DOUBLE, 0, newworld);

    if(my_rank == 0){
        //printf("The max of all ranks is %lf.\n", maxTime);
        cout << my_rank << ": Total time " << maxTime << endl << flush;
        cout << "-----------------------------------------------" << endl << flush;

        ofstream myfile;
        myfile.open ("/home/selin/3way-Par-Results/" + matrix_names[inputType] + "/result.txt", ios::out | ios::trunc);
        // cout << "Writing to output... " << endl;
        for (int i=0; i<n; i++) {
            myfile << std::fixed << std::setprecision(dbl::max_digits10) << output[i] << '\t';
        }
        myfile.close();
        //cout << "Completed output... " << endl;


        delete [] x;
        delete[] matrixRowptr;
        delete[] matrixRowDiff;
        delete[] matrixColind;
        delete[] matrixOffDiagonal;
        delete[] matrixDiagonal;
        delete[] NNZs2;
        delete[] pieceSizeArr;
    }


    // Finalize MPI
    // This must always be called after all other MPI functions
    delete [] myX;
    delete [] myColInd;
    delete [] myOffDiags;
    delete [] myDiags;
    delete [] myRowDiff;

    //cout << "Rank: " << my_rank << " Bitti." << endl << flush;

    MPI_Finalize();

    return 0;
}