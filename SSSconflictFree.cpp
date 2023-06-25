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
int readCSRFormat(int z, int outerBW) {
    fs::path matrixFolder;
    matrixFolder = "/home/selin/Split-Data/" + matrix_names[z] + "/middle/CSR-Data";
    for(auto const& dir_entry: fs::directory_iterator{matrixFolder}){
        std::fstream myfile(dir_entry.path(), std::ios_base::in);
        if(dir_entry.path().stem() == ("1.000000-" + to_string((double) bandwithSize[z]-outerBW) + "-row")) {
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
        else if(dir_entry.path().stem() == ("1.000000-" + to_string( (double) bandwithSize[z]-outerBW) + "-col")) {
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
        else if(dir_entry.path().stem() == ("1.000000-" + to_string((double) bandwithSize[z]-outerBW) + "-val")){
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
        if (dir_entry.path().stem() == ("coordinate-1.000000-" + to_string((double) bandwithSize[z] - outerBW) + "-row")) {
            int tempValInt;
            while (myfile >> tempValInt) {
                outer_row.push_back(tempValInt);
            }
            cout << dir_entry.path() << " has been read with size: " << outer_row.size() << endl;
            myfile.close();
        } else if (dir_entry.path().stem() == ("coordinate-1.000000-" + to_string((double) bandwithSize[z] - outerBW) + "-col")) {
            int tempValInt;
            while (myfile >> tempValInt) {
                outer_col.push_back(tempValInt);
            }
            cout << dir_entry.path() << " has been read with size: " << outer_col.size() << endl;
            myfile.close();
        } else if (dir_entry.path().stem() == ("coordinate-1.000000-" + to_string((double) bandwithSize[z] - outerBW) + "val")) {
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

    double time=0.0, end_time, start_time;
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
        readCSRFormat(atoi(argv[1]), 3);

        n = matrixSize[atoi(argv[1])];
        inputType = atoi(argv[1]);
        pieceSize = ceil(n / world_size);
        nnz = colindSize[0];
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
        for (size_t i = 0; i < n; i++) x[i] = 1.0;
    }
    pieceSizeArr = new int[world_size];
    displs = new int[world_size];
    int approx_nnz;
    if(my_rank==0)  {

        matrixRowDiff = new int[n];
        for (size_t i = 0; i < n; i++) {
            matrixRowDiff[i] = matrixRowptr[i+1] - matrixRowptr[i];
        }
        approx_nnz = nnz/world_size + 0.5; //  ceiling


        int temp_nnz, temp_limit=0;
        for (int i = 0; i < world_size; i++) {
            temp_nnz=0;
            for (size_t j = temp_limit; j < n; j++) {
                temp_nnz += matrixRowDiff[j];
                if(i==world_size-1){
                    if(j == n-1 ){
                        pieceSizeArr[i] = j-temp_limit+1;
                        temp_limit = j+1;
                        break;
                    }
                }
                else {
                    if (temp_nnz >= approx_nnz) {
                        pieceSizeArr[i] = j - temp_limit+1;
                        temp_limit = j+1;
                        break;
                    }
                }
            }
        }
        remainingPieceSize= pieceSizeArr[world_size-1];

        for (int i = 0; i < world_size; i++) {
            cout << "row count for process " << i << ": " << pieceSizeArr[i] << '\n' << flush;
        }

        displs[0] = 0;
        for (int i=1; i<world_size; i++)
            displs[i] = displs[i-1] + pieceSizeArr[i-1];

    }
    myPieceSize = displs[my_rank+1] - displs[my_rank];
    //start_time = MPI_Wtime();
    MPI_Bcast( (void*) displs, world_size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast( (void*) pieceSizeArr, world_size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&remainingPieceSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // scatter pieceSizeArr.
    MPI_Scatter(pieceSizeArr, 1, MPI_INT, &myPieceSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // broadcast expected global piece size.
    MPI_Bcast(&pieceSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    //end_time = MPI_Wtime();
    //time += end_time-start_time;

    assert (myPieceSize!=0);
    myX = new double[myPieceSize];
    myRowDiff = new int[myPieceSize];
    myDiags = new double[myPieceSize];

    //start_time = MPI_Wtime();
    // scatter x.
    MPI_Scatterv(x, pieceSizeArr, displs, MPI_DOUBLE, myX, myPieceSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // scatter matrixRowptr.
    MPI_Scatterv(matrixRowDiff, pieceSizeArr, displs, MPI_INT, myRowDiff, myPieceSize, MPI_INT, 0, MPI_COMM_WORLD);
    // scatter matrixDiagonal.
    MPI_Scatterv(matrixDiagonal, pieceSizeArr, displs, MPI_DOUBLE, myDiags, myPieceSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //end_time = MPI_Wtime();
    //time += end_time-start_time;

    if(my_rank==0) NNZs = new int[world_size];

    myNNZ=0;
    for (int i = 0; i < myPieceSize; i++) {
        myNNZ += myRowDiff[i];
    }

   // start_time = MPI_Wtime();
    MPI_Gather(&myNNZ, 1, MPI_INT, NNZs, 1, MPI_INT, 0, MPI_COMM_WORLD);
    //end_time = MPI_Wtime();
    //time += end_time-start_time;
    cout << my_rank<< ": myNNZ-" << myNNZ << '\n' << flush;

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
    for (size_t i = 0; i < myPieceSize; i++) y[i] = 0.0;

    int neighbourSize;
    double* closingYs;
    double *output;
    if(!my_rank){
        cout << "firstEmptyProcess: " << firstEmptyProcess << '\t' << flush;
        output = new double[n];
        for (size_t i = 0; i < n; i++) output[i]=0.0;
    }
    MPI_Request request;
    MPI_Status status;
    int totalSize=0;
    if(firstEmptyProcess < world_size) {
        MPI_Win window2;
        MPI_Win_create(output, 1*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &window2);
        MPI_Win_fence(0, window2);
        if(my_rank >= firstEmptyProcess){
            // process diagonals before quitting
            start_time = MPI_Wtime();
            for (size_t i = 0; i < (displs[my_rank+1] - displs[my_rank]); i++) y[i] = myDiags[i] * myX[i];
            // accumulate diag results onto output vector.
            MPI_Accumulate(y, (displs[my_rank+1] - displs[my_rank]), MPI_DOUBLE, 0, displs[my_rank],
                           (displs[my_rank+1] - displs[my_rank]), MPI_DOUBLE, MPI_SUM, window2);
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
    //start_time = MPI_Wtime();
    MPI_Bcast( (void*) displs, firstEmptyProcess, MPI_INT, 0, newworld);
    MPI_Scatterv(matrixColind, NNZs2, displs, MPI_INT, myColInd, myNNZ, MPI_INT, 0, newworld);
    // scatter matrixOffDiagonal.
    MPI_Scatterv(matrixOffDiagonal, NNZs2, displs, MPI_DOUBLE, myOffDiags, myNNZ, MPI_DOUBLE, 0, newworld);
    //end_time = MPI_Wtime();
    //time += end_time-start_time;


    displs[0] = 0;
    for (int i=1; i<firstEmptyProcess+1; i++)
        displs[i] = displs[i-1] + pieceSizeArr[i-1];


    // detect conflicts and their square ids.
    vector<int> confSquares;
    int squareIndex=-1;
    int accum_colInd=0, colInd, lastConfSquare=-1;
    for (size_t i = 0; i < myPieceSize; i++) {
        for (size_t j =accum_colInd ; j< accum_colInd + myRowDiff[i]; j++) {
            colInd = myColInd[j] - 1;
            for (int i=0; i<firstEmptyProcess; i++){
                if(colInd >= displs[i] && colInd<displs[i+1]) {
                    if(i!=my_rank) squareIndex = i;
                    break;
                }
            }
            if(lastConfSquare!=squareIndex && find(confSquares.begin(), confSquares.end(), squareIndex) == confSquares.end() ) {
                lastConfSquare = squareIndex;
                confSquares.push_back(squareIndex);
            }
        }
        accum_colInd+=myRowDiff[i];
    }
    for (size_t i = 0; i < confSquares.size(); i++)  cout << "Rank: " << my_rank << " confs: " << confSquares[i] << '\t' << flush;

    double *neighbourX;
    //start_time = MPI_Wtime();
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
    //end_time = MPI_Wtime();
    //time += end_time-start_time;

    //cout << my_rank << ": done with sending neigbours " <<  endl << flush;
    // send ids of requested X pieces to root, for non-neighbour confs.
    int tobesent, temp_tobesent;
    vector<double*> Xsquares_in_process;
    for(int i=0; i<firstEmptyProcess; i++){
        double *ptr = nullptr;
        Xsquares_in_process.push_back(ptr);
    }
    // store first process id then needed square id that needs another X.
    //start_time = MPI_Wtime();
    vector<int> send_Xids;
    if(my_rank){
        tobesent = confSquares.size() -1 ;
        MPI_Send(&tobesent, 1, MPI_INT, 0, my_rank, newworld);
        if(tobesent) {
            for (size_t i = 0; i < confSquares.size(); i++) {
                if (confSquares[i] == my_rank - 1) continue;
                else {
                    cout << my_rank << " sending to: " << confSquares[i] << endl << flush;
                    MPI_Send(&(confSquares[i]), 1, MPI_INT, 0, my_rank, newworld);
                }
            }
            for (size_t i = 0; i < confSquares.size(); i++) {
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
            for(size_t j=0; j< send_Xids.size(); j++){
                MPI_Send(x+displs[send_Xids[j]], displs[send_Xids[j]+1]-displs[send_Xids[j]], MPI_DOUBLE, i, 0, newworld);
            }
            send_Xids.clear();
        }
    }
    //end_time = MPI_Wtime();
    //time += end_time-start_time;

    //cout << my_rank << ": done with sending far neighbours " <<  endl << flush;

    double val;
    // new row ptr
    int *rowPtr = new int[myPieceSize+1];
    rowPtr[0] = 0;
    for(size_t i=1; i<myPieceSize+1; i++){
        rowPtr[i] = rowPtr[i-1] + myRowDiff[i-1];
    }

    int colIndModulo, processID, colIndModuloInterP;
    int *colIndOffset, *colIndProcessID, *colIndOffsetInter;
    colIndOffset = new int[myNNZ];
    colIndProcessID = new int[myNNZ];
    colIndOffsetInter = new int[myNNZ];
    for (size_t i = 0; i < myPieceSize; i++) {
        for (size_t j = rowPtr[i]; j < rowPtr[i + 1]; j++) {
            colInd = myColInd[j] - 1;
            if(colInd < displs[1]) {
                colIndProcessID[j] = 0;
                colIndOffsetInter[j] = colInd;
                colIndOffset[j] = colInd;
            }
            else {
                for (int i=0; i<firstEmptyProcess; i++){
                    if(colInd >= displs[i] && colInd<displs[i+1]) {
                        colIndProcessID[j] = i;
                        break;
                    }
                }
                colIndModulo = (int) fmod(colInd - displs[colIndProcessID[j]],myPieceSize );
                colIndOffsetInter[j] = colIndModulo;
                colIndOffset[j] = colIndModulo;
            }
        }
    }
    //cout << my_rank << ": donw with indices for middle " <<  endl << flush;
    // pre-determine outer's offsets
    int *outer_transposedProcessID, *outer_transposedOffset;
    int *outer_rowIndOffset, *outer_rowIndProcessID;
    if(my_rank==0) {
        int rowInd;
        outer_transposedOffset = new int[outer_col.size()];
        outer_transposedProcessID = new int[outer_col.size()];
        outer_rowIndOffset = new int[outer_col.size()];
        outer_rowIndProcessID = new int[outer_col.size()];
        for (size_t id = 0; id < outer_col.size(); id++) {
            colInd = outer_col[id] - 1;
            rowInd = outer_row[id] - 1;
            for (int i=0; i<firstEmptyProcess; i++){
                if(colInd >= displs[i] && colInd<displs[i+1]) {
                    outer_transposedProcessID[id] = i;
                    break;
                }
            }
            outer_transposedOffset[id] = (int) (colInd - displs[outer_transposedProcessID[id]]);

            for (int i=0; i<firstEmptyProcess; i++){
                    if(rowInd >= displs[i] && rowInd<displs[i+1]) {
                        outer_rowIndProcessID[id] = i;
                        break;
                    }
            }
            outer_rowIndOffset[id] = (int) (rowInd - displs[outer_rowIndProcessID[id]]);
        }
    }
    cout << my_rank << ": done with indices for outer " <<  endl << flush;

    double maxTime = 0.0;
    MPI_Barrier(newworld);

    // do the actual multiplication for non-conflicting region and store results in y's.
    // transpose of non-conflciting region is also handled here.
    // conflicting ones are multiplied and stored in array y to be communicated later on.
    vector<double>  conflicting_mult;
    vector<int> conflicting_targetID, conflicting_targetOffset;
    start_time = MPI_Wtime();
    for (size_t i = 0; i < myPieceSize; i++) {
        val = myDiags[i] * myX[i];
        for (size_t j = rowPtr[i]; j < rowPtr[i+1]; j++) {
            processID = colIndProcessID[j];
            if ((myColInd[j] - 1) <  displs[my_rank]) {
                if (processID == my_rank - 1) {
                    val += myOffDiags[j] * neighbourX[colIndOffsetInter[j]];
                }
                else {
                    val += myOffDiags[j] * (Xsquares_in_process[processID])[colIndOffsetInter[j]];
                }
                conflicting_targetID.push_back(processID);
                conflicting_targetOffset.push_back(colIndOffsetInter[j]);
                conflicting_mult.push_back(myOffDiags[j] * myX[i]);
                // for older solution with mPI_accumualte and it did not work (now P2P send/receives using above vectors).
                //(Ysquares_in_process[processID])[colIndOffsetInter[j]] -= myOffDiags[j] * myX[i];
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
    if(my_rank==1) cout << my_rank << " at [1]:" << y[1]   << endl << flush;
    int confSize = conflicting_targetID.size();


    double *conf_arry_mult = new double[conflicting_mult.size()];
    int *conf_arry_targetID = new int[conflicting_targetID.size()];
    int *conf_arry_targetOffset = new int[conflicting_targetID.size()];
    for(size_t i=0; i<conflicting_mult.size(); i++) conf_arry_mult[i] = conflicting_mult[i];
    for(size_t i=0; i<conflicting_targetID.size(); i++) conf_arry_targetID[i] = conflicting_targetID[i];
    for(size_t i=0; i<conflicting_targetOffset.size(); i++) conf_arry_targetOffset[i] = conflicting_targetOffset[i];


    if(conflicting_mult.size() != conflicting_targetID.size()) cout << "ERROR-conflicting_mult"  << endl << flush;
    if(conflicting_targetOffset.size() != conflicting_targetID.size()) cout << "ERROR-conflicting_targetOffset"  << endl << flush;


    // gather just to confirm code works.
   // start_time = MPI_Wtime();

    //end_time = MPI_Wtime();
    //cout << my_rank << ": communication - gather output vector " << end_time - start_time << endl << flush;
    //time += end_time - start_time;

    // accumulate conflicting results in array y's in each process.
    start_time = MPI_Wtime();
    for(int comm_index=1; comm_index<firstEmptyProcess; comm_index++ ) {
        // send only
        if (my_rank==firstEmptyProcess-1) {
            if(my_rank-comm_index >= 0 ) {
                MPI_Send(&confSize, 1, MPI_INT, my_rank - comm_index, 0, newworld);
                if (confSize) {
                    MPI_Send(conf_arry_mult, conflicting_mult.size(), MPI_DOUBLE, my_rank-comm_index, 0, newworld);
                    MPI_Send(conf_arry_targetID, conflicting_targetID.size(), MPI_INT, my_rank-comm_index, 0, newworld);
                    MPI_Send(conf_arry_targetOffset, conflicting_targetOffset.size(), MPI_INT, my_rank-comm_index, 0, newworld);
                }
            }
        }
        // interm. processes other than root and last
        else if (my_rank && my_rank!=firstEmptyProcess-1) {
            if (my_rank + comm_index < firstEmptyProcess) {
                double *val, mult;
                int *offset, *targetID, offs, id, tempSize;
                MPI_Recv(&tempSize, 1, MPI_INT, my_rank + comm_index, 0, newworld, &status);
                if (tempSize) {
                    val = new double[tempSize];
                    offset = new int[tempSize];
                    targetID = new int[tempSize];
                    MPI_Recv(val, tempSize, MPI_DOUBLE, my_rank + comm_index, 0, newworld, &status);
                    MPI_Recv(targetID, tempSize, MPI_INT, my_rank + comm_index, 0, newworld, &status);
                    MPI_Recv(offset, tempSize, MPI_INT, my_rank + comm_index, 0, newworld, &status);
                    // cout << my_rank << " receiving conflicting from: " << my_rank + comm_index
                    //     << " with data checks val[0] : " << val[0]
                    //     << endl << flush;

                    int ncf = 0;
                    cout << "proc1 confs: " ;
                    for (size_t j = 0; j < tempSize; j++) {
                        if(targetID[j]==my_rank) {
                            y[offset[j]] -= val[j];
                            ncf++;
                            if(my_rank ==1 && offset[j]==1) cout << -val[j] << '\t' << flush;
                        }
                    }
                    //cout << endl << my_rank << " processed " << ncf << " conflicts"  << endl << flush;

                    delete[] val;
                    delete[] offset;
                    delete[] targetID;
                }
            }
            // both does receive and send
            if(my_rank-comm_index >= 0) {
                MPI_Send(&confSize, 1, MPI_INT, my_rank-comm_index, 0, newworld);
                if (confSize) {
                    MPI_Send(conf_arry_mult, conflicting_mult.size(), MPI_DOUBLE, my_rank-comm_index, 0, newworld);
                    MPI_Send(conf_arry_targetID, conflicting_targetID.size(), MPI_INT, my_rank-comm_index, 0, newworld);
                    MPI_Send(conf_arry_targetOffset, conflicting_targetOffset.size(), MPI_INT, my_rank-comm_index, 0, newworld);
                }
            }
        }
        // just receive
        else if(my_rank==0) {
            if (my_rank + comm_index < firstEmptyProcess) {
                double *val, mult;
                int *offset, *targetID, offs, id, tempSize;
                MPI_Recv(&tempSize, 1, MPI_INT, my_rank + comm_index, 0, newworld, &status);
                if (tempSize) {
                    val = new double[tempSize];
                    offset = new int[tempSize];
                    targetID = new int[tempSize];
                    MPI_Recv(val, tempSize, MPI_DOUBLE, my_rank + comm_index, 0, newworld, &status);
                    MPI_Recv(targetID, tempSize, MPI_INT, my_rank + comm_index, 0, newworld, &status);
                    MPI_Recv(offset, tempSize, MPI_INT, my_rank + comm_index, 0, newworld, &status);
                    cout << my_rank << " receiving conflicting from: " << my_rank + comm_index
                         << " with data checks val[0] : " << val[0]
                         << endl << flush;
                    int ncf = 0;
                    for (size_t j = 0; j < tempSize; j++) {
                        if(targetID[j]==my_rank) {
                            y[offset[j]] -= val[j];
                            ncf++;
                        }
                    }
                    cout << my_rank << " processed " << ncf << " conflicts"  << endl << flush;
                    delete[] val;
                    delete[] offset;
                    delete[] targetID;
                }
            }
        }
    }
    end_time = MPI_Wtime();
    time += end_time-start_time;
    cout << my_rank << ": communication - conflicts " << end_time - start_time << endl << flush;

    MPI_Gatherv(y, myPieceSize, MPI_DOUBLE, output, pieceSizeArr, displs, MPI_DOUBLE, 0, newworld);

    // accumulate serially outer region's sparse results.
    if(my_rank == 0) {
        int val1, val2;
        cout << my_rank << "-size of outer " << outer_col.size() << endl << flush;
        start_time = MPI_Wtime();
        for (size_t i = 0; i < outer_col.size(); i++) {
            val1 = -outer_val[i] * x[outer_row[i] - 1];  // transposed output
            val2 = outer_val[i] * x[outer_col[i] - 1];
            output[displs[outer_transposedProcessID[i]] + outer_transposedOffset[i] ]  += val1;
            output[displs[outer_rowIndProcessID[i]] + outer_rowIndOffset[i] ]  += val2;
            cout << my_rank << ": computing output elm with processID " << outer_rowIndProcessID[i] << endl << flush;
            cout << my_rank << ": computing output elm with transposed processID " << outer_transposedProcessID[i] << endl << flush;

        }
        end_time = MPI_Wtime();
        time += end_time - start_time;
        //cout << my_rank << ": computation - outer " << end_time - start_time << endl << flush;
    }

    MPI_Reduce(&time, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, newworld);

    if(my_rank == 0){
        //printf("The max of all ranks is %lf.\n", maxTime);
        cout << "-----------------------------------------------" << endl << flush;
        cout << my_rank << ": Total time " << maxTime << endl << flush;

        ofstream myfile;
        myfile.open ("/home/selin/3way-Par-Results/" + matrix_names[inputType] + "/result.txt", ios::out | ios::trunc);
        // cout << "Writing to output... " << endl;
        for (size_t i=0; i<n; i++) {
            myfile << std::fixed << std::setprecision(dbl::max_digits10) << output[i] << '\t';
        }
        myfile.close();
        //cout << "Completed output... " << endl;

        // Clean-up root specific data
        delete [] x;
        delete[] matrixRowptr;
        delete[] matrixRowDiff;
        delete[] matrixColind;
        delete[] matrixOffDiagonal;
        delete[] matrixDiagonal;
        delete[] NNZs2;
        delete[] pieceSizeArr;
    }


    // Clean-up all Processes and Finalize MPI
    delete [] myX;
    delete [] myColInd;
    delete [] myOffDiags;
    delete [] myDiags;
    delete [] myRowDiff;

    //cout << "Rank: " << my_rank << " Exiting" << endl << flush;

    MPI_Finalize();

    return 0;
}