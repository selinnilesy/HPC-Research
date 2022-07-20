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

    // Obtain the group of processes in the world communicator
    MPI_Group world_group, new_group;
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);
    MPI_Comm newworld;

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

    int firstEmptyProcess = world_size;
    if(my_rank==0) {
        for (int i=1; i<world_size; i++) {
            displs[i] = displs[i - 1] + NNZs[i - 1];
            if( NNZs[i] == 0) {
                if(firstEmptyProcess > i) firstEmptyProcess = i;
            }
        }
    }
    MPI_Bcast(&firstEmptyProcess, 1, MPI_INT, 0, MPI_COMM_WORLD);


    if(firstEmptyProcess < world_size) {
        int ranges[1][3] = {{firstEmptyProcess, world_size - 1, 1}};
        MPI_Group_range_excl(world_group, 1, ranges, &new_group);

        // Create a new communicator
        MPI_Comm_create(MPI_COMM_WORLD, new_group, &newworld);
        if (newworld == MPI_COMM_NULL) {
            MPI_Finalize();
            exit(0);
        }
    }
    else{
        MPI_Comm_create(MPI_COMM_WORLD, world_group, &newworld);
    }
    MPI_Comm_rank(newworld, &my_rank);
    cout << "my rank: " << my_rank << endl;
    //assert (myNNZ!=0);
    myColInd = new int[myNNZ];
    myOffDiags = new double[myNNZ];
    // scatter matrixColind.
    MPI_Scatterv(matrixColind, NNZs, displs, MPI_INT, myColInd, myNNZ, MPI_INT, 0, newworld);
    // scatter matrixOffDiagonal.
    MPI_Scatterv(matrixOffDiagonal, NNZs, displs, MPI_DOUBLE, myOffDiags, myNNZ, MPI_DOUBLE, 0, newworld);

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
                    double *temp = new double[pieceSize];
                    MPI_Recv(temp, pieceSize, MPI_DOUBLE, 0, 0, newworld, &status);
                    cout << "Rank: " << my_rank << " received X[" << confSquares[i] << "]" << endl;
                    Xsquares_in_process[confSquares[i]]=temp;
                }
            }
        }
    }
    else{
        for(int i=1; i<firstEmptyProcess; i++){
            MPI_Recv(&tobesent, 1, MPI_INT, i, i, newworld, &status);
            while(tobesent--){
                send_Xids.push_back(i);
                MPI_Recv(&temp_tobesent, 1, MPI_INT, i, i, newworld, &status);
                send_Xids.push_back(temp_tobesent);
            }
        }
        for (int i = 0; i < send_Xids.size(); i+=2){
            MPI_Send(x+send_Xids[i+1]*pieceSize, pieceSize, MPI_DOUBLE, send_Xids[i], 0, newworld);
            cout << "Rank: " << my_rank << " sent X[" << send_Xids[i+1] << "] to process: " << send_Xids[i] << endl;
        }
    }


    double val;
    int accumIndex=0;
    vector<double*> Ysquares_in_process;
    if(my_rank) {
        for (int i = 0; i < firstEmptyProcess; i++) {
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
                if(colInd / pieceSize == my_rank-1) {
                    val += myOffDiags[j] * neighbourX[colIndModulo];
                    if(my_rank==3 && i==0)cout << "my rank: " << my_rank << " - accumulating for y[238050] " << myOffDiags[j] << " by " << neighbourX[colIndModulo]  << endl;
                }
                else{
                    val += myOffDiags[j] * (Xsquares_in_process[colInd / pieceSize])[colIndModulo];
                    if(my_rank==3 && i==0) cout << "my rank: " << my_rank << " -- accumulating for y[238050] "  << myOffDiags[j] << " by " << (Xsquares_in_process[colInd / pieceSize])[colIndModulo] << endl;
                }
                (Ysquares_in_process[colInd / pieceSize])[colIndModulo] -= myOffDiags[j] * myX[i];
                //if(my_rank==3 &&  colInd==0) cout << "my rank: " << my_rank << " computed transposed y " << colInd / pieceSize << " at " << colIndModulo << " " <<  myOffDiags[j] * myX[i] << " by: " << myOffDiags[j]  << " x " <<  myX[i] << endl;
            }
            else{
                y[colIndModulo] -= myOffDiags[j] * myX[i];
                val += myOffDiags[j] * myX[colIndModulo];
                if(my_rank==3 && i==0) cout << "my rank: " << my_rank << " --- accumulating for y[238050] " << myOffDiags[j] << " by " << myX[colIndModulo] << endl;
            }
        }
        y[i] += val;
        if(my_rank==3 && i==0){
            cout << "my rank: " << my_rank << " result for y[238050] " << val << endl;
        }
        accumIndex+=myRowDiff[i];
    }
    double end_time = MPI_Wtime();
    //if(!my_rank) printf("It took me %f seconds for parallel run.\n", end_time-start_time);
    //if(my_rank==1) cout << "my rank: " << my_rank << " computed : " << y[0]<< endl;

    MPI_Win window;
    MPI_Win_create(y, pieceSize*sizeof(double), sizeof(double), MPI_INFO_NULL, newworld, &window);
    MPI_Win_fence(0, window);
    // accumulate y results.
    if(my_rank) {
        for (int i = 0; i < confSquares.size(); i++) {
            {
                MPI_Accumulate(Ysquares_in_process[confSquares[i]], pieceSize, MPI_DOUBLE, confSquares[i], 0, pieceSize,
                               MPI_DOUBLE, MPI_SUM, window);
                //if(confSquares[i]==1) cout << "my rank: " << my_rank << " is sending : " << Ysquares_in_process[confSquares[i]][0]<< endl;
            }
        }
    }
    MPI_Win_fence(0, window);

    // Destroy the window
    MPI_Win_free(&window);


    double *output;
    if(my_rank==0) {
        output = new double[n];
        displs[0] = 0;
        for (int i=1; i<firstEmptyProcess; i++)
            displs[i] = displs[i-1] + pieceSizeArr[i-1];
    }

    MPI_Gatherv(y, myPieceSize, MPI_DOUBLE, output, pieceSizeArr, displs, MPI_DOUBLE, 0, newworld);

    // termination
    if(my_rank==0) {
        // process outer region
        for (int i = 0; i < outer_col.size(); i++) {
            output[outer_col[i]-1] -= outer_val[i] * x[outer_row[i] -1];
            output[outer_row[i]-1] += outer_val[i] * x[outer_col[i] -1];
            if(outer_col[i]-1 == 686171) cout << "outer val:" << -outer_val[i] * x[outer_row[i] -1]<< " added. " << endl;
            else if(outer_row[i]-1 == 686171) cout << "outer val:" << outer_val[i] * x[outer_col[i] -1] << " added." << endl;

        }

        ofstream myfile;
        myfile.open ("/home/selin/3way-Par-Results/" + matrix_names[inputType] + "/result.txt", ios::out | ios::trunc);
        cout << "Writing to output... " << endl;
        for (int i=0; i<n; i++) {
            myfile << std::fixed << std::setprecision(dbl::max_digits10) << output[i] << '\t';
        }
        myfile.close();
        cout << "Completed output... " << endl;

        delete [] x;
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