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
            if(dir_entry.path().stem() == "coordinate-row") {
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
            else if(dir_entry.path().stem() == "coordinate-col") {
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
            else if(dir_entry.path().stem() == "coordinate-val"){
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

    /*
    cout << "previous boneS10-CSRout_vals:" <<  endl << endl;
    fstream myfile("/home/selin/CSR-Data/" + matrix_names[atoi(argv[1])] + "/banded/CSRout_val.txt", std::ios_base::in);
    double x;
    for(int i=0; i<50; i++){
        myfile >> x;
        cout << x << '\t' ;
    }
    myfile.close();
    cout << "boneS10-CSRout_vals:" <<  endl << endl;
    myfile.open("/home/selin/HPC-Research/boneS10-CSRout_val.txt", std::ios_base::in);
    for(int i=0; i<50; i++){
        myfile >> x;
        cout << x << '\t' ;
    }
    myfile.close();

    cout << "Coordinate vals:" << endl << endl;
    myfile.open("/home/selin/CSR-Data/boneS10/coordinate-val.txt", std::ios_base::in);
    for(int i=0; i<50; i++){
        myfile >> x;
        cout << x << '\t' ;
    }
    myfile.close();
    return 0 ;
     */

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


    int *banded_coordRow = coord_row;
    int *banded_coordCol = coord_col;
    double *banded_coordval = coord_val;

    int *banded_csrRow = new int[matrixSize[inputType]+1];
    int *banded_csrCol = new int[nonzerosSize[inputType]];
    double *banded_csrval = new double[nonzerosSize[inputType]];

    int nnz = nonzerosSize[inputType];
    int nrow=matrixSize[inputType];

    std::cout  <<  " starts computing coocsr... " << endl;
    coocsr_(&nrow,  &nnz, banded_coordval, banded_coordRow, banded_coordCol, banded_csrval, banded_csrCol, banded_csrRow);
    std::cout  <<  " FINISHED computing coocsr... " << banded_csrval[10] << " " << banded_csrCol[10] << " " << banded_csrRow[10] << endl;

    ofstream myfile1, myfile2, myfile3;
    myfile1.open ("/home/selin/CSR-Data/" + matrix_names[inputType] + "/banded/CSRout_row.txt", ios::out | ios::trunc);
    myfile2.open ("/home/selin/CSR-Data/" + matrix_names[inputType] + "/banded/CSRout_col.txt", ios::out | ios::trunc);
    myfile3.open ("/home/selin/CSR-Data/" + matrix_names[inputType] + "/banded/CSRout_val.txt", ios::out | ios::trunc);

    cout << "Writing to " << "-CSRout_row.txt"  << endl;
    for (int i=0; i<matrixSize[inputType]+1; i++) {
        myfile1 << banded_csrRow[i] << '\t';
    }
    cout << "Writing to " << "-CSRout_col.txt"  << endl;
    for (int i=0; i<nonzerosSize[inputType]; i++) {
        myfile2 << banded_csrCol[i] << '\t';
    }
    cout << "Writing to " << "-CSRout_val.txt"  << endl;
    for (int i=0; i<nonzerosSize[inputType]; i++){
        myfile3 << banded_csrval[i] << '\t';
    }
    myfile1.close();
    myfile2.close();
    myfile3.close();
}

/*
int main(int argc, char **argv) {
    int my_rank, world_size;

    // Initialize MPI
    MPI_Init(&argc, &argv);

    // Get the number of processes in MPI_COMM_WORLD
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of this process in MPI_COMM_WORLD
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    double time1,time2;
    int n, rowLimit;

    if(my_rank==0) {
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
        rowLimit = n / world_size + 0.5; // ceiling function

        // double time1 = MPI_Wtime()
        std::cout << "Rank:" << my_rank << " Matrix: " << matrix_names[inputType]  << endl;

        int lowerLimit, upperLimit, pieceSize, elementCount, offset;
        for(int t=1; t<world_size; t++){
            //std::cout << "Rank:" << my_rank << " will send to: " << t  << endl;
            lowerLimit = rowLimit*(t);
            upperLimit = min( lowerLimit + rowLimit , n );
            pieceSize= upperLimit - lowerLimit;
            // send it to child^th process
            MPI_Send(&n, 1, MPI_INT, t, 0, MPI_COMM_WORLD);
            MPI_Send(&rowLimit, 1, MPI_INT, t, 0, MPI_COMM_WORLD);

            MPI_Send(&inputType, 1, MPI_INT, t, 0, MPI_COMM_WORLD);
            MPI_Send(&pieceSize, 1, MPI_INT, t, 0, MPI_COMM_WORLD);

            elementCount = matrixRowptr[upperLimit] - matrixRowptr[lowerLimit] ;
            MPI_Send(matrixRowptr + lowerLimit, pieceSize+1, MPI_INT, t, 0, MPI_COMM_WORLD);
            offset = (matrixRowptr[lowerLimit]);

            MPI_Send(&elementCount, 1, MPI_INT, t, 0, MPI_COMM_WORLD);

            MPI_Send(matrixColind+ offset, elementCount, MPI_INT, t, 0, MPI_COMM_WORLD);
            MPI_Send(matrixOffDiagonal+ offset, elementCount, MPI_DOUBLE, t, 0, MPI_COMM_WORLD);

            //std::cout << "Rank:" << my_rank << " has sent to: " << t  << endl;
        }
        //std::cout  << "Rank: " << my_rank << " sent all data" << endl;
        int elmCountPerRow, rowBegin = 0;
        int confCount=0, colInd, rowEnd = rowLimit;
        time1=MPI_Wtime();
        //std::cout  << "Rank: " << my_rank << " starts computing... " << endl;
        // include remaining rows

        int *banded_coordRow = coord_row;
        int *banded_coordCol = coord_col;
        double *banded_coordval = coord_val;

        std::cout  <<  " starts computing coocsr... " << endl;
        int *banded_csrRow = new int[matrixSize[inputType]+1];
        int *banded_csrCol = new int[nonzerosSize[inputType]];
        double *banded_csrval = new double[nonzerosSize[inputType]];
        // output
        int nnz = nonzerosSize[inputType];
        int nrow=matrixSize[inputType];
        coocsr_(&nrow,  &nnz, banded_coordval, banded_coordRow, banded_coordCol, banded_csrval, banded_csrCol, banded_csrRow);
        std::cout  <<  " FINISHED computing coocsr... " << banded_csrval[10] << " " << banded_csrCol[10] << " " << banded_csrRow[10] << endl;



        for (int i = 0; i < rowEnd; i++) {
            elmCountPerRow = matrixRowptr[i + 1] - matrixRowptr[i];
            for (int j = 0; j < elmCountPerRow; j++) {
                colInd = matrixColind[j];
                if (colInd < rowBegin) { confCount++;}
            }
        }
        /*
        vector<pair<int,int>> pieceVector;
        for(int t=1; t<world_size; t++){
            vector<pair<int,int>> pieceVector_temp;
            MPI_Recv(&nwSize, 1, MPI_INT, t, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            pieceVector_temp.resize(nwSize);
            MPI_Recv((void*)  pieceVector_temp.data(), sizeof(pair<int,int>) * nwSize, MPI_BYTE, t, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            std::copy(pieceVector_temp.begin(), pieceVector_temp.begin()+nwSize , back_inserter(pieceVector));
            if(t==1) continue;
            cout << "Rank 0 received computed vectors from: " << t << endl;
            double t1 = MPI_Wtime();
            removeDuplicateEdges(pieceVector, pieceVector.size());
            double t2 = MPI_Wtime();
            printf( "Elapsed time for joint vectors cleansing: %f\n", t2 - t1 );
        }

        time2=MPI_Wtime();
        std::cout  << "Rank: " << my_rank << " computed # Conflicts: "  << confCount << " Time: " << time2-time1 << endl;


        delete [] matrixOffDiagonal;
        delete [] matrixColind;
        delete [] matrixRowptr;

    }
    else {

        int inputType, pieceSize, elementCount;
        int *piece_rowptr, *piece_colind;
        double *piece_offdiags;
        // ! index piece_rowptr with pieceSize
        // ! index piece_colind,piece_offdiags with elementCount
        MPI_Recv(&n, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
        MPI_Recv(&rowLimit, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);


        MPI_Recv(&inputType, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
        MPI_Recv(&pieceSize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);

        piece_rowptr = new int[pieceSize+1];
        MPI_Recv(piece_rowptr, pieceSize+1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);

        MPI_Recv(&elementCount, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
        piece_colind = new int[elementCount];
        piece_offdiags = new double[elementCount];
        MPI_Recv(piece_colind, elementCount, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
        MPI_Recv(piece_offdiags, elementCount, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);

        //std::cout  << "Rank: " << my_rank << " received all data with piece-size: " << pieceSize << endl;

        int elmCountPerRow, rowBegin = my_rank*rowLimit;
        int confCount=0, colInd, rowEnd = pieceSize;
        // include remaining rows
        //std::cout  << "Rank: " << my_rank << " starts computing... " << endl;
        time1= MPI_Wtime();
        for (int i = 0; i < rowEnd; i++) {
            elmCountPerRow = piece_rowptr[i + 1] - piece_rowptr[i];
            for (int j = 0; j < elmCountPerRow; j++) {
                colInd = piece_colind[j];
                if (colInd < rowBegin) {
                    confCount++;
                }
            }
        }
        time2= MPI_Wtime();
        std::cout  << "Rank: " << my_rank << " computed # Conflicts: "  << confCount << " Time: " << time2-time1 << endl;
        delete [] piece_rowptr;
        delete [] piece_colind;
        delete [] piece_offdiags;
    }

    // Finalize MPI
    // This must always be called after all other MPI functions
    MPI_Finalize();

    return 0;
}
*/