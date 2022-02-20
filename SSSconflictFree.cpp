//
// Created by Selin Yıldırım on 7.01.2022.
//

#include <iostream>
#include  "header.h"
#include <mpi.h>

using namespace std;

int p;
vector<int> x;
vector<int> y;
vector<int> elm_row; // use this to find out row_ptr values
#define MATRIX_COUNT 6
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
        const fs::path matrixFolder{"/home/selin/Paper-Implementation/CSR-Data/" + matrix_names[z]};
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
            while (myfile >> tempVal) {
                tempVec.push_back(tempVal);
            }

            if(dir_entry.path().stem() == "dvalues"){
                dvaluesPtrs.push_back(new double[tempVec.size()]);
                double *temp = dvaluesPtrs[0];
                for(int i=0; i<tempVec.size(); i++) temp[i]=tempVec[i];
                dvaluesSize.push_back(tempVec.size());
            }
            else if(dir_entry.path().stem() == "nond_values"){
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
int main(int argc, char **argv) {
    int my_rank, world_size;

    // Initialize MPI
    MPI_Init(&argc, &argv);

    // Get the number of processes in MPI_COMM_WORLD
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of this process in MPI_COMM_WORLD
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    bitset<914904> *graph0 = nullptr;
    bitset<923136> *graph1 = nullptr;
    bitset<952208> *graph2 = nullptr;
    bitset<503632> *graph3 = nullptr;
    bitset<1391352> *graph4 = nullptr;
    bitset<943696> *graph5 = nullptr;

    double time1,time2;
    int n, rowLimit, inputType, bitsetSize;

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

        delete [] dvaluesPtrs[0];

        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&inputType, 1, MPI_INT, 0, MPI_COMM_WORLD);

            int lowerLimit, upperLimit, pieceSize, elementCount, offset;
            for(int t=1; t<world_size; t++){
                //std::cout << "Rank:" << my_rank << " will send to: " << t  << endl;
                lowerLimit = rowLimit*(t);
                if (t==world_size-1) {upperLimit=n; pieceSize= n - lowerLimit;}
                else {upperLimit = lowerLimit + rowLimit; pieceSize= rowLimit;}
                // send it to child^th process
                MPI_Send(&pieceSize, 1, MPI_INT, t, 0, MPI_COMM_WORLD);
                //std::cout << "Rank:" << my_rank << " sent first send." << endl;
                elementCount = matrixRowptr[upperLimit] - matrixRowptr[lowerLimit] ;
                MPI_Send(matrixRowptr + lowerLimit, pieceSize+1, MPI_INT, t, 0, MPI_COMM_WORLD);
                offset = (matrixRowptr[lowerLimit]);

                MPI_Send(&elementCount, 1, MPI_INT, t, 0, MPI_COMM_WORLD);
                //std::cout << "Rank:" << my_rank << " sent 2nd send." << endl;

                MPI_Send(matrixColind+ offset, elementCount, MPI_INT, t, 0, MPI_COMM_WORLD);
                //std::cout << "Rank:" << my_rank << " sent 3rd send." << endl;
                MPI_Send(matrixOffDiagonal+ offset, elementCount, MPI_DOUBLE, t, 0, MPI_COMM_WORLD);
                //std::cout << "Rank:" << my_rank << " sent 4th send." << endl;
                //std::cout << "Rank:" << my_rank << " has sent to: " << t  << endl;
                }
            delete [] matrixOffDiagonal;

        std::cout  << "Rank: " << my_rank << " starts computing... " << endl;
        switch(inputType) {
            case 0 : {
                graph0 = new bitset<914904>[rowLimit];
                break;
            }
            case 1 : {
                graph1 = new bitset<923136>[rowLimit];
                break;
            }
            case 2 : {
                graph2 = new bitset<952208>[rowLimit];
                break;
            }
            case 3 : {
                graph3 = new bitset<503632>[rowLimit];
                break;
            }
            case 4 : {
                graph4 = new bitset<1391352>[rowLimit];
                break;
            }
            case 5 : {
                graph5 = new bitset<943696>[rowLimit];
                break;
            }
        }

            //std::cout  << "Rank: " << my_rank << " sent all data" << endl;
            int elmCountPerRow, rowBegin = 0;
            int confCount=0, colInd, rowEnd = rowLimit;
            time1=MPI_Wtime();
            // include remaining rows
            for (int i = 0; i < rowEnd; i++) {
                elmCountPerRow = matrixRowptr[i + 1] - matrixRowptr[i];
                for (int j = 0; j < elmCountPerRow; j++) {
                    colInd = matrixColind[j];
                    if (colInd < rowBegin) {
                        switch(inputType) {
                            case 0 : {
                                   graph0[i].set(colInd);
                                    confCount++;
                                break;
                            }
                            case 1 : {
                                    graph1[i].set(colInd);
                                    confCount++;

                                break;
                            }
                            case 2 : {
                                    graph2[i].set(colInd);
                                    confCount++;

                                break;
                            }
                            case 3 : {
                                   graph3[i].set(colInd);
                                    confCount++;

                                break;
                            }
                            case 4 : {
                                    graph4[i].set(colInd);
                                    confCount++;

                                break;
                            }
                            case 5 : {
                                    graph5[i].set(colInd);
                                    confCount++;

                                break;
                            }
                        }
                    }
                }
            }
        time2=MPI_Wtime();
        std::cout  << "Rank: " << my_rank << " computed # Conflicts: "  << confCount << " Time: " << time2-time1 << endl;

        delete [] graph0;
        delete [] matrixColind;
        delete [] matrixRowptr;

        bitset<914904>** graph0_double;
        bitset<923136>** graph1_double;
        bitset<952208>** graph2_double;
        bitset<503632>** graph3_double;
        bitset<1391352>** graph4_double;
        bitset<943696>** graph5_double;
        switch(inputType) {
            case 0 : {
                graph0_double= new bitset<914904>*[world_size-1];
                break;
            }
            case 1 : {
                graph1_double= new bitset<923136>*[world_size-1];
                break;
            }
            case 2 : {
                graph2_double= new bitset<952208>*[world_size-1];
                break;
            }
            case 3 : {
                graph3_double= new bitset<503632>*[world_size-1];
                break;
            }
            case 4 : {
                graph4_double= new bitset<1391352>*[world_size-1];
                break;
            }
            case 5 : {
                graph5_double= new bitset<943696>*[world_size-1];
                break;
            }
        }

        // receive computed conflict map.
        for(int t=1; t<world_size; t++){
            if (t==world_size-1) pieceSize= n - lowerLimit;
            else pieceSize= rowLimit;
            switch(inputType) {
                case 0 : {
                    graph0_double[t-1] = new bitset<914904>[pieceSize];
                    MPI_Recv(graph0_double[t-1], pieceSize * (914904/8), MPI_CHAR, t, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
                    break;
                }
                case 1 : {
                    graph1_double[t-1] = new bitset<923136>[pieceSize];
                    MPI_Recv(graph1_double[t-1], pieceSize * (923136/8), MPI_CHAR, t, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
                    break;
                }
                case 2 : {
                    graph2_double[t-1] = new bitset<952208>[pieceSize];
                    MPI_Recv(graph2_double[t-1], pieceSize * (952208/8), MPI_CHAR, t, 0, MPI_COMM_WORLD,MPI_STATUSES_IGNORE);
                    break;
                }
                case 3 : {
                    graph3_double[t-1] = new bitset<503632>[pieceSize];
                    MPI_Recv(graph3_double[t-1], pieceSize * (503632/8), MPI_CHAR, t, 0, MPI_COMM_WORLD,MPI_STATUSES_IGNORE);
                    break;
                }
                case 4 : {
                    graph4_double[t-1] = new bitset<1391352>[pieceSize];
                    MPI_Recv(graph4_double[t-1], pieceSize * (1391352/8), MPI_CHAR, t, 0, MPI_COMM_WORLD,MPI_STATUSES_IGNORE);
                    break;
                }
                case 5 : {
                    graph5_double[t-1] = new bitset<943696>[pieceSize];
                    MPI_Recv(graph5_double[t-1], pieceSize * (943696/8), MPI_CHAR, t, 0, MPI_COMM_WORLD,MPI_STATUSES_IGNORE);
                    break;
                }
            }
            std::cout << "Rank:" << my_rank << " has received computed conflict map from rank: " << t  << endl;
        }

        // compute conflict column indexes here
        vector<vector<vector<int>>> conflictNeighbours;
        // communication for each row will increase parallel overhead.
        // therefore, compute row's conflicts here, all at once.
        string trickyVar;
        int indexCount;
        for(int t=1; t<world_size; t++){
            vector<vector<int>> processorConfs;
            if (t==world_size-1) pieceSize= n - rowLimit*(t);
            else pieceSize= rowLimit;
            for(int i=0; i<pieceSize; i++) {
                indexCount=0;
                vector<int> rowConfs;
                switch (inputType) {
                    case 0 : {
                        trickyVar = (string) graph0_double[t - 1][i].to_string();
                        break;
                    }
                    case 1 : {
                        trickyVar = (string) graph1_double[t - 1][i].to_string();
                        break;
                    }
                    case 2 : {
                        trickyVar = (string) graph2_double[t - 1][i].to_string();
                        break;
                    }
                    case 3 : {
                        trickyVar = (string) graph3_double[t - 1][i].to_string();
                        break;
                    }
                    case 4 : {
                        trickyVar = (string) graph4_double[t - 1][i].to_string();
                        break;
                    }
                    case 5 : {
                        trickyVar = (string) graph5_double[t - 1][i].to_string();
                        break;
                    }
                }
                // fetch conflict columns below
                while(trickyVar.size()) {
                    if (trickyVar.back()=='1') rowConfs.push_back(n-1-indexCount);
                    trickyVar.pop_back();
                    indexCount++;
                }
                processorConfs.push_back(rowConfs);
            }
            conflictNeighbours.push_back(processorConfs);
            std::cout << "Rank:" << my_rank << " has computed conflicts of processor rank: " << t  << endl;
        }

    }
    else {

        int pieceSize, elementCount;
        int *piece_rowptr, *piece_colind;
        double *piece_offdiags;
        // ! index piece_rowptr with pieceSize
        // ! index piece_colind,piece_offdiags with elementCount
        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
        rowLimit = n / world_size + 0.5;
        MPI_Bcast(&inputType, 1, MPI_INT, 0, MPI_COMM_WORLD);

        MPI_Recv(&pieceSize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
        piece_rowptr = new int[pieceSize+1];
        MPI_Recv(piece_rowptr, pieceSize+1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);

        MPI_Recv(&elementCount, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
        piece_colind = new int[elementCount];
        piece_offdiags = new double[elementCount];
        MPI_Recv(piece_colind, elementCount, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
        MPI_Recv(piece_offdiags, elementCount, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
        switch(inputType) {
            case 0 : {
                graph0 = new bitset<914904>[pieceSize];
                break;
            }
            case 1 : {
                graph1 = new bitset<923136>[pieceSize];
                break;
            }
            case 2 : {
                graph2 = new bitset<952208>[pieceSize];
                break;
            }
            case 3 : {
                graph3 = new bitset<503632>[pieceSize];
                break;
            }
            case 4 : {
                graph4 = new bitset<1391352>[pieceSize];
                break;
            }
            case 5 : {
                graph5 = new bitset<943696>[pieceSize];
                break;
            }
        }
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
                    switch(inputType) {
                        case 0 : {
                                graph0[i].set(colInd);
                                confCount++;
                            break;
                        }
                        case 1 : {
                            graph1[i].set(colInd);
                                confCount++;
                            break;
                        }
                        case 2 : {
                            graph2[i].set(colInd);
                                confCount++;
                            break;
                        }
                        case 3 : {
                            graph3[i].set(colInd);
                                confCount++;
                            break;
                        }
                        case 4 : {
                              graph4[i].set(colInd);
                                confCount++;
                            break;
                        }
                        case 5 : {
                               graph5[i].set(colInd);
                                confCount++;
                            break;
                        }
                    }
                }
            }
        }
        time2= MPI_Wtime();
        std::cout  << "Rank: " << my_rank << " computed # Conflicts: "  << confCount << " Time: " << time2-time1 << endl;


        switch(inputType) {
            case 0 : {
                MPI_Send(graph0, pieceSize * (914904/8), MPI_CHAR, 0, 0, MPI_COMM_WORLD);
                break;
            }
            case 1 : {
                MPI_Send(graph1, pieceSize * (923136/8), MPI_CHAR, 0, 0, MPI_COMM_WORLD);
                break;
            }
            case 2 : {
                MPI_Send(graph2, pieceSize * (952208/8), MPI_CHAR, 0, 0, MPI_COMM_WORLD);
                break;
            }
            case 3 : {
                MPI_Send(graph3, pieceSize * (503632/8), MPI_CHAR, 0, 0, MPI_COMM_WORLD);
                break;
            }
            case 4 : {
                MPI_Send(graph4, pieceSize * (1391352/8), MPI_CHAR, 0, 0, MPI_COMM_WORLD);
                break;
            }
            case 5 : {
                MPI_Send(graph5, pieceSize * (943696/8), MPI_CHAR, 0, 0, MPI_COMM_WORLD);
                break;
            }
        }
        //std::cout << "Rank:" << my_rank << " has sent to: " << t  << endl;

        delete [] piece_rowptr;
        delete [] piece_colind;
        delete [] piece_offdiags;
    }

    // any ranked process down below :

    switch(inputType) {
        case 0 : {
            if(graph0) delete [] graph0;
            break;
        }
        case 1 : {
            if(graph1) delete [] graph1;
            break;
        }
        case 2 : {
            if(graph2) delete [] graph2;
            break;
        }
        case 3 : {
            if(graph3) delete [] graph3;
            break;
        }
        case 4 : {
            if(graph4) delete [] graph4;
            break;
        }
        case 5 : {
            if(graph5) delete [] graph5;
            break;
        }
    }

    // Finalize MPI
    // This must always be called after all other MPI functions
    MPI_Finalize();

    return 0;
}