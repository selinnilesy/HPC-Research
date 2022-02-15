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

bool sizeCheck(int i){
    return (rowptrSize[0]==matrixSize[i]+1);
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

    bitset<914898> *graph0 = nullptr;
    bitset<923136> *graph1 = nullptr;
    bitset<952203> *graph2 = nullptr;
    bitset<503625> *graph3 = nullptr;
    bitset<1391349> *graph4 = nullptr;
    bitset<943695> *graph5 = nullptr;

    double time1,time2;
    int n, rowLimit, remainderRows;

    if(my_rank==0) {
        cout << "i call readSSSFormat. " << endl;
        //init();
        if(!argv[1]){
            cout << "please provide input matrix index (int): boneS10, Emilia_923, ldoor, af_5_k101, Serena, audikw_1" << endl;
            return -1;
        }
        readSSSFormat(atoi(argv[1]));
        int inputType = atoi(argv[1]);
        if( rowptrSize[0]!=matrixSize[inputType]+1) {
            cout << "size check failed. quitting..." << endl;
            return -1;
        }
        n = matrixSize[inputType];
        double *matrixOffDiagonal = valuesPtrs[0];
        int *matrixColind = colindPtrs[0];
        int *matrixRowptr= rowptrPtrs[0];
        rowLimit = n / world_size + 0.5; // ceiling function for keeping all processes full. only last one have larger
        remainderRows =  n % world_size;

        switch(inputType) {
            case 0 : {
                graph0 = new bitset<914898>[914898]; break;
            }
            case 1 : {
                graph1 = new bitset<923136>[923136]; break;
            }
            case 2 : {
                graph2 = new bitset<952203>[952203]; break;
            }
            case 3 : {
                graph3 = new bitset<503625>[503625]; break;
            }
            case 4 : {
                graph4 = new bitset<1391349>[1391349]; break;
            }
            case 5 : {
                graph5 = new bitset<943695>[943695]; break;
            }
        }
        // double time1 = MPI_Wtime()
        std::cout << "Rank:" << my_rank << " Matrix: " << matrix_names[inputType]  << endl;
        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&inputType, 1, MPI_INT, 0, MPI_COMM_WORLD);

        //MPI_Scatter(matrixRowptr, rowLimit+1, MPI_INT, void* buffer_recv, rowLimit+1, MPI_INT, 0, MPI_COMM_WORLD);

            int lowerLimit, upperLimit, pieceSize, elementCount, offset;
            for(int t=1; t<world_size; t++) {
                //std::cout << "Rank:" << my_rank << " will send to: " << t  << endl;
                offset = matrixRowptr[t*rowLimit];

                if(t==world_size-1) {
                    elementCount = matrixRowptr[n] - matrixRowptr[t*rowLimit];
                    MPI_Send(matrixRowptr + t*rowLimit, rowLimit + remainderRows + 1, MPI_INT, t, 0, MPI_COMM_WORLD);
                }
                else{
                    elementCount = matrixRowptr[t*rowLimit + rowLimit] - matrixRowptr[t*rowLimit];
                    MPI_Send(matrixRowptr + t*rowLimit, rowLimit + 1, MPI_INT, t, 0, MPI_COMM_WORLD);
                }

                MPI_Send(matrixColind + offset, elementCount, MPI_INT, t, 0, MPI_COMM_WORLD);
                MPI_Send(matrixOffDiagonal + offset, elementCount, MPI_DOUBLE, t, 0, MPI_COMM_WORLD);
                std::cout << "Rank:" << my_rank << " has sent to: " << t  << endl;
            }
        bitset<914898> *temp = new bitset<914898>[rowLimit];
        switch (inputType) {
            case 0 : {
                MPI_Scatter(graph0, rowLimit * 914898 / sizeof(char), MPI_CHAR, temp, rowLimit * 914898 / sizeof(char), MPI_CHAR, 0, MPI_COMM_WORLD);
                break;
            }
            case 1 : {
                MPI_Scatter(graph1, rowLimit * 923136 / sizeof(char), MPI_CHAR, graph1, rowLimit * 923136 / sizeof(char), MPI_CHAR, 0, MPI_COMM_WORLD);
                break;
            }
            case 2 : {
                MPI_Scatter(graph2, rowLimit * 952203 / sizeof(char), MPI_CHAR, graph2, rowLimit * 952203 / sizeof(char), MPI_CHAR, 0, MPI_COMM_WORLD);
                break;
            }
            case 3 : {
                MPI_Scatter(graph3, rowLimit * 503625 / sizeof(char), MPI_CHAR, graph3, rowLimit * 503625 / sizeof(char), MPI_CHAR, 0, MPI_COMM_WORLD);
                break;
            }
            case 4 : {
                MPI_Scatter(graph4, rowLimit * 1391349 / sizeof(char), MPI_CHAR, graph4, rowLimit * 1391349 / sizeof(char), MPI_CHAR, 0, MPI_COMM_WORLD);
                break;
            }
            case 5 : {
                MPI_Scatter(graph5, rowLimit * 943695 / sizeof(char), MPI_CHAR, graph5, rowLimit * 943695 / sizeof(char), MPI_CHAR, 0, MPI_COMM_WORLD);
                break;
            }
        }
        std::cout << "Rank:" << my_rank << " scatter ok. Now starts last part send" << endl;
                switch (inputType) {
                    case 0 : {
                        MPI_Send(graph0 + (world_size - 1) * rowLimit + rowLimit,
                                 remainderRows * 914898 / sizeof(char), MPI_CHAR,
                                 world_size - 1, 0, MPI_COMM_WORLD);
                        break;
                    }
                    case 1 : {
                        MPI_Send(graph1 + (world_size - 1) * rowLimit + rowLimit, remainderRows * 923136 / sizeof(char),
                                 MPI_CHAR, world_size - 1, 0, MPI_COMM_WORLD);
                        break;
                    }
                    case 2 : {
                        MPI_Send(graph2 + (world_size - 1) * rowLimit + rowLimit, remainderRows * 952203 / sizeof(char),
                                 MPI_CHAR, world_size - 1, 0, MPI_COMM_WORLD);
                        break;
                    }
                    case 3 : {
                        MPI_Send(graph3 + (world_size - 1) * rowLimit + rowLimit, remainderRows * 503625 / sizeof(char),
                                 MPI_CHAR, world_size - 1, 0, MPI_COMM_WORLD);
                        break;
                    }
                    case 4 : {
                        MPI_Send(graph4 + (world_size - 1) * rowLimit + rowLimit, remainderRows * 1391349 / sizeof(char),
                                 MPI_CHAR, world_size - 1, 0, MPI_COMM_WORLD);
                        break;
                    }
                    case 5 : {
                        MPI_Send(graph5 + (world_size - 1) * rowLimit + rowLimit, remainderRows * 943695 / sizeof(char),
                                 MPI_CHAR, world_size - 1, 0, MPI_COMM_WORLD);
                        break;
                    }
                }
            //std::cout  << "Rank: " << my_rank << " sent all data" << endl;
            int elmCountPerRow, rowBegin = 0;
            int confCount=0, colInd, rowEnd= rowLimit;
            time1=MPI_Wtime();
            //std::cout  << "Rank: " << my_rank << " starts computing... " << endl;
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
        switch(inputType) {
            case 0 : {
                delete [] graph0;
                break;
            }
            case 1 : {
                delete [] graph1;
                break;
            }
            case 2 : {
                delete [] graph2;
                break;
            }
            case 3 : {
                delete [] graph3;
                break;
            }
            case 4 : {
                delete [] graph4;
                break;
            }
            case 5 : {
                delete [] graph5;
                break;
            }
        }
            delete [] matrixOffDiagonal;
            delete [] matrixColind;
            delete [] matrixRowptr;

    }
    else {

        int inputType, elementCount;
        int *piece_rowptr, *piece_colind;
        double *piece_offdiags;
        // ! index piece_rowptr with pieceSize
        // ! index piece_colind,piece_offdiags with elementCount
        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&inputType, 1, MPI_INT, 0, MPI_COMM_WORLD);
        rowLimit = n / world_size + 0.5;
        remainderRows = n % world_size;

        if(my_rank == world_size-1){
            piece_rowptr = new int[rowLimit+remainderRows+1];
            MPI_Recv(piece_rowptr, rowLimit+remainderRows+1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
            elementCount = piece_rowptr[rowLimit+remainderRows] - piece_rowptr[0];
        }
        else{
            piece_rowptr = new int[rowLimit+1];
            MPI_Recv(piece_rowptr, rowLimit+1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
            elementCount = piece_rowptr[rowLimit] - piece_rowptr[0];
        }
        piece_colind = new int[elementCount];
        piece_offdiags = new double[elementCount];
        MPI_Recv(piece_colind, elementCount, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
        MPI_Recv(piece_offdiags, elementCount, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
        std::cout << "Rank:" << my_rank << " received until scatter " << endl;

        switch(inputType) {
            case 0 : {
                if(my_rank!=world_size-1) graph0 = new bitset<914898>[rowLimit];
                else graph0 = new bitset<914898>[rowLimit+ remainderRows];
                MPI_Scatter(graph0, rowLimit * 943695 / sizeof(char), MPI_CHAR, graph0, rowLimit * 943695 / sizeof(char), MPI_CHAR, 0, MPI_COMM_WORLD);
                break;
            }
            case 1 : {
                if(my_rank!=world_size-1) graph1 = new bitset<923136>[rowLimit];
                else graph1 = new bitset<923136>[rowLimit+ remainderRows];
                MPI_Scatter(graph1, rowLimit * 923136 / sizeof(char), MPI_CHAR, graph1, rowLimit * 923136 / sizeof(char), MPI_CHAR, 0, MPI_COMM_WORLD);
                break;
            }
            case 2 : {
                if(my_rank!=world_size-1) graph2 = new bitset<952203>[rowLimit];
                else graph2 = new bitset<952203>[rowLimit+ remainderRows];
                MPI_Scatter(graph2, rowLimit * 952203 / sizeof(char), MPI_CHAR, graph2, rowLimit * 952203 / sizeof(char), MPI_CHAR, 0, MPI_COMM_WORLD);
                break;
            }
            case 3 : {
                if(my_rank!=world_size-1) graph3 = new bitset<503625>[rowLimit];
                else graph3 = new bitset<503625>[rowLimit+ remainderRows];
                MPI_Scatter(graph3, rowLimit * 503625 / sizeof(char), MPI_CHAR, graph3, rowLimit * 503625 / sizeof(char), MPI_CHAR, 0, MPI_COMM_WORLD);
                break;
            }
            case 4 : {
                if(my_rank!=world_size-1) graph4 = new bitset<1391349>[rowLimit];
                else graph4 = new bitset<1391349>[rowLimit+ remainderRows];
                MPI_Scatter(graph4, rowLimit * 1391349 / sizeof(char), MPI_CHAR, graph4, rowLimit * 1391349 / sizeof(char), MPI_CHAR, 0, MPI_COMM_WORLD);
                break;
            }
            case 5 : {
                if(my_rank!=world_size-1) graph5 = new bitset<943695>[rowLimit];
                else graph5 = new bitset<943695>[rowLimit+ remainderRows];
                MPI_Scatter(graph5, rowLimit * 943695 / sizeof(char), MPI_CHAR, graph5, rowLimit * 943695 / sizeof(char), MPI_CHAR, 0, MPI_COMM_WORLD);
                break;
            }
        }
        if(my_rank==world_size-1 && remainderRows) {
            switch (inputType) {
                case 0 : {
                    MPI_Recv(graph0 + rowLimit,
                             remainderRows * 914898 / sizeof(char), MPI_CHAR,
                             0, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
                    break;
                }
                case 1 : {
                    MPI_Recv(graph1 + rowLimit, remainderRows * 923136 / sizeof(char),
                             MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
                    break;
                }
                case 2 : {
                    MPI_Recv(graph2 + rowLimit, remainderRows * 952203 / sizeof(char),
                             MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
                    break;
                }
                case 3 : {
                    MPI_Recv(graph3 + rowLimit, remainderRows * 503625 / sizeof(char),
                             MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
                    break;
                }
                case 4 : {
                    MPI_Recv(graph4 + rowLimit, remainderRows * 1391349 / sizeof(char),
                             MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
                    break;
                }
                case 5 : {
                    MPI_Recv(graph5 + rowLimit, remainderRows * 943695 / sizeof(char),
                             MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
                    break;
                }
            }
        }
        //std::cout  << "Rank: " << my_rank << " received all data with piece-size: " << pieceSize << endl;

        int elmCountPerRow, rowBegin = my_rank*rowLimit;
        int confCount=0, colInd, rowEnd;
        rowEnd = (my_rank==world_size-1 && remainderRows) ? rowLimit+remainderRows : rowLimit;
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
        delete [] piece_rowptr;
        delete [] piece_colind;
        delete [] piece_offdiags;
        switch(inputType) {
            case 0 : {
                delete [] graph0;
                break;
            }
            case 1 : {
                delete [] graph1;
                break;
            }
            case 2 : {
                delete [] graph2;
                break;
            }
            case 3 : {
                delete [] graph3;
                break;
            }
            case 4 : {
                delete [] graph4;
                break;
            }
            case 5 : {
                delete [] graph5;
                break;
            }
        }
    }

    // Finalize MPI
    // This must always be called after all other MPI functions
    MPI_Finalize();

    return 0;
}