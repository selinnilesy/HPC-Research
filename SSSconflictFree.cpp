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
    for(int i=0; i<6; i++) {
        valuesPtrs.push_back(new double);
        dvaluesPtrs.push_back(new double);
        colindPtrs.push_back(new int);
        rowptrPtrs.push_back(new int);
    }
}

struct FooHasher {
    size_t operator()(const pair<int, int>&) const {
        return 1;
    }
};
struct Equal_to {
    bool operator()(const pair<int, int> &lhs, const pair<int, int> &rhs) const {
        return (lhs.first == rhs.first && lhs.second == rhs.second) || (lhs.second == rhs.first && lhs.first == rhs.second);
    }
};
void removeDuplicateEdges(std::vector<pair<int, int>> &v, int vecsize)
{
    std::vector<pair<int, int>>::iterator itr = v.begin();
    std::unordered_set<pair<int, int>, FooHasher, Equal_to> s;

    for (auto curr = v.begin(); curr != v.begin()+vecsize; ++curr)
    {
        if (s.insert(*curr).second) {
            *itr++ = pair<int,int>(curr->first, curr->second);
        }
    }

    v.erase(itr, v.begin()+vecsize);
}
int readSSSFormat() {
    double tempVal;
    vector<double> tempVec;
    for (int i=0; i<1; i++) {
        const fs::path matrixFolder{"/home/selin/Paper-Implementation/CSR-Data/" + matrix_names[i]};
        for(auto const& dir_entry: fs::directory_iterator{matrixFolder}){
            std::fstream myfile(dir_entry.path(), std::ios_base::in);
            if(dir_entry.path().stem() == "rowptr") {
                int tempValInt;
                vector<int> tempVecInt;
                while (myfile >> tempValInt) {
                    tempVecInt.push_back(tempValInt);
                }
                rowptrPtrs[i] = new int[tempVecInt.size()];
                int *temp = rowptrPtrs[i];
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
                colindPtrs[i] = new int[tempVecInt.size()];
                int *temp = colindPtrs[i];
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
                dvaluesPtrs[i] = new double[tempVec.size()];
                double *temp = dvaluesPtrs[i];
                for(int i=0; i<tempVec.size(); i++) temp[i]=tempVec[i];
                dvaluesSize.push_back(tempVec.size());
            }
            else if(dir_entry.path().stem() == "nond_values"){
                valuesPtrs[i] = new double[tempVec.size()];
                double *temp = valuesPtrs[i];
                for(int i=0; i<tempVec.size(); i++) temp[i]=tempVec[i];
                valuesSize.push_back(tempVec.size());
            }
            else cout << "unexpected file name: " << dir_entry.path() << endl;
            cout << dir_entry.path() << " has been read." << endl;
            tempVec.clear();
            myfile.close();
        }
    }
    return 0;
}
int main(int argc, char **argv) {
    int my_rank,  world_size;

    // Initialize MPI
    // This must always be called before any other MPI functions
    MPI_Init(&argc, &argv);

    // Get the number of processes in MPI_COMM_WORLD

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of this process in MPI_COMM_WORLD

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    unsigned long long size;
    if(my_rank==0) {
        cout << "i call readSSSFormat. " << endl;
        init();
        readSSSFormat();
        /*for(int i=0; i<MATRIX_COUNT; i++){
            size =rowptrSize[i];
            std::cout << "Rank: " << my_rank << "Size: " << size << std::endl;
            MPI_Bcast(&size, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
            MPI_Bcast(rowptrPtrs[i], size, MPI_INT, 0, MPI_COMM_WORLD);
            size =colindSize[i];
            // std::cout << "Rank: " << my_rank << "Size: " << size << std::endl;
            MPI_Bcast(&size, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
            MPI_Bcast(colindPtrs[i], size, MPI_INT, 0, MPI_COMM_WORLD);
            size =valuesSize[i];
            // std::cout << "Rank: " << my_rank << "Size: " << size << std::endl;
            MPI_Bcast(&size, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
            MPI_Bcast(valuesPtrs[i], size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            size =dvaluesSize[i];
            // std::cout << "Rank: " << my_rank << "Size: " << size << std::endl;
            MPI_Bcast(&size, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
            MPI_Bcast(dvaluesPtrs[i], size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }*/
        double *matrixOffDiagonal;
        int *matrixColind, *matrixRowptr;
        vector <pair<int, int>> globalConf_BoneS10, globalConf_Emilia, globalConf_ldoor, globalConf_af_5_k101, globalConf_Serena, globalConf_audikw_1;
        vector <pair<int, int>> *vecPtr = new vector <pair<int, int>>[6];
        vecPtr[0]=globalConf_BoneS10;
        vecPtr[1]=globalConf_Emilia;
        vecPtr[2]=globalConf_ldoor;
        vecPtr[3]=globalConf_af_5_k101;
        vecPtr[4]=globalConf_Serena;
        vecPtr[5]=globalConf_audikw_1;

        int rowLimit;
        vector<pair<int, int>>::iterator it, it_symmetric;
        for(int z=0; z<1; z++) {
            bool **graph = new bool*[dvaluesSize[z]];
            for(int x = 0; x < dvaluesSize[z]; ++x)
                graph[x] = new bool[dvaluesSize[z]]();
            //vector<pair<int, int> > conflicts;
            matrixOffDiagonal = valuesPtrs[z];
            matrixColind = colindPtrs[z];
            matrixRowptr = rowptrPtrs[z];
            // matrix size : boneS10_DiagonalSize x boneS10_DiagonalSize
            rowLimit = (dvaluesSize[z]) / 4;
            int rowInd, colInd;
            int global_OffDCount = 0;
            int confSize=0;
            std::cout << "Rank: " << my_rank << "Matrix: " << matrix_names[z] << " rowLimit: " << rowLimit << endl;
            for(int k=0; k<world_size; k++) {
                int rowBegin = k*rowLimit;
                int rowEnd = (k+1)*rowLimit;
                if(k==world_size) rowEnd = dvaluesSize[z];
                for (int i = rowBegin; i < rowEnd; i++) {
                    int elmCountPerRow = matrixRowptr[i + 1] - matrixRowptr[i];
                    for (int j = 0; j < elmCountPerRow; j++) {
                        colInd = matrixColind[ global_OffDCount++ ];
                        if (colInd > rowEnd || colInd < rowBegin) {
                            //if(!graph[i][colInd] && !graph[colInd][i]){
                                graph[i][colInd] = 1;
                                //confSize++;
                            //}
                        }
                    }
                }
                std::cout  << "Row piece " << k << " end --------------------------------------------------------  # Conflicts: " << endl;
            }
            /*vecPtr[z] = conflicts;

            int pieceVecSize= vecPtr[z].size()/(world_size-1) + 0.5; // ceiling function
            int nwSize;
            for(int t=1; t<world_size; t++){
                vector<pair<int,int>> pieceVector;
                int upperLimit = min( t*pieceVecSize , (int) vecPtr[z].size() );
                std::copy(vecPtr[z].begin()+((t-1)*pieceVecSize), vecPtr[z].begin()+upperLimit , back_inserter(pieceVector));
                nwSize=pieceVector.size();
                // send it to child^th process
                cout << "Rank 0 sending Child " << t << " conf piece of size: " << nwSize << endl;
                MPI_Send(&nwSize, 1, MPI_INT, t, 0, MPI_COMM_WORLD);
                MPI_Send((void*) pieceVector.data(), sizeof(pair<int,int>) * nwSize, MPI_BYTE, t, 0, MPI_COMM_WORLD);
            }
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


            cout << "Result: " << pieceVector.size() <<  endl;
             */
        }
    }
    else {
        vector<double*> Local_valuesPtrs, Local_dvaluesPtrs;
        vector<int*> Local_colindPtrs, Local_rowptrPtrs;
        /*for (int i = 0; i < MATRIX_COUNT; i++) {
            MPI_Bcast(&size, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
            int *temp = new int[size];
            MPI_Bcast(temp, size, MPI_INT, 0, MPI_COMM_WORLD);
            Local_rowptrPtrs.push_back(temp);

            MPI_Bcast(&size, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
            int *temp2 = new int[size];
            MPI_Bcast(temp2, size, MPI_INT, 0, MPI_COMM_WORLD);
            Local_colindPtrs.push_back(temp2);

            MPI_Bcast(&size, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
            double* temp3 = new double[size];
            MPI_Bcast(temp3, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            Local_valuesPtrs.push_back(temp3);

            MPI_Bcast(&size, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
            double* temp4 = new double[size];
            MPI_Bcast(temp4, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            Local_dvaluesPtrs.push_back(temp4);


            std::cout <<"Rank: " << my_rank << "Completed Matrix Type: " << i << std::endl;
        }
         x
        int vecsize;
        vector<pair<int,int>> conflicts_bones10;
        MPI_Recv(&vecsize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        conflicts_bones10.resize(vecsize);
        MPI_Recv((void*)  conflicts_bones10.data(), sizeof(pair<int,int>) * vecsize, MPI_BYTE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // cleansing phase. keep unique edges as graph is undirectional.
        /*for(int z=0; z<1; z++) {
            removeDuplicateEdges(globalConflicts_bones10);
            std::cout  << "---------------------------------------------------------------------------------------------------------------- # Total conflicts globalConflicts_bones10" <<  globalConflicts_bones10.size() << endl;
        }
        removeDuplicateEdges(conflicts_bones10, vecsize);
        cout << "Rank " << my_rank << " received conf piece of size: " << conflicts_bones10.size() << " finished cleansing." << endl;

        vecsize = conflicts_bones10.size();
        MPI_Send(&vecsize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send((void*) conflicts_bones10.data(), sizeof(pair<int,int>) * vecsize, MPI_BYTE, 0, 0, MPI_COMM_WORLD);
        //cout << "Rank " << my_rank << " sent back: " << vecsize <<  endl;
        */
    }

    // Finalize MPI
    // This must always be called after all other MPI functions
    MPI_Finalize();

    return 0;
}