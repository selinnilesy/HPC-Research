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
        double time1 = MPI_Wtime();
        for(int z=0; z<1; z++) {
            bitset<914898> *graph0;
            bitset<923136> *graph1;
            bitset<952203> *graph2;
            bitset<503625> *graph3;
            bitset<1391349> *graph4;
            bitset<943695> *graph5;
            int n = dvaluesSize[z];
            if(z==0) graph0 = new bitset<914898>[914898];
            else if(z==1) graph1 = new bitset<923136>[923136];
            else if(z==2) graph2 = new bitset<952203>[952203];
            else if(z==3) graph3 = new bitset<503625>[503625];
            else if(z==4) graph4 = new bitset<1391349>[1391349];
            else if(z==5) graph5 = new bitset<943695>[943695];

            //vector<pair<int, int> > conflicts;
            matrixOffDiagonal = valuesPtrs[z];
            matrixColind = colindPtrs[z];
            matrixRowptr = rowptrPtrs[z];
            // matrix size : boneS10_DiagonalSize x boneS10_DiagonalSize
            rowLimit = (n) / world_size + 0.5; // ceiling function
            int rowInd, colInd;
            int global_OffDCount = 0;
            int confSize=0;
            std::cout << "Rank: " << my_rank << "Matrix: " << matrix_names[z] << " broadcasted rowLimit: " << rowLimit << endl;
            MPI_Bcast(&rowLimit, 1, MPI_INT, t, MPI_COMM_WORLD);

            int nwSize;

            for(int t=1; t<world_size; t++){
                int lowerLimit = rowLimit*(t+1);
                int upperLimit = min( (t+2)*rowLimit , (int) n );
                nwSize= upperLimit - lowerLimit;
                // send it to child^th process
                cout << "Rank 0 sending Child " << t << " conf piece of size: " << nwSize << endl;
                MPI_Send(&nwSize, 1, MPI_INT, t, 0, MPI_COMM_WORLD);
                MPI_Send(&n, 1, MPI_INT, t, 0, MPI_COMM_WORLD);

                if(z==0) MPI_Send((void*) graph0[lowerLimit], (914898*1 / sizeof(byte)) * nwSize, MPI_BYTE, t, 0, MPI_COMM_WORLD);
                else if(z==1) MPI_Send((void*) graph0[lowerLimit], (923136*1 / sizeof(byte)) * nwSize, MPI_BYTE, t, 0, MPI_COMM_WORLD);
                else if(z==2) MPI_Send((void*) graph0[lowerLimit], (952203*1 / sizeof(byte)) * nwSize, MPI_BYTE, t, 0, MPI_COMM_WORLD);
                else if(z==3) MPI_Send((void*) graph0[lowerLimit], (503625*1 / sizeof(byte)) * nwSize, MPI_BYTE, t, 0, MPI_COMM_WORLD);
                else if(z==4) MPI_Send((void*) graph0[lowerLimit], (1391349*1 / sizeof(byte)) * nwSize, MPI_BYTE, t, 0, MPI_COMM_WORLD);
                else if(z==5) MPI_Send((void*) graph0[lowerLimit], (943695*1 / sizeof(byte)) * nwSize, MPI_BYTE, t, 0, MPI_COMM_WORLD);
            }


            /*int rowBegin = rowLimit;
            int rowEnd = 2*rowLimit;
            // include remaining rows
            for (int i = rowBegin; i < rowEnd; i++) {
                int elmCountPerRow = matrixRowptr[i + 1] - matrixRowptr[i];
                for (int j = 0; j < elmCountPerRow; j++) {
                    colInd = matrixColind[ global_OffDCount++ ];
                    if (colInd < rowBegin || colInd > rowEnd) {
                        if(z==0)  {if(!graph0[i].test(colInd) && !graph0[colInd].test(i)) graph0[i].set(colInd);}
                        else if(z==1) {if(!graph1[i].test(colInd) && !graph1[colInd].test(i)) graph1[i].set(colInd);}
                        else if(z==2) {if(!graph2[i].test(colInd) && !graph2[colInd].test(i)) graph2[i].set(colInd);}
                        else if(z==3) {if(!graph3[i].test(colInd) && !graph3[colInd].test(i)) graph3[i].set(colInd);}
                        else if(z==4) {if(!graph4[i].test(colInd) && !graph4[colInd].test(i)) graph4[i].set(colInd);}
                        else if(z==5) {if(!graph5[i].test(colInd) && !graph5[colInd].test(i)) graph5[i].set(colInd);}
                    }
                }
            }
            double time2 = MPI_Wtime();
            // count found conflicts.
            if(z==0){
                for(int x=0; x<n; x++)
                    confSize += graph0[x].count();
            }
            if(z==1){
                for(int x=0; x<n; x++)
                    confSize += graph1[x].count();
            }
            else if(z==2){
                for(int x=0; x<n; x++)
                    confSize += graph2[x].count();
            }
            else if(z==3){
                for(int x=0; x<n; x++)
                    confSize += graph3[x].count();
            }
            else if(z==4){
                for(int x=0; x<n; x++)
                    confSize += graph4[x].count();
            }
            else if(z==5){
                for(int x=0; x<n; x++)
                    confSize += graph5[x].count();
            }
            std::cout  << "Rank: " << my_rank << " -----------------------------------------------------------  # Conflicts: " << confSize << " Elapsed Time" << time2-time1 << endl;



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

         */

        int rowLimit,n;
        MPI_Bcast(&rowLimit, 1, MPI_INT, t, MPI_COMM_WORLD); // receive row partition size
        MPI_Bcast(&n, 1, MPI_INT, t, MPI_COMM_WORLD); // receive row partition size
        int rowBegin = (my_rank+1)*rowLimit;
        int rowEnd = (my_rank+2)*rowLimit;
        if(my_rank==world_size) rowEnd = n;
        for (int i = rowBegin; i < rowEnd; i++) {
            int elmCountPerRow = matrixRowptr[i + 1] - matrixRowptr[i];
            for (int j = 0; j < elmCountPerRow; j++) {
                colInd = matrixColind[ global_OffDCount++ ];
                if (colInd < rowBegin || colInd > rowEnd) {
                    if(z==0)  {if(!graph0[i].test(colInd) && !graph0[colInd].test(i)) graph0[i].set(colInd);}
                    else if(z==1) {if(!graph1[i].test(colInd) && !graph1[colInd].test(i)) graph1[i].set(colInd);}
                    else if(z==2) {if(!graph2[i].test(colInd) && !graph2[colInd].test(i)) graph2[i].set(colInd);}
                    else if(z==3) {if(!graph3[i].test(colInd) && !graph3[colInd].test(i)) graph3[i].set(colInd);}
                    else if(z==4) {if(!graph4[i].test(colInd) && !graph4[colInd].test(i)) graph4[i].set(colInd);}
                    else if(z==5) {if(!graph5[i].test(colInd) && !graph5[colInd].test(i)) graph5[i].set(colInd);}
                }
            }
        }

        if(z==0){
            for(int x=0; x<n; x++)
                confSize += graph0[x].count();
        }
        if(z==1){
            for(int x=0; x<n; x++)
                confSize += graph1[x].count();
        }
        else if(z==2){
            for(int x=0; x<n; x++)
                confSize += graph2[x].count();
        }
        else if(z==3){
            for(int x=0; x<n; x++)
                confSize += graph3[x].count();
        }
        else if(z==4){
            for(int x=0; x<n; x++)
                confSize += graph4[x].count();
        }
        else if(z==5){
            for(int x=0; x<n; x++)
                confSize += graph5[x].count();
        }
        std::cout  << "Row piece " << k << " end --------------------------------------------------------  # Conflicts: " << confSize << endl;

        double time2 = MPI_Wtime();
        std::cout  << "Elapsed Time: " << time2-time1 << endl;

    }

    // Finalize MPI
    // This must always be called after all other MPI functions
    MPI_Finalize();

    return 0;
}