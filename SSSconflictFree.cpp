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

int readSSSFormat() {
    double tempVal;
    vector<double> tempVec;
    for (int i=0; i<matrix_names.size(); i++) {
        const fs::path matrixFolder{"/home/selin/Paper-Implementation/CSR-Data/" + matrix_names[i]};
        for(auto const& dir_entry: fs::directory_iterator{matrixFolder}){
            std::fstream myfile(dir_entry.path(), std::ios_base::in);
            if(dir_entry.path().stem() == "rowptr") {
                int tempValInt;
                vector<int> tempVecInt;
                while (myfile >> tempValInt) {
                    tempVecInt.push_back(tempValInt);
                }
                rowptr_read.push_back(tempVecInt);
                cout << dir_entry.path() << " has been read." << endl;
                myfile.close();
                continue;
            }
            /*
            else if(dir_entry.path().stem() == "col") {
                int tempValInt;
                vector<int> tempVecInt;
                while (myfile >> tempValInt) {
                    tempVecInt.push_back(tempValInt);
                }
                colind_read.push_back(tempVecInt);
                cout << dir_entry.path() << " has been read." << endl;
                myfile.close();
                continue;
            }
            while (myfile >> tempVal) {
                tempVec.push_back(tempVal);
            }
            if(dir_entry.path().stem() == "dvalues") dvalues_read.push_back(tempVec);
            else if(dir_entry.path().stem() == "nond_values") values_read.push_back(tempVec);
            else cout << "unexpected file name: " << dir_entry.path() << endl;
            */
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

    vector<vector<double>> values_read_local, dvalues_read_local;
    vector<vector<int>> rowptr_read_local, colind_read_local;
    unsigned long long size;
    if(my_rank==0) {
        cout << "i call readSSSFormat. " << endl;
        readSSSFormat();
        for(int i=0; i<MATRIX_COUNT; i++){
            size =rowptr_read[i].size();
            cout << size << endl;
            MPI_Bcast(&size, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
            MPI_Bcast(&rowptr_read[i], size, MPI_INT, 0, MPI_COMM_WORLD);
        }
    }
    else {
        for(int i=0; i<MATRIX_COUNT; i++){
            MPI_Bcast(&size, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
            vector<int> temp;
            temp.resize(size);
            MPI_Bcast(&temp.front(), size, MPI_INT, 0, MPI_COMM_WORLD);
            rowptr_read_local.push_back(temp);

        }

        std::cout << "World Size: " << world_size << "   Rank: " << my_rank << std::endl;
        std::cout << "rowptr_read_local Size: " << rowptr_read_local.size() << std::endl;
        //std::cout << "colind_read Size: " << colind_read.size() << std::endl;
        //std::cout << "values_read Size: " << values_read.size() << std::endl;
        //std::cout << "dvalues_read Size: " << dvalues_read.size() << std::endl;
    }

    // Finalize MPI
    // This must always be called after all other MPI functions
    MPI_Finalize();

    return 0;
}