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
        for(int i=0; i<MATRIX_COUNT; i++){
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
        }
    }
    else {
        vector<double*> Local_valuesPtrs, Local_dvaluesPtrs;
        vector<int*> Local_colindPtrs, Local_rowptrPtrs;
        for (int i = 0; i < MATRIX_COUNT; i++) {
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
    }

    // Finalize MPI
    // This must always be called after all other MPI functions
    MPI_Finalize();

    return 0;
}