#include <iostream>
#include <vector>
#include  "header.h"

using namespace std;
vector<double> tempVec;

int readSSSFormat(int z) {
    fs::path matrixFolder;
    matrixFolder = "/home/selin/HPC-Research/" ;
    for(auto const& dir_entry: fs::directory_iterator{matrixFolder}){
        std::fstream myfile(dir_entry.path(), std::ios_base::in);
        if(dir_entry.path().stem() == (matrix_names[z] + "-nnzs")){
            double tempVal;
            while (myfile >> tempVal) {
                tempVec.push_back(tempVal);
            }
            cout << dir_entry.path() << " has been read with size: " <<  tempVec.size() << endl;
        }
    }
    return 0;
}

int func(int N, int P, int* row_index, int* row_size, int** alloc , int* rows_allocated, int* bucket_size){
    int size, greedy_proc, k, i, j;
    for(int i=0; i< P; i++){
        bucket_size[i] = row_size[i];
        rows_allocated[i] =1;
        alloc[i][0] = row_index[i];
    }
    for(int i=P; i< N; i++){
        size = bucket_size[0];
        greedy_proc = 0;
        for(int j=1; j< P; j++){
            if(bucket_size[j] < size)
            {   size = bucket_size[j];
                greedy_proc = j;
            }
        }
        bucket_size[greedy_proc] += row_size[i];
        k=rows_allocated[greedy_proc];
        alloc[greedy_proc][k] = row_index[i];
        rows_allocated[greedy_proc]++;
    }
    cout << "bitti" << endl << flush;
    for(int i=0; i< P; i++){
        cout << "Process: " << endl;
        for(int j=0; j< rows_allocated[i]; j++){
            if(alloc[i][j]) cout << alloc[i][j] << endl;
        }
    }
 
}

bool compareFunc(std::pair<int, int> &a, std::pair<int, int> &b)
{
    return a.first > b.first;
}

int main()
{
    cout<<"Hello World";
    int N=5;
    int P=4;
    int *row_index = new int[N];
    int* row_size = new int[N];
    int** alloc = new int*[P];
    for(int i=0; i< P; i++){
        alloc[i] = new int[N];
    }
    int* rows_allocated = new int[P];
    int* bucket_size = new int[P];

    /*
    readSSSFormat(1);

    vector<pair<int, int>> valInd;
    for(int i=0; i< tempVec.size(); i++){
        std::pair<int, int> temp;
        temp.first = tempVec[i];
        temp.second = i;
        valInd.push_back(temp);
    }

    std::sort(valInd.begin(), valInd.end(), compareFunc);
    */
    for(int i=0; i< N; i++){
        row_size[i] = 5;
        row_index[i] = i;
    }
    
    func( N, P, row_index, row_size, alloc, rows_allocated, bucket_size );
    return 0;
}
