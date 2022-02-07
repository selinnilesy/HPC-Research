#include "header.h"
// don't use this file in SSSconflictFree.cpp as mpi compiler peerceives as double inclusion of variables in header. (due to this file being cpp file I suppose.)
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
                myfile.close();
                continue;
            }
            else if(dir_entry.path().stem() == "col") {
                int tempValInt;
                vector<int> tempVecInt;
                while (myfile >> tempValInt) {
                    tempVecInt.push_back(tempValInt);
                }
                colind_read.push_back(tempVecInt);
                myfile.close();
                continue;
            }
            while (myfile >> tempVal) {
                tempVec.push_back(tempVal);
            }
            if(dir_entry.path().stem() == "dvalues") dvalues_read.push_back(tempVec);
            else if(dir_entry.path().stem() == "nond_values") values_read.push_back(tempVec);
            else cout << "unexpected file name: " << dir_entry.path() << endl;
            tempVec.clear();
            myfile.close();
        }
    }
    return 0;
}