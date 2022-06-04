#include <stdio.h>
#include <iostream>
#include <string.h>
#include <string>
#include <cmath>
#include "header.h"

using namespace std;
vector<double> dgbmv_res;
vector<double> dsbmv_res;

int readSSBMVResult(int z, double ratio) {
    cout <<  " start readSSBMVResult..." << endl;
    double doubleVal;
    int intVal;
    const fs::path matrixFolder{"/home/selin/Outputs/" + matrix_names[z]};
    for(auto const& dir_entry: fs::directory_iterator{matrixFolder}) {
        if (dir_entry.path().stem() == ("inner-" + to_string(ratio))) {

            fstream myfile;
            myfile.open(dir_entry.path().string(), std::ios_base::in);

            while (myfile >> doubleVal) {
                dsbmv_res.push_back(doubleVal);
            }
            myfile.close();
            cout << dir_entry.path().string() << " has been read with size: " << dsbmv_res.size() << endl;
        }
    }
    return 0;
}
int readDGBMVResult(int z, double ratio) {
    cout <<  " start readDGBMVResult..." << endl;
    double doubleVal;
    int intVal;
    const fs::path matrixFolder{"/home/selin/Outputs/" + matrix_names[z]};
    for(auto const& dir_entry: fs::directory_iterator{matrixFolder}) {
        if ( dir_entry.path().stem() == ("dgbmv-inner-" + to_string(ratio))) {

            fstream myfile;
            myfile.open(dir_entry.path().string(), std::ios_base::in);

            while (myfile >> doubleVal) {
                dgbmv_res.push_back(doubleVal);
            }
            myfile.close();
            cout << dir_entry.path().string() << " has been read with size: " << dgbmv_res.size() << endl;
        }
    }
    return 0;
}

int main(int argc, char **argv)
{

    int n = matrixSize[atoi(argv[1])];
    int inputType = atoi(argv[1]);
    double inputRatio = atof(argv[2]);

    cout << "input ratio: " << inputRatio << endl;
    // inner read = 1 , middle read = 0 !!!
    readSSBMVResult(inputType, inputRatio);
    readDGBMVResult(inputType, inputRatio);


    if(dgbmv_res.size() != dsbmv_res.size()){
        cout << "not equal size !!! " << endl;
        return -1;
    }


    for(int i=0; i<dgbmv_res.size(); i++){
        cout << "not equal index: " << i << ". vals: " << dgbmv_res[i] << "," << dsbmv_res[i] << endl;
        return -1;
    }
    cout << "all equal. " << endl;
    return 0;
}