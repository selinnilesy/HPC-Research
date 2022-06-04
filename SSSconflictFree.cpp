#include <stdio.h>
#include <iostream>
#include <string.h>
#include <string>
#include <cmath>
#include "header.h"

using namespace std;
vector<double> serialBanded_res;
vector<double> dsbmv_res;

int readDSBMVResult(int z, double ratio) {
    cout <<  " start readDSBMVResult..." << endl;
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
           //cout << dsbmv_res[dsbmv_res.size()-3] << " " << dsbmv_res[dsbmv_res.size()-2] << " "<< dsbmv_res[dsbmv_res.size()-1];

            myfile.close();
            cout << dir_entry.path().string() << " has been read with size: " << dsbmv_res.size() << endl;
        }
    }
    return 0;
}
int readSerialBandedResult(int z, double ratio) {
    cout <<  " start readSerialBandedResult..." << endl;
    double doubleVal;
    int intVal;
    const fs::path matrixFolder{"/home/selin/Seq-Results/" + matrix_names[z] +"/banded"};
    for(auto const& dir_entry: fs::directory_iterator{matrixFolder}) {
        if ( dir_entry.path().stem() == "result") {

            fstream myfile;
            myfile.open(dir_entry.path().string(), std::ios_base::in);

            while (myfile >> doubleVal) {
                serialBanded_res.push_back(doubleVal);
            }
            myfile.close();
            cout << dir_entry.path().string() << " has been read with size: " << serialBanded_res.size() << endl;
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
    readDSBMVResult(inputType, inputRatio);
    readSerialBandedResult(inputType, inputRatio);


    if(serialBanded_res.size() != dsbmv_res.size()){
        cout << "not equal size !!! " << endl;
        return -1;
    }


    for(int i=0; i<dsbmv_res.size(); i++){
        if(serialBanded_res[i] != dsbmv_res[i]){
            cout << "not equal index: " << i << ". vals: " << serialBanded_res[i] << "," << dsbmv_res[i] << endl;
            return -1;
        }
    }
    cout << "all equal. " << endl;
    return 0;
}