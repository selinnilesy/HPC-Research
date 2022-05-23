#include <stdio.h>
#include <iostream>
#include <string.h>
#include <string>
#include <cmath>
#include "header.h"

using namespace std;

vector<double> lowerRes,upperRes, totalRes,diag,  sum;
int nnz;

int readResult(int z, double ratio, double mratio, int innerB) {
    cout <<  " start reading multip result lower/upper..." << endl;

    double doubleVal;
    string fileName;
    // middle read not yet IMPLEMENTED !!!
    // TO DO
    fileName = "/home/selin/Split-Data/" + matrix_names[z]  + "/middle/CSR-Data/upper/" + to_string(ratio)+ "-" + to_string(mratio) + "-result.txt";
    std::fstream myfile(fileName, std::ios_base::in);
    while (myfile >> doubleVal) {
        upperRes.push_back(doubleVal);
    }
    myfile.close();
    cout << fileName << " has been read with size: " << upperRes.size() << endl;

    fileName = "/home/selin/Split-Data/" + matrix_names[z]  + "/middle/CSR-Data/" + to_string(ratio)+ "-" + to_string(mratio)+ "-result.txt";
    myfile.open(fileName, std::ios_base::in);
    while (myfile >> doubleVal) {
        lowerRes.push_back(doubleVal);
    }
    myfile.close();
    nnz=lowerRes.size();
    cout << fileName << " has been read with size: " << lowerRes.size() << endl;

    fileName = "/home/selin/SSS-Data/" + matrix_names[z]  + "/diag.txt";
    myfile.open(fileName, std::ios_base::in);
    while (myfile >> doubleVal) {
        diag.push_back(doubleVal);
    }
    myfile.close();
    cout << fileName << " has been read with size: " << diag.size() << endl;

    fileName = "/home/selin/Outputs/" + matrix_names[z]  + "/middle-" + to_string(ratio)+ "-" + to_string(mratio) + ".txt";
    myfile.open(fileName, std::ios_base::in);
    while (myfile >> doubleVal) {
        totalRes.push_back(doubleVal);
    }
    myfile.close();
    cout << fileName << " has been read with size: " << totalRes.size() << endl;

    return 0;
}


int main(int argc, char **argv)
{
    int n = matrixSize[atoi(argv[1])];
    int inputType = atoi(argv[1]);
    double inputRatio = atof(argv[2]);
    double middleRatio = atof(argv[3]);
    cout << "input ratio: " << inputRatio << endl;
    cout << "middle ratio: " << middleRatio << endl;

    int innerBandwith = (int) (nnz_n_Ratios[inputType]*bandwithProportions[inputType] * inputRatio);
    int middleBandwith = (int) ((bandwithSize[inputType] - innerBandwith)*middleRatio);
    readResult(inputType, inputRatio, middleRatio, innerBandwith);
    cout << "lowerRes: " << lowerRes[0] << endl;
    cout << "upperRes: " << upperRes[0] << endl;
    cout << "diag: " << diag[0] << endl;
    cout << "cblas_dsbmv: " << totalRes[0] << endl;


    for(int i=0; i<lowerRes.size(); i++){
        sum.push_back(lowerRes[i] + upperRes[i] + diag[i]);
    }
    for(int i=0; i<sum.size(); i++){
        if(abs(sum[i] - totalRes[i]) > 0.01) {
            cout << "not equal index: " << i << " you found: " <<sum[i] << " banded computed: " << totalRes[i];
            return -1;
        }
    }
    return 0;

}