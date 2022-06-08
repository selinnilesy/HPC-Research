#include <stdio.h>
#include <iostream>
#include <string.h>
#include <string>
#include <cmath>
#include "header.h"

using namespace std;

vector<double> lowerRes,upperRes, dsbmvRes,diag, dgbmvRes, serialRes;
vector<double> inner, inner_outer;
int nnz;

int readResult(int z, double ratio) {
    cout <<  " start reading multip result lower/upper..." << endl;

    double doubleVal;
    string fileName;
    // middle read not yet IMPLEMENTED !!!
    // TO DO
    fileName = "/home/selin/Split-Data/" + matrix_names[z]  + "/inner/CSR-Data/upper/" + to_string(ratio)+ "-result.txt";
    std::fstream myfile(fileName, std::ios_base::in);
    while (myfile >> doubleVal) {
        upperRes.push_back(doubleVal);
    }
    myfile.close();
    cout << fileName << " has been read with size: " << upperRes.size() << endl;

    fileName = "/home/selin/Split-Data/" + matrix_names[z]  + "/inner/CSR-Data/" + to_string(ratio)+ "-result.txt";
    myfile.open(fileName, std::ios_base::in);
    while (myfile >> doubleVal) {
        lowerRes.push_back(doubleVal);
    }
    myfile.close();
    nnz=lowerRes.size();
    cout << fileName << " has been read with size: " << lowerRes.size() << endl;

    fileName = "/home/selin/Outputs/" + matrix_names[z]  + "/inner-" + to_string(ratio)+ ".txt";
    myfile.open(fileName, std::ios_base::in);
    while (myfile >> doubleVal) {
        dsbmvRes.push_back(doubleVal);
    }
    myfile.close();
    cout << fileName << " has been read with size: " << dsbmvRes.size() << endl;

    fileName = "/home/selin/Outputs/" + matrix_names[z]  + "/dgbmv-inner-" + to_string(ratio)+ ".txt";
    myfile.open(fileName, std::ios_base::in);
    while (myfile >> doubleVal) {
        dgbmvRes.push_back(doubleVal);
    }
    myfile.close();
    cout << fileName << " has been read with size: " << dgbmvRes.size() << endl;

    fileName = "/home/selin/SSS-Data/" + matrix_names[z]  + "/diag.txt";
    myfile.open(fileName, std::ios_base::in);
    while (myfile >> doubleVal) {
        diag.push_back(doubleVal);
    }
    myfile.close();
    cout << fileName << " has been read with size: " << diag.size() << endl;

    fileName = "/home/selin/Split-Data/" + matrix_names[z]  + "/inner/CSR-Data/serial-result.txt";
    myfile.open(fileName, std::ios_base::in);
    while (myfile >> doubleVal) {
        serialRes.push_back(doubleVal);
    }
    myfile.close();
    cout << fileName << " has been read with size: " << serialRes.size() << endl;

    /*
    fileName = "/home/selin/Split-Data/" + matrix_names[z]  + "/inner/coordinate-" +to_string(ratio) + "-0.300000-col.txt";
    myfile.open(fileName, std::ios_base::in);
    while (myfile >> doubleVal) {
        inner.push_back(doubleVal);
    }
    myfile.close();
    cout << fileName << " has been read with size: " << inner.size() << endl;

    fileName = "/home/selin/Split-Data/" + matrix_names[z]  + "/inner-outer-equal/inner/coordinate-" +to_string(ratio) + "-col.txt";
    myfile.open(fileName, std::ios_base::in);
    while (myfile >> doubleVal) {
        inner_outer.push_back(doubleVal);
    }
    myfile.close();
    cout << fileName << " has been read with size: " << inner_outer.size() << endl;
    for(int i=0; i<inner_outer.size(); i++){
        if(inner_outer[i] != inner[i]) {
            cout << "not equal inner-outer and inner data" << endl;
            break;
        }
    }
    inner_outer.clear();
    inner.clear();

    fileName = "/home/selin/Split-Data/" + matrix_names[z]  + "/inner/coordinate-" +to_string(ratio) + "-0.300000-row.txt";
    myfile.open(fileName, std::ios_base::in);
    while (myfile >> doubleVal) {
        inner.push_back(doubleVal);
    }
    myfile.close();
    cout << fileName << " has been read with size: " << inner.size() << endl;

    fileName = "/home/selin/Split-Data/" + matrix_names[z]  + "/inner-outer-equal/inner/coordinate-" +to_string(ratio) + "-row.txt";
    myfile.open(fileName, std::ios_base::in);
    while (myfile >> doubleVal) {
        inner_outer.push_back(doubleVal);
    }
    myfile.close();
    cout << fileName << " has been read with size: " << inner_outer.size() << endl;
    for(int i=0; i<inner_outer.size(); i++){
        if(inner_outer[i] != inner[i]) {
            cout << "not equal inner-outer and inner data" << endl;
            break;
        }
    }
    inner_outer.clear();
    inner.clear();

    fileName = "/home/selin/Split-Data/" + matrix_names[z]  + "/inner/coordinate-" +to_string(ratio) + "-0.300000-val.txt";
    myfile.open(fileName, std::ios_base::in);
    while (myfile >> doubleVal) {
        inner.push_back(doubleVal);
    }
    myfile.close();
    cout << fileName << " has been read with size: " << inner.size() << endl;

    fileName = "/home/selin/Split-Data/" + matrix_names[z]  + "/inner-outer-equal/inner/coordinate-" +to_string(ratio) + "-val.txt";
    myfile.open(fileName, std::ios_base::in);
    while (myfile >> doubleVal) {
        inner_outer.push_back(doubleVal);
    }
    myfile.close();
    cout << fileName << " has been read with size: " << inner_outer.size() << endl;
    for(int i=0; i<inner_outer.size(); i++){
        if(inner_outer[i] != inner[i]) {
            cout << "not equal inner-outer and inner data" << endl;
            break;
        }
    }
    inner_outer.clear();
    inner.clear();
     */
    return 0;
}


int main(int argc, char **argv)
{
    int n = matrixSize[atoi(argv[1])];
    int inputType = atoi(argv[1]);
    double inputRatio = atof(argv[2]);
    cout << "input ratio: " << inputRatio << endl;

    readResult(inputType, inputRatio);

    cout << endl;
    cout << "-----------checking DSBMV by lower/upper AMUX: " << endl;
    std::cout.precision(10);
    for(int i=0; i<dsbmvRes.size(); i++){
        if(abs( (lowerRes[i] + upperRes[i] + diag[i] )- dsbmvRes[i]) > 0.1 ) {
            cout << "not equal - index: " << i << " correct result: " <<  (lowerRes[i] + upperRes[i] + diag[i] ) << " dsbmv computed: " << dsbmvRes[i] << " with difference: " << abs( (lowerRes[i] + upperRes[i] + diag[i] )- dsbmvRes[i]) << endl;
            break;
        }
    }
    cout << "-----------checking DSBMV by Serial Multiplication: " << endl;
    for(int i=0; i<dsbmvRes.size(); i++){
        if(abs( serialRes[i]+ diag[i]- dsbmvRes[i]) > 0.1 ) {
            cout << "not equal - index: " << i << " correct result: " <<  serialRes[i]+ diag[i] << " dsbmv computed: " << dsbmvRes[i] << " with difference: " << abs( serialRes[i]+ diag[i]- dsbmvRes[i]) << endl;
            break;
        }
    }
    cout << "-----------checking DGBMV: " << endl;
    for(int i=0; i<dgbmvRes.size(); i++){
        if(abs( (lowerRes[i] + upperRes[i] + diag[i] )- dgbmvRes[i]) > 0.1) {
            cout << "not equal - index: " << i << " correct result: " <<  (lowerRes[i] + upperRes[i] + diag[i] ) << " dgbmv computed: " << dgbmvRes[i] << " with difference: " << abs( (lowerRes[i] + upperRes[i] + diag[i] )- dgbmvRes[i]) << endl;
            break;
        }
    }
    return 0;

}