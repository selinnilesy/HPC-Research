#include <stdio.h>
#include <iostream>
#include <string.h>
#include <string>
#include <cmath>
#include "header.h"

using namespace std;

vector<double> val;
vector<int> row,col;
int *rowPtr, *colPtr;
double *valPtr;
int nnz;

extern "C" {
    extern void amux_(int *n, double* x, double *y, double *a, int *ja, int *ia);
}

int readMiddleCSRStorage(int z, double ratio, double middleR, bool lower) {
    cout <<  " start reading Upper CSR Storage..." << endl;

    double doubleVal;
    int intVal;
    string fileName, rowfile, colfile, valfile;
    int counter=0;
    // middle read not yet IMPLEMENTED !!!
    // TO DO
    if(!lower) fileName = "/home/selin/Split-Data/" + matrix_names[z]  + "/middle/CSR-Data/upper";
    else fileName = "/home/selin/Split-Data/" + matrix_names[z] + "/middle/CSR-Data" ;

    rowfile = fileName + "/" + to_string(ratio) + "-" +to_string(middleR) + "-row.txt";
    colfile = fileName + "/" + to_string(ratio) + "-" +to_string(middleR) + "-col.txt";
    valfile = fileName + "/" + to_string(ratio) + "-" +to_string(middleR) + "-val.txt";

    std::fstream myfile(rowfile, std::ios_base::in);
    while (myfile >> intVal) {
        row.push_back(intVal);
    }
    myfile.close();
    cout << rowfile << " has been read with size: " << row.size() << endl;

    myfile.open(colfile, std::ios_base::in);
    while (myfile >> intVal) {
        col.push_back(intVal);
    }
    myfile.close();

    nnz=col.size();
    cout << colfile << " has been read with size: " << col.size() << endl;

    myfile.open(valfile, std::ios_base::in);
    while (myfile >> doubleVal) {
        val.push_back(doubleVal);
    }
    myfile.close();
    cout << valfile << " has been read with size: " << val.size() << endl;
    return 0;
}

int readCSRStorage(int z, double ratio, bool lower) {
    cout <<  " start reading Upper CSR Storage..." << endl;

    double doubleVal;
    int intVal;
    string fileName, rowfile, colfile, valfile;
    int counter=0;
    // middle read not yet IMPLEMENTED !!!
    // TO DO
    if(!lower) fileName = "/home/selin/Split-Data/" + matrix_names[z]  + "/inner-outer-equal/inner/CSR-Data/upper";
    else fileName = "/home/selin/Split-Data/" + matrix_names[z] + "/inner-outer-equal/inner/CSR-Data" ;

    rowfile = fileName + "/" + to_string(ratio) + "-row.txt";
    colfile = fileName + "/" + to_string(ratio) + "-col.txt";
    valfile = fileName + "/" + to_string(ratio) + "-val.txt";

    std::fstream myfile(rowfile, std::ios_base::in);
    while (myfile >> intVal) {
        row.push_back(intVal);
    }
    myfile.close();
    cout << rowfile << " has been read with size: " << row.size() << endl;

    myfile.open(colfile, std::ios_base::in);
    while (myfile >> intVal) {
        col.push_back(intVal);
    }
    myfile.close();

    nnz=col.size();
    cout << colfile << " has been read with size: " << col.size() << endl;

    myfile.open(valfile, std::ios_base::in);
    while (myfile >> doubleVal) {
        val.push_back(doubleVal);
    }
    myfile.close();
    cout << valfile << " has been read with size: " << val.size() << endl;
    return 0;
}


int main(int argc, char **argv)
{
    int n = matrixSize[atoi(argv[1])];
    int inputType = atoi(argv[1]);
    double inputRatio = atof(argv[2]);
    double middleRatio = atof(argv[3]);

    bool lower = atoi(argv[4]);
    bool inner = atoi(argv[5]);
    cout << "input ratio: " << inputRatio << endl;
    cout << "middle ratio: " << middleRatio << endl;
    cout << "lower bool: " << lower << endl;

    if(inner) readCSRStorage( inputType, inputRatio, lower);
    else readMiddleCSRStorage( inputType, inputRatio,middleRatio, lower);

    rowPtr = new int[matrixSize[inputType] + 1];
    colPtr = new int[nnz];
    valPtr = new double[nnz];


    for(int i=0; i<row.size(); i++) rowPtr[i] = row[i];
    for(int i=0; i<col.size(); i++) colPtr[i] = col[i];
    for(int i=0; i<val.size(); i++) valPtr[i] = val[i];

    cout << "test: " << rowPtr[100] << endl;
    cout << "test: " << colPtr[100] << endl;
    cout << "test: " << valPtr[100] << endl;
    int nrow= matrixSize[inputType];

    double* X = new double[nrow];
    for(int i=0; i<nrow; i++) X[i] = 1.0;

    double* Y = new double[nrow];
    for(int i=0; i<nrow; i++) Y[i] = 0.0;

    cout << "starts computing amux..." << endl;
    amux_(&nrow, X, Y, valPtr, colPtr, rowPtr);
    std::cout  <<  " finished computing amux_... " << endl;


    ofstream myfile;
    string output;
    if(inner){
        if(lower) output = "/home/selin/Split-Data/" + matrix_names[inputType]  + "/inner-outer-equal/inner/CSR-Data/" + to_string(inputRatio) + "-result.txt";
        else if(!lower) output = "/home/selin/Split-Data/" + matrix_names[inputType]  + "/inner-outer-equal/inner/CSR-Data/upper/" + to_string(inputRatio) + "-result.txt";
    }
    else{
        if(lower) output = "/home/selin/Split-Data/" + matrix_names[inputType]  + "/middle/CSR-Data/" + to_string(inputRatio) + "-" +to_string(middleRatio) + "-result.txt";
        else if(!lower) output = "/home/selin/Split-Data/" + matrix_names[inputType]  + "/middle/CSR-Data/upper/" + to_string(inputRatio) + "-" +to_string(middleRatio)  + "-result.txt";
    }

    myfile.open(output, ios::out | ios::trunc);

    cout << "Writing result: " << endl;
    for(int i=0; i<nrow; i++) {
        myfile << Y[i] << " " ;
    }
    myfile.close();
    cout << "Output completed." << endl;

    delete [] X;
    delete [] Y;
    delete [] rowPtr;
    delete [] colPtr;
    delete [] valPtr;

}