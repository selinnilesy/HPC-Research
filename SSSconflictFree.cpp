//
// Created by Selin Yıldırım on 7.01.2022.
//

#include <iostream>
#include <omp.h>
#include  "header.h"
#include <limits>
#include <iomanip>
//#include <mpi.h>
//#include "rcmtest.cpp"
//#include "geeks.cpp"

using namespace std;
#define MATRIX_COUNT 6

typedef std::numeric_limits< double > dbl;

int readSSSFormat(int z, bool banded) {
    fs::path matrixFolder;
    if(!banded) matrixFolder = "/home/selin/SSS-Data/" + matrix_names[z] + "/unbanded" ;
    else matrixFolder = "/home/selin/SSS-Data/" + matrix_names[z];
    for(auto const& dir_entry: fs::directory_iterator{matrixFolder}){
        std::fstream myfile(dir_entry.path(), std::ios_base::in);
        if(dir_entry.path().stem() == "rowptr") {
            int tempValInt;
            vector<int> tempVecInt;
            while (myfile >> tempValInt) {
                tempVecInt.push_back(tempValInt);
            }
            rowptrPtrs.push_back(new int[tempVecInt.size()]);
            int *temp = rowptrPtrs[0];
            for(int i=0; i<tempVecInt.size(); i++) temp[i]=tempVecInt[i];
            rowptrSize.push_back(tempVecInt.size());
            cout << dir_entry.path() << " has been read with size: " <<tempVecInt.size() << endl;
            myfile.close();
        }
        else if(dir_entry.path().stem() == "colind") {
            int tempValInt;
            vector<int> tempVecInt;
            while (myfile >> tempValInt) {
                tempVecInt.push_back(tempValInt);
            }
            colindPtrs.push_back(new int[tempVecInt.size()]);
            int *temp = colindPtrs[0];
            for(int i=0; i<tempVecInt.size(); i++) temp[i]=tempVecInt[i];
            colindSize.push_back(tempVecInt.size());
            cout << dir_entry.path() << " has been read with size: " <<  tempVecInt.size() << endl;
            myfile.close();
        }
        else if(dir_entry.path().stem() == "diag"){
            double tempVal;
            vector<double> tempVec;
            while (myfile >> tempVal) {
                tempVec.push_back(tempVal);
            }
            dvaluesPtrs.push_back(new double[tempVec.size()]);
            double *temp = dvaluesPtrs[0];
            for(int i=0; i<tempVec.size(); i++) temp[i]=tempVec[i];
            dvaluesSize.push_back(tempVec.size());
            cout << dir_entry.path() << " has been read with size: " <<  tempVec.size() << endl;
            myfile.close();
        }
        else if(dir_entry.path().stem() == "vals"){
            double tempVal;
            vector<double> tempVec;
            while (myfile >> tempVal) {
                tempVec.push_back(tempVal);
            }
            valuesPtrs.push_back(new double[tempVec.size()]);
            double *temp = valuesPtrs[0];
            for(int i=0; i<tempVec.size(); i++) temp[i]=tempVec[i];
            valuesSize.push_back(tempVec.size());
            cout << dir_entry.path() << " has been read with size: " <<  tempVec.size() << endl;
            myfile.close();
        }
        else cout << "unexpected file name: " << dir_entry.path() << endl;
    }
    return 0;
}
/*
    A RNxN : matrix in SSS format
    x RN : input vector
    y RN : output vector

 */

int main(int argc, char **argv){
    int n;
    //init();
    if(!argv[1]){
        cout << "please provide input matrix index (int): boneS10, Emilia_923, ldoor, af_5_k101, Serena, audikw_1" << endl;
        return -1;
    }
    if(!argv[2]){
        cout << "please provide bool for banded" << endl;
        return -1;
    }
    if(!argv[3]){
        cout << "please state if its openmp or serial" << endl;
        return -1;
    }
    bool banded = atoi(argv[2]);
    bool parallel = atoi(argv[3]);
    
    readSSSFormat(atoi(argv[1]) , banded);

    n = matrixSize[atoi(argv[1])];
    int inputType = atoi(argv[1]);


    double *matrixOffDiagonal = valuesPtrs[0];
    double *matrixDiagonal = dvaluesPtrs[0];
    int *matrixColind = colindPtrs[0];
    int *matrixRowptr= rowptrPtrs[0];
    double *x = new double[n];
    double *y = new double[n];
    for(int i=0; i<n; i++) x[i] = 1.0;
    for(int i=0; i<n; i++) y[i] = 0.0;
    ofstream myfile1;

    int colInd;
    if(!parallel){
        cout << "start computing serial SSS mv..." << endl;
        double row_i,row_e, val, itime, ftime;
        itime = omp_get_wtime();
        for (int run = 0; run < 1000; run++) {
            for (int i = 0; i < n; i++) {
                val=0.0;
                row_i = matrixRowptr[i] - 1;
                row_e = matrixRowptr[i+1] - 1;
                val += matrixDiagonal[i] * x[i];
                for (int j =row_i; j < row_e; j++) {
                    colInd = matrixColind[j] - 1;
                    val -= matrixOffDiagonal[j] * x[colInd];
                    y[colInd] += matrixOffDiagonal[j] * x[i];
                }
                y[i] += val;
            }
        }
        ftime=omp_get_wtime();
        printf ("It took me %f seconds for 1000-times serial run.\n", ftime-itime);

        if(!banded) myfile1.open ("/home/selin/Seq-Results/" + matrix_names[inputType] + "/unbanded/result.txt", ios::out | ios::trunc);
        else myfile1.open ("/home/selin/Seq-Results/" + matrix_names[inputType] + "/banded/result.txt", ios::out | ios::trunc);
    }
    else{
        int threadCount = atoi(argv[4]);
        if(!banded) myfile1.open ("/home/selin/Seq-Results/" + matrix_names[inputType] + "/unbanded/result_omp.txt", ios::out | ios::trunc);
        else myfile1.open ("/home/selin/Seq-Results/" + matrix_names[inputType] + "/banded/result_omp.txt", ios::out | ios::trunc);
        int  perThread=0, row_elm, conflictSize, col_elm, pieceSize, istart, imax, jmax;
        pieceSize = (n/threadCount);
        cout << "start computing parallel SSS mv..." << endl;
        double itime, ftime, val, row_i,row_e;
        int *rowStart = new int[threadCount+1];
        for(int z=0; z<threadCount; z++){
            rowStart[z] = z*pieceSize;
        }
        rowStart[threadCount] = n;
        for(int z=1; z<threadCount; z++){
            conflictSize += rowStart[z];
        }
        double *conflicts = new double[(threadCount-1)*n];
        for(int z=0; z<(threadCount-1)*n; z++){
            conflicts[z] = 0.0;
        }
        //for (int tx = 0; tx < 7; tx++) {
        //threadCount *= 2;
        cout << "Threads: " << threadCount << endl;
        //block = n/threadCount;
        double x_time = omp_get_wtime();
        omp_set_dynamic(0);     // Explicitly disable dynamic teams
        omp_set_num_threads(threadCount); // Use 4 threads for all consecutive parallel regions
        for (int run = 0; run < 1000; run++) {
            //itime = omp_get_wtime();
            #pragma omp parallel for private(val, colInd, row_i, row_e)
            for (int i = 0; i < n; i++) {
                //perThread++;
                row_i = matrixRowptr[i] - 1;
                row_e = matrixRowptr[i + 1] - 1;
                val = matrixDiagonal[i] * x[i];
                for (int j = row_i; j < row_e; j++) {
                    colInd = matrixColind[j] - 1;
                    // lower part is going to be negated for DKEW SYMMETRIC.
                    val -= matrixOffDiagonal[j] * x[colInd];
                    // upper part's sign will be kept
                    if(colInd < rowStart[omp_get_thread_num()]){
                        conflicts[((omp_get_thread_num()-1)*n) + colInd] += matrixOffDiagonal[j] * x[i];
                    }
                        //#pragma omp atomic update
                    else y[colInd] += matrixOffDiagonal[j] * x[i];
                }
                //#pragma omp atomic update
                y[i] += val;
            }
            #pragma omp barrier
            omp_set_num_threads(threadCount-1); // Use 3 threads for below
            #pragma omp parallel for
            for (int px = 0; px < threadCount-1; px++) {
                for (int nx = 0; nx < n; nx++) {
                    y[nx] += conflicts[(omp_get_thread_num()*n) + nx];
                    //cout << "y[" << nx << "]: " << y[nx] << " after conf: " << conflicts[(omp_get_thread_num()*n) + nx] << " \t\t";
                }
            }
            //ftime = omp_get_wtime();
            //myfile1 << ftime - itime << " ";
        }
        ftime = omp_get_wtime();
        printf ("It took me %f seconds for 1000-times serial run.\n", ftime-x_time);
        // }

        myfile1.close();
        cout << "Completed output... " << endl;
        }

    /*
    cout << "Writing to output... " << endl;
    for (int i=0; i<n; i++) {
        myfile1 << std::fixed << std::setprecision(dbl::max_digits10) << y[i] << '\t';
    }
    myfile1.close();
    cout << "Completed output... " << endl;
    */



    delete [] x;
    delete [] y;
    delete [] matrixRowptr;
    delete [] matrixColind;
    delete [] matrixOffDiagonal;
    delete [] matrixDiagonal;
    return 0;
}