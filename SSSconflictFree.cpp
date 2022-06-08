#include <cblas.h>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <string>
#include <cmath>
#include "header.h"
#include <time.h>

using namespace std;
vector<int> inner_colVec;
vector<int> inner_rowVec;
vector<double> inner_valVec;
vector<double> inner_diagVec;

vector<int> middle_colVec;
vector<int> middle_rowVec;
vector<double> middle_valVec;
vector<double> middle_diagVec;

int readCooFormat(int z, double ratio, double middleRatio, bool inner) {
    cout <<  " start reading coo files..." << endl;
    double doubleVal;
    int intVal;
    const fs::path matrixFolder{"/home/selin/Split-Data/" + matrix_names[z]};
    for(auto const& dir_entry: fs::directory_iterator{matrixFolder}) {
        if (inner && dir_entry.path().stem() == "inner") {
            string rowFileName = dir_entry.path().string() + "/coordinate-" + to_string(ratio)+ "-"+ to_string(middleRatio) + "-row.txt";
            string colFileName = dir_entry.path().string()+ "/coordinate-" + to_string(ratio)+ "-" + to_string(middleRatio)+ "-col.txt";
            string valFileName = dir_entry.path().string() + "/coordinate-" + to_string(ratio) + "-" + to_string(middleRatio)+  "-val.txt";

            std::fstream myfile(rowFileName, std::ios_base::in);
            while (myfile >> intVal) {
                inner_rowVec.push_back(intVal);
            }
            myfile.close();
            cout << rowFileName << " has been read with size: " << inner_rowVec.size() << endl;

            myfile.open(colFileName, std::ios_base::in);
            while (myfile >> intVal) {
                inner_colVec.push_back(intVal);
            }
            myfile.close();
            cout << colFileName << " has been read with size: " << inner_colVec.size() << endl;

            myfile.open(valFileName, std::ios_base::in);
            while (myfile >> doubleVal) {
                inner_valVec.push_back(doubleVal);
            }
            myfile.close();
            cout << valFileName << " has been read with size: " << inner_valVec.size() << endl;
        }
        else if (!inner && dir_entry.path().stem() == "middle") {
            string rowFileName = dir_entry.path().string() + "/coordinate-" + to_string(ratio)+ "-" + to_string(middleRatio) + "-row.txt";
            string colFileName = dir_entry.path().string()+ "/coordinate-" + to_string(ratio)+ "-" + to_string(middleRatio) + "-col.txt";
            string valFileName = dir_entry.path().string() + "/coordinate-" + to_string(ratio)+ "-" + to_string(middleRatio) + "-val.txt";

            std::fstream myfile(rowFileName, std::ios_base::in);
            // else, start reading doubles.
            while (myfile >> intVal) {
                middle_rowVec.push_back(intVal);
            }
            myfile.close();
            cout << rowFileName << " has been read with size: " << middle_rowVec.size() << endl;

            myfile.open(colFileName, std::ios_base::in);
            while (myfile >> intVal) {
                middle_colVec.push_back(intVal);
            }
            myfile.close();
            cout << colFileName << " has been read with size: " << middle_colVec.size() << endl;

            myfile.open(valFileName, std::ios_base::in);
            while (myfile >> doubleVal) {
                middle_valVec.push_back(doubleVal);
            }
            myfile.close();
            cout << valFileName << " has been read with size: " << middle_valVec.size() << endl;
        }
    }
    return 0;
}
int readDiag(int z) {
    cout <<  " start reading diag file..." << endl;
    double doubleVal;
    string diagFile = "/home/selin/SSS-Data/" + matrix_names[z] + "/diag.txt" ;

    std::fstream myfile(diagFile, std::ios_base::in);
    // else, start reading doubles.
    while (myfile >> doubleVal) {
        inner_diagVec.push_back(doubleVal);
    }
    myfile.close();
    cout << diagFile << " has been read." << endl;
    return 0;
}

int main(int argc, char **argv)
{

    int n = matrixSize[atoi(argv[1])];
    int inputType = atoi(argv[1]);
    double inputRatio = atof(argv[2]);
    double middleRatio = atof(argv[3]);
    bool inner = atoi(argv[4]);
    cout << "input ratio: " << inputRatio << endl;
    cout << "middle ratio: " << middleRatio << endl;
    cout << "inner: " << inner << endl;
    // inner read = 1 , middle read = 0 !!!
    readCooFormat(inputType, inputRatio, middleRatio, inner);
    readDiag(inputType);


    if(inner_rowVec.size() != inner_valVec.size()){
        cout << "not equal size - rowVec and valVec (NNZ) !!! " << endl;
        return -1;
    }
    if(n != inner_diagVec.size()){
        cout << "not equal size - diagVec (n) !!! " << endl;
        return -1;
    }


    double innerBandwith, middleBandwith;
    innerBandwith = (int) (nnz_n_Ratios[inputType]*bandwithProportions[inputType] * inputRatio);
    middleBandwith = (int) ((bandwithSize[inputType] - innerBandwith)*middleRatio);
    //innerBandwith=4;
    //middleBandwith=0;
    cout << "inner bandwith: " << innerBandwith << endl;
    cout << "middle bandwith: " << middleBandwith << endl;
    //cout << "total bandwith: " << bandwithSize[inputType] << endl;

    int i,j;
    int size= n;
    int size_1,size_2;
    int k,lda;
    double* A, **A_middle;
    int row, col, diff, neededCol, x;
    double val;
    int kl, ku;

    if(inner) {
        k = innerBandwith;
        kl=ku=k;
        //lda = k+1;
        lda =kl+ku+1;
        //size= n+kl;
        size_1=n+kl;
        size_2 = lda+1;
        A = new double[size_1*size_2];
        /*
        for( i=0; i<size_1; i++) {
            A[i]  = new double[size_2];
        }
         */
        for( i=0; i<size_1; i++) {
            for( j=0; j<size_2; j++) {
                A[i*size_2 + j]  = 0.0;
            }
        }
    }
    else if(!inner){
        k = middleBandwith;
        // BE CAREFUL WITH THIS !!!!!!!!!!!!!!!!!
        lda = k+1;
        size_1 = size-innerBandwith;
        size_2=middleBandwith;
        A_middle = new double*[size_1];
        for( i=0; i<size_1; i++) {
            A_middle[i]  = new double[size_2];
        }

        for( i=0; i<size_1; i++) {
            for( j=0; j<size_2; j++) {
                A_middle[i][j]  = 0.0;
            }
        }
    }
    cout << "A/A-middle initialized to 0.0" <<  endl;

    // i keeps track of whole element count.
    // cout << "start generating banded for inner" <<  endl;

    int rowDiff;
    for( i=0; i<inner_rowVec.size(); i++) {
        row = inner_rowVec[i] - 1;
        col = inner_colVec[i];
        val = inner_valVec[i];
        if(row <= innerBandwith){
            neededCol = 1;
            diff=col-neededCol;
            A[row*size_2 + ((int)innerBandwith) - row + (diff)] =  -val;
            // for dbmv
            rowDiff=  row-diff ;
            A[(row-rowDiff) * size_2 + (lda-1)  -  (((int)innerBandwith) - row + (diff))   ] =  val;
            for(x=0; x<row-1; x++){
                if( inner_rowVec[i + (x+1)]-1 != row) break;
                val = inner_valVec[i + (x+1)];
                diff =  (inner_colVec[i + (x+1)]) - neededCol;
                rowDiff=  row-diff ;
                A[row* size_2  +  ((int)innerBandwith) - row + (diff)] =  -val;
                // for dbmv
                A[(row-rowDiff)*size_2 + (lda-1)  -  (((int)innerBandwith) - row + (diff))   ] =  val;
            }
            i+=x;
        }
        else{
            neededCol = (row+1) - innerBandwith;
            diff=col-neededCol;
            A[row * size_2 + diff] =  -val;
            // for dbmv
            A[(row-((int)innerBandwith-diff)) * size_2 + (lda-1)  - diff] =  val;
        }
    }
    cout << "writing also diag onto inner-A..." << endl;
    for(i=0; i<inner_diagVec.size(); i++){
        A[i*size_2 + (int) innerBandwith] = inner_diagVec[i];
    }
    cout << "completed band storage." << endl;

    middle_colVec.clear();
    middle_rowVec.clear();
    middle_valVec.clear();
    middle_diagVec.clear();

    /*
    cout << "start generating banded for middle-A" <<  endl;
    int row, col, x, diff, neededCol;
    double val;
    x=0;
    for( i=0; i<middle_rowVec.size(); i++) {
        row = middle_rowVec[i] - innerBandwith - 1;
        col = middle_colVec[i];
        val = middle_valVec[i];
        //cout << "row: " << row << " col: " << col <<  " val: " << val << endl;
        if(row <= middleBandwith){
            neededCol = 1;
            diff=col-neededCol;
            // insert first element
            A_middle[row][((int)size_2) - row + (diff)] =  val;
            //cout << "inserted: " << val << endl;
            for(x=0; x<row-1; x++){
                if( middle_rowVec[i + (x+1)]-innerBandwith-1 != row) break;
                val = middle_valVec[i + (x+1)];
                diff =  (middle_colVec[i + (x+1)]) - neededCol;
                A_middle[row][((int)size_2) - row + diff] =  val;
            }
            i+=x;
        }
        else{
            neededCol = row - middleBandwith + 1;
            diff=col-neededCol;
            A_middle[row][diff] =  val;
        }
    }
    cout << "algos finished. Formed A/A_middle." << endl;
     */


    ofstream myfile;
    string output;
    /*

       if(inner) output = "/home/selin/Split-Data/" + matrix_names[inputType] + "/inner/dgbmv-inner-banded-A"  + to_string(inputRatio) + ".txt";
       if(!inner) output =  "/home/selin/Split-Data/" + matrix_names[inputType] + "/middle-banded-A"  + to_string(inputRatio) + "-" + to_string(middleRatio) + ".txt";
       myfile.open(output, ios::out | ios::trunc);
       cout << "writing A/A_middle ..." << endl;
       for( i=0; i<size_1; i++) {
           for( j=0 ; j<size_2; j++){
               myfile << A[i*size_2 + j] << " " ;
           }
           myfile <<  endl;
           myfile <<  endl;
       }
       cout << "written A/A_middle." << endl;
       myfile.close();
       */



    double* X = new double[n];
    for(i=0; i<n; i++) X[i] = 1.0;
    double* Y = new double[n];
    for(i=0; i<n; i++) Y[i] = 0.0;

    double alpha = 1.0;
    double beta = 0.0;
    int incx = 1;
    int incy = 1;


    clock_t t;
    cout << "Call cblas_dgbmv... " << endl ;
    t = clock();
    //for(int i=0; i<0000; i++)
    //cblas_dsbmv(CblasColMajor, CblasUpper, n, k, alpha, B, lda, X, incx, beta, Y, incy);

    // for dgbmv, A is one dimensional anyway.
    cblas_dgbmv(CblasColMajor, CblasNoTrans , n, n, kl, ku, alpha, A, lda, X, incx, beta, Y, incy);
    t = clock() - t;
    printf ("It took me %f seconds.\n",((float)t)/CLOCKS_PER_SEC);

    if(inner) output = "/home/selin/Outputs/" + matrix_names[inputType] + "/dgbmv-inner-"  + to_string(inputRatio) + ".txt";
    if(!inner) output =  "/home/selin/Outputs/" + matrix_names[inputType] + "/middle-"  + to_string(inputRatio) + "-" + to_string(middleRatio) + ".txt";
    myfile.open(output, ios::out | ios::trunc);

    cout << "Writing Y: " << endl;
    for( i=0; i<n; i++) {
        myfile << Y[i] << " " ;
    }
    myfile.close();
    cout << "Output completed." << output << endl;
     


    if(X) {
        delete [] X;
        cout << "deleted X." << endl;
    }
    if(Y) {
        delete [] Y;
        cout << "deleted Y." << endl;
    }

    if(inner) delete [] A;
    else delete [] A_middle;
}