#include <cblas.h>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <string>
#include <cmath>
#include "header.h"

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
            // else, start reading doubles.
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
    //int size=8;
    // this is except diagonal.
    int size_1,size_2;
    int k,lda;
    double** A, ** A_middle;
    int row, col, diff, neededCol, x;
    double val;

    if(inner) {
        k = innerBandwith;
        lda = k+1;
        size= n;
        size_1 = size;
        size_2=lda;
        A = new double*[size];
        for( i=0; i<size; i++) {
            A[i]  = new double[lda];
        }
        for( i=0; i<size; i++) {
            for( j=0; j<lda; j++) {
                A[i][j]  = 0.0;
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

   for( i=0; i<inner_rowVec.size(); i++) {
       row = inner_rowVec[i] - 1;
       col = inner_colVec[i];
       val = inner_valVec[i];
       //cout << "row: " << row << " col: " << col <<  " val: " << val << endl;
       if(row <= innerBandwith){
           // insert found row-element
           neededCol = 1;
           diff=col-neededCol;
           A[row][((int)innerBandwith) - row + (diff)] =  val;
           //cout << "inserted: " << val << endl;
           for(x=0; x<row-1; x++){
               // do not assume all entries in the row(sub-diags) are filled.
               if( inner_rowVec[i + (x+1)]-1 != row) break;
               val = inner_valVec[i + (x+1)];
               diff =  (inner_colVec[i + (x+1)]) - neededCol;
               A[row][((int)innerBandwith) - row + (diff)] =  val;
           }
           i+=x;
       }
       // soldan saga doldurmaya baslayabilirsin artik.
       else{
           neededCol = (row+1) - innerBandwith;
           diff=col-neededCol;
           A[row][diff] =  val;
           //cout << "wrote " << val <<endl;
       }
      // cout << i << " ";
   }
   cout << "writing also diag onto inner-A..." << endl;
   for(i=0; i<inner_diagVec.size(); i++){
       A[i][(int) innerBandwith] = inner_diagVec[i];
   }



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
    if(inner) output = "/home/selin/Split-Data/" + matrix_names[inputType] + "/inner-outer-equal/inner/inner-banded-A2"  + to_string(inputRatio) + ".txt";
    if(!inner) output =  "/home/selin/Split-Data/" + matrix_names[inputType] + "/middle-banded-A"  + to_string(inputRatio) + "-" + to_string(middleRatio) + ".txt";
    myfile.open(output, ios::out | ios::trunc);
    cout << "writing A/A_middle ..." << endl;
    for( i=0; i<size_1; i++) {
        for( j=0 ; j<size_2; j++){
            myfile << A[i][j] << " " ;
        }
        myfile <<  endl;
        myfile <<  endl;
    }
    cout << "written A/A_middle." << endl;
    myfile.close();
     */


    double* X = new double[n];
    for(int i=0; i<n; i++) X[i] = 1.0;
    double* Y = new double[n];
    for(int i=0; i<n; i++) Y[i] = 0.0;
    double alpha = 1;
    double beta = 0;
    int incx = 1;
    int incy = 1;

    vector<double> copy_banded;
    for(int i=0; i<size; i++) {
        for(int j=0; j<lda; j++) {
            copy_banded.push_back(A[i][j]);
        }
    }
    if(inner) delete [] A;
    else delete [] A_middle;

    double *B = new double[size*lda];
    for(int i=0; i<size*lda; i++) B[i] = copy_banded[i];

    /*
    if(inner){
        B = *A;
    }
    else if(!inner){
        B = *A_middle;
    }
     */

    cout << "Call cblas_ssbmv. " << endl ;
    // BE CAREFUL WITH K=LDA CASE WHEN USING MIDDLE = !INNER
    cblas_dsbmv(CblasColMajor, CblasUpper, n, k, alpha, B, lda, X, incx, beta, Y, incy);

    if(inner) output = "/home/selin/Outputs/" + matrix_names[inputType] + "/inner-"  + to_string(inputRatio) + ".txt";
    if(!inner) output =  "/home/selin/Outputs/" + matrix_names[inputType] + "/middle-"  + to_string(inputRatio) + "-" + to_string(middleRatio) + ".txt";
    myfile.open(output, ios::out | ios::trunc);

    cout << "Writing Y: " << endl;
    for( i=0; i<n; i++) {
        myfile << Y[i] << " " ;
    }
    myfile.close();
    cout << "Output completed." << output << endl;

    delete [] X;
    delete [] Y;
    //if(inner) delete [] A;
    //else delete [] A_middle;
}