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

int readCooFormat(int z, double ratio, bool inner) {
    cout <<  " start reading coo files..." << endl;
    double doubleVal;
    int intVal;
    const fs::path matrixFolder{"/home/selin/Split-Data/" + matrix_names[z]};
    for(auto const& dir_entry: fs::directory_iterator{matrixFolder}) {
        if (inner && dir_entry.path().stem() == "inner") {
            string rowFileName = dir_entry.path().string() + "/coordinate-" + to_string(ratio) + "-row.txt";
            string colFileName = dir_entry.path().string()+ "/coordinate-" + to_string(ratio) + "-col.txt";
            string valFileName = dir_entry.path().string() + "/coordinate-" + to_string(ratio) + "-val.txt";

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
            string rowFileName = dir_entry.path().string() + "/coordinate-" + to_string(ratio) + "-row.txt";
            string colFileName = dir_entry.path().string()+ "/coordinate-" + to_string(ratio) + "-col.txt";
            string valFileName = dir_entry.path().string() + "/coordinate-" + to_string(ratio) + "-val.txt";

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
    cout << "input ratio: " << inputRatio << endl;
    // inner read = 1 , middle read = 0 !!!
    readCooFormat(inputType, inputRatio, 1);
    readDiag(inputType);
    /*
     * [1 0 0 0 0  0 0 0    ]
     * [8 1 0 0 0  0 0 0    ]
     * [9 5 1 0 0  0 0 0    ]
     * [5 7 6 1 0  0 0 0   ]
     * [3 5 4 3 1  0 0 0   ]
     * [0 3 5 3 4  1 0 0   ]
     * [0 0 3 5 9  9 1 0   ]
     * [0 0 0 3 5  8 2 1   ]
     *

    rowVec.push_back(2);
    rowVec.push_back(3);

    rowVec.push_back(3);
    rowVec.push_back(4);
    rowVec.push_back(4);
    rowVec.push_back(4);

    rowVec.push_back(5);
    rowVec.push_back(5);
    rowVec.push_back(5);
    rowVec.push_back(5);

    rowVec.push_back(6);
    rowVec.push_back(6);
    rowVec.push_back(6);
    rowVec.push_back(6);

    rowVec.push_back(7);
    rowVec.push_back(7);
    rowVec.push_back(7);
    rowVec.push_back(7);

    rowVec.push_back(8);
    rowVec.push_back(8);
    rowVec.push_back(8);
    rowVec.push_back(8);

    // ---
    valVec.push_back(8);
    valVec.push_back(9);
    valVec.push_back(5);
    valVec.push_back(5);

    valVec.push_back(7);
    valVec.push_back(6);
    valVec.push_back(3);
    valVec.push_back(5);

    valVec.push_back(4);
    valVec.push_back(3);
    valVec.push_back(3);
    valVec.push_back(5);

    valVec.push_back(3);
    valVec.push_back(4);
    valVec.push_back(3);
    valVec.push_back(5);

    valVec.push_back(9);
    valVec.push_back(9);
    valVec.push_back(3);
    valVec.push_back(5);

    valVec.push_back(8);
    valVec.push_back(2);
 */


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
    middleBandwith = bandwithSize[inputType] - innerBandwith;
    //innerBandwith=4;
    //middleBandwith=0;
    cout << "inner bandwith: " << innerBandwith << endl;
    cout << "middle bandwith: " << middleBandwith << endl;
    //cout << "total bandwith: " << bandwithSize[inputType] << endl;

    int i,j;
    int size= n;
    //int size=8;
    // this is except diagonal.
    int k = innerBandwith;
    int lda = k+1;
    int size_1,size_2;

    size_1 = size;
    size_2=lda;
   float** A = new float*[size];
   for( i=0; i<size; i++) {
       A[i]  = new float[lda];
   }

   for( i=0; i<size; i++) {
       for( j=0; j<lda; j++) {
           A[i][j]  = 0.0;
       }
   }


    /*
   size_1 = size-innerBandwith-1;
   size_2=middleBandwith-innerBandwith;
   float** A_middle = new float*[size_1];
   for( i=0; i<size-innerBandwith-1; i++) {
       A_middle[i]  = new float[size_2];
   }

   for( i=0; i<size-innerBandwith-1; i++) {
       for( j=0; j<middleBandwith-innerBandwith; j++) {
           A_middle[i][j]  = 0.0;
       }
   }
     */


    cout << "A and A-middle initialized to 0.0" <<  endl;
    int row, col, counter=0, x, diff, neededCol, newdiff;
    double val;
    // i keeps track of whole element count.
   // cout << "start generating banded for inner" <<  endl;

   for( i=0; i<inner_rowVec.size(); i++) {
       row = inner_rowVec[i] - 1;
       col = inner_colVec[i];
       val = inner_valVec[i];
       //cout << "row: " << row << " col: " << col <<  " val: " << val << endl;
       //if(counter==innerBandwith) counter=0;

       if(row <= innerBandwith){
           // insert first element
           neededCol = 1;
           diff=col-neededCol;
           // be careful about col diff here as well !!
           A[row][((int)innerBandwith) - row + (diff)] =  val;
           //cout << "inserted: " << val << endl;
           for(x=0; x<row-1; x++){
               // do not assume all entries in the row(sub-diags) are filled.
               if( inner_rowVec[i + (x+1)]-1 != row) break;
               val = inner_valVec[i + (x+1)];
               diff =  (inner_colVec[i + (x+1)]) - neededCol;
               A[row][((int)innerBandwith) - row + (diff)] =  val;
               //cout << "before bw wrote " << val <<endl;
           }
           i+=x;
       }
       // soldan saga doldurmaya baslayabilirsin artik.
       else{
           neededCol = (row+1) - innerBandwith;
           diff=col-neededCol;
           //cout << "finished first part" << " " ;
           A[row][diff] =  val;
           //cout << "wrote " << val <<endl;
           //counter++;
       }
      // cout << i << " ";
   }
   cout << "writing also diag onto inner-A..." << endl;
   for(i=0; i<inner_diagVec.size(); i++){
       A[i][(int) innerBandwith] = inner_diagVec[i];
   }


   /*
    cout << "start generating banded for middle-A" <<  endl;
    x=counter=0;
    for( i=0; i<middle_rowVec.size(); i++) {
        row = middle_rowVec[i] - innerBandwith - 1;
        //col = colVec[i] - 1;
        val = middle_valVec[i];
        //cout << "row: " << row << " col: " << col <<  " val: " << val << endl;
        if(counter==size_2) counter=0;

        if((row+innerBandwith) <= middleBandwith){
            // insert first element
            A_middle[row-1][((int)size_2) - row] =  val;
            //cout << "inserted: " << val << endl;
            for(x=0; x<row-1; x++){
                val = middle_valVec[i + (x+1)];
                A_middle[row-1][((int)size_2) - row + (x+1)] =  val;
                //cout << "before bw wrote " << val <<endl;
            }
            i+=x;
        }
        // soldan saga doldurmaya baslayabilirsin artik.
        else{
            //cout << "finished first part" << " " ;
            A_middle[row-1][counter] =  val;
            //cout << "wrote " << val <<endl;
            counter++;
        }
        // cout << i << " ";
    }
    */

    cout << "algos finished." << endl;

    ofstream myfile;
    string output =  "/home/selin/Split-Data/"+ matrix_names[inputType]+ "/inner-banded-A" + to_string(inputRatio)+ ".txt";
    myfile.open(output, ios::out | ios::trunc);


    cout << "Formed A: " << endl;
    for( i=0; i<size_1; i++) {
        for( j=0 ; j<size_2; j++){
            myfile << A[i][j] << " " ;
        }
        myfile <<  endl;
        myfile <<  endl;
    }

    /* needed data for openblas ssbmv
   float* X = new float[size];
   for(int i=0; i<size; i++) X[i] = 1.0;
   float* Y = new float[size];
   for(int i=0; i<size; i++) Y[i] = 0.0;
   float alpha = 1;
   float beta = 0;
   int incx = 1;
   int incy = 1;
    */

    //cout << "Call cblas_ssbmv. " << endl ;

    // column major : upper verince kernel lower'a giriyor. lower verince upper'a.
    // column major : CblasLower -> ( ! kernel'de upper) 10x10 1 bandwith icin calisiyor.
    //cblas_ssbmv(CblasColMajor, CblasUpper, size, k, alpha, A, lda, X, incx, beta, Y, incy);

    //cout << "Output: " << endl;
    //for(int j=0; j<size; j++) {
    //    cout << Y[j] << endl;
    //}
}
