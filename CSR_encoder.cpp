//
// Created by Selin Yıldırım on 23.12.2021.
//
#include "header.h"

vector<int> row_ptr;
// diag values
vector<double> dvalues;
// coordinates
vector<int> unordered_rows, rows;
vector<int> unordered_cols, cols;
// non-diag values
vector<double> unordered_values, values;

vector<int> sort_indexes(const vector<int> &v) {

    vector<int> idx(v.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    stable_sort(idx.begin(), idx.end(),
                [&v]( int i1,  int i2) {return v[i1] < v[i2];});

    return idx;
}

char matrixTypeRecognizer(string m){
    if(m=="boneS10" ) return 0;
    else if(m=="Emilia_923") return 1;
    else if(m=="ldoor") return 2;
    else if(m=="af_5_k101") return 3;
    else if(m=="Serena") return 4;
    else if(m=="audikw_1") return 5;
    else{ cout << "DID NOT RECOGNIZE MATRIX TYPE!"; return 0;}
}

void printToFiles(string m, const vector<double> &dvalues_i, const vector<double> &values_i, const vector<int> &rowptr_i, const vector<int> &cols_i, char matrixTypeIndex){
    std::ofstream dval_file, nondval_file, col_file, rowptr_file;
    dval_file.open("./CSR-Data/" + m + "/dvalues.txt", ios::out | ios::trunc );
    nondval_file.open("./CSR-Data/" + m + "/nond_values.txt", ios::out | ios::trunc );
    col_file.open("./CSR-Data/" + m + "/col.txt", ios::out | ios::trunc );
    rowptr_file.open("./CSR-Data/" + m + "/rowptr.txt", ios::out | ios::trunc );

    if (dval_file.is_open()) {
        for (auto val : dvalues_i) {
            dval_file << val << "\t";
            dval_file.flush();
            if (dval_file.bad()) {
                throw std::runtime_error("Failed to write to outfile!");
            }

        }
        cout << "Writing dvalues of " + m + " is complete." << endl;
        dval_file.close();
    }

   if (nondval_file.is_open()) {
        for (auto val : values_i) {
            nondval_file << val << "\t";
            nondval_file.flush();
            if (nondval_file.bad()) {
                throw std::runtime_error("Failed to write to outfile!");
            }

        }
        cout << "Writing non-diagonal values of " + m + " is complete." << endl;
        nondval_file.close();
    }

    if (col_file.is_open()) {
        for (auto col_x : cols_i) {
            col_file << col_x << "\t";
            col_file.flush();
            if (col_file.bad()) {
                throw std::runtime_error("Failed to write to outfile!");
            }

        }
        cout << "Writing columns of " + m + " is complete." << endl;
        col_file.close();
    }

    if (rowptr_file.is_open()) {
        for (auto val : rowptr_i) {
            rowptr_file << val << "\t";
            rowptr_file.flush();
            if (rowptr_file.bad()) {
                throw std::runtime_error("Failed to write to outfile!");
            }

        }
        cout << "Writing rowptr of " + m + " is complete." << endl;
        rowptr_file.close();
    }
}

int main() {
   std::string path = "/home/selin/Paper-Implementation/matrices/";
   std::string ext(".mat");
   double num = 0.0;
   char matrixTypeIndex;

   /* not needed anymore, hardcoded.
   for (auto &p : fs::recursive_directory_iterator(path))
   {
       if (p.path().extension() == ext)
           matrix_names.push_back(p.path().stem().string());
   }
    */

   for (auto &m : matrix_names) {
       cout << m << endl;
       matrixTypeIndex = matrixTypeRecognizer(m);

       std::ifstream diag_f(path + m + "-diag.txt", std::ios::in);
       std::ifstream row_f(path + m + "-row_except_diag.txt", std::ios::in);
       std::ifstream col_f(path + m + "-col_except_diag.txt", std::ios::in);
       std::ifstream nondiag_f(path + m + "-values_except_diag.txt", std::ios::in);
       std::ifstream rowptr_f(path + m + "-rowptr.txt", std::ios::in);

       //check to see that the file was opened correctly:
       if (!diag_f.is_open() || !row_f.is_open() || !col_f.is_open() || !nondiag_f.is_open() || !rowptr_f.is_open()) {
           std::cerr << "There was a problem opening the input file!\n";
           exit(1);//exit or do additional error checking
       }

       //keep storing values from the text file so long as data exists:
       while (row_f >> num) {
           unordered_rows.push_back(num);
       }
       cout << "reading file of unordered rows complete." << endl;while (diag_f >> num) {
           dvalues.push_back(num);
       }
       cout << "reading file of diag values complete." << endl;
       while (col_f >> num) {
           unordered_cols.push_back(num);
       }
       cout << "reading file of unordered col indices complete." << endl;
       while (nondiag_f >> num) {
           unordered_values.push_back(num);
       }
       cout << "reading file of nondiag values complete." << endl;
       while (rowptr_f >> num) {
           row_ptr.push_back(num);
       }
       cout << "reading file of rowptr complete." << endl;
       diag_f.close();
       row_f.close();
       col_f.close();
       nondiag_f.close();
       rowptr_f.close();
       // conform to the order of CSR encoding
       // by reordering MATLAB's nonzero find function's output.
       if(unordered_cols.size() != unordered_rows.size() || unordered_values.size()!=unordered_rows.size() || unordered_values.size()!=unordered_cols.size()){
           cout << "INEQUAL !! " << endl;
           continue;
       }
       for (auto i: sort_indexes(unordered_rows)) {
           // CAREFUL ! MATLAB indexes were 1-based.
           rows.push_back(unordered_rows[i]);
           cols.push_back(unordered_cols[i]);
           values.push_back(unordered_values[i]);
       }

       // clean unordered vectors, not used anymore.
       unordered_cols.clear();
       unordered_rows.clear();
       unordered_values.clear();

       cout << "started printing..." << endl;
       printToFiles(m, dvalues, values, row_ptr, cols, matrixTypeIndex);


        // cleaning for the next matrix
        dvalues.clear();
        cols.clear();
        values.clear();
        rows.clear();
        row_ptr.clear();

    }
    return 0;
}