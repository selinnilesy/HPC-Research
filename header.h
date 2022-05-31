//
// Created by Selin Yıldırım on 8.01.2022.
//

#ifndef PAPER_IMPLEMENTATION_HEADER_H
#define PAPER_IMPLEMENTATION_HEADER_H

#include <iostream>
#include <vector>
#include <numeric>      // iota
#include <fstream>
#include <cstdlib>
#include <iterator>
#include <bitset>
#include <cmath>
#include <cstring>
#include <unordered_set>
#include <algorithm>    // std::sort, std::stable_sort
#include <experimental/filesystem>
#include <cstdint> // include this header for uint64_t


using namespace std;
namespace fs = std::experimental::filesystem;

vector<string> matrix_names = {"boneS10", "Emilia_923", "ldoor", "af_5_k101", "Serena", "audikw_1"};
// ORDER: boneS10, Emilia_923, ldoor, af_5_k101, Serena, audikw_1
vector<int> matrixSize{914898, 923136, 952203 , 503625, 1391349, 943695};
vector<int> nonzerosSize{40878708, 40373538 , 42493817 , 17550675, 64131971, 77651847};

// matrix SSS storage vectors
vector<double*> valuesPtrs, dvaluesPtrs;
vector<int*> colindPtrs, rowptrPtrs;
vector<int> valuesSize, dvaluesSize, colindSize, rowptrSize;
double *values_boneS10, *values_Emilia_923, *values_ldoor, *values_af_5_k101, *values_audikw_1;
double *dvalues_boneS10, *dvalues_Emilia_923, *dvalues_ldoor, *dvalues_af_5_k101, *dvalues_audikw_1;
int *rowptr_boneS10, *rowptr_Emilia_923, *rowptr_ldoor, *rowptr_af_5_k101, *rowptr_audikw_1;
int *colind_boneS10, *colind_Emilia_923, *colind_ldoor, *colind_af_5_k101, *colind_audikw_1;

#endif //PAPER_IMPLEMENTATION_HEADER_H
