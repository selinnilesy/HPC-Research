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
#include <unordered_set>
#include <algorithm>    // std::sort, std::stable_sort
#include <experimental/filesystem>


using namespace std;
namespace fs = std::experimental::filesystem;

vector<string> matrix_names = {"boneS10", "Emilia_923", "ldoor", "af_5_k101", "Serena", "audikw_1"};


// Bandwiths :
// Emilia =  14672, Serena = 87872, af_5_k101 = 1274
// audikw_1 = 35102 , boneS10 = 13727, ldoor = 8707
// serena > audikw1 > emilia > bone_s10 > ldoor = ad_5_k_101

// ORDER: boneS10, Emilia_923, ldoor, af_5_k101, Serena, audikw_1
// nnz/n = 44, 43 , 44 , 34, 46, 82
vector<int> matrixSize{914898, 923136, 952203 , 503625, 1391349, 943695};
vector<int> nonzerosSize{40878708, 40373538 , 42493817 , 17550675, 64131971, 77651847};
vector<int> bandwithSize{13727, 14672 , 8707 , 1274, 87872, 35102};
vector<int> bandwithProportions{137, 146 , 87 , 12, 878, 351};
vector<int> nnz_n_Ratios{44, 43 , 44 , 34, 46, 82};

// matrix SSS storage vectors
vector<double*> valuesPtrs, dvaluesPtrs;
vector<int*> colindPtrs, rowptrPtrs;
vector<int> valuesSize, dvaluesSize, colindSize, rowptrSize;

#endif //PAPER_IMPLEMENTATION_HEADER_H
