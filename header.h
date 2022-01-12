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
#include <algorithm>    // std::sort, std::stable_sort
#include <experimental/filesystem>

namespace fs = std::experimental::filesystem;
using namespace std;
using namespace fs;

vector<string> matrix_names = {"boneS10", "Emilia_923", "ldoor", "af_5_k101", "Serena", "audikw_1"};
// ORDER: boneS10, Emilia_923, ldoor, af_5_k101, Serena, audikw_1
vector<int> matrixSize{914898, 923136, 952203 , 503625, 1391349, 943695};

// matrix SSS storage vectors
vector<vector<double>> values_read, dvalues_read;
vector<vector<int>> rowptr_read, colind_read;

#endif //PAPER_IMPLEMENTATION_HEADER_H
