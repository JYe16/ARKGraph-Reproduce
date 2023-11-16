//
// Created by Jackie on 10/4/2023.
// 数据读取: vector<vector<float>>
// ifstream很可能不行

#include <vector>
#include "iostream"
#include "string"
#include "fstream"

using namespace std;

void load_data(const string& filename, int N, int dim);

int main() {
    string data_path = "data/siftsmall_base.fvecs";

    load_data(data_path, 1000, 98);

    return 0;
}

void load_data(const string& filename, int N, int dim) {
    ifstream fvecsFile;

    fvecsFile.open(filename);
    vector<vector<float>> results;

    if (!fvecsFile.is_open()) {
        cout << "File not found!" << endl;
        return;
    }

    while (!fvecsFile.eof()) {
        cout << "found" << endl;
    }

    fvecsFile.close();
}

