//
// Created by Jackie on 10/4/2023.
// 数据读取: vector<vector<float>>
// ifstream很可能不行

#include <vector>
#include "iostream"
#include "string"
#include "fstream"
#include <math.h>
#include <unordered_set>
#include <algorithm>
#include "limits.h"

using namespace std;

void load_ivecs_data(const string filename, vector<vector<float> > &raw_data, unsigned &num, unsigned &dim);

float calcEulideanDistance(vector<float> v1, vector<float> v2);

bool distanceComparison(const pair<int, float> &a, const pair<int, float> &b);

vector<vector<int>> buildBruteForceIndex(vector<vector<float>> raw_data, int k);

int main() {
    string data_path = "data/siftsmall_base.fvecs";
    vector<vector<float> > raw_data;
    unsigned dim, num;
    load_ivecs_data(data_path, raw_data, num, dim);
    //raw_data.resize(10);
    vector<vector<int>> bruteForceIndex = buildBruteForceIndex(raw_data, 10);

    return 0;
}


void load_ivecs_data(const string filename, vector<vector<float> > &raw_data, unsigned &num, unsigned &dim) {
    ifstream in(filename, std::ios::binary);
    if (!in.is_open()) {
        std::cout << "open file error" << std::endl;
        exit(-1);
    }
    in.read((char *) &dim, 4);


    in.seekg(0, std::ios::end);
    std::ios::pos_type ss = in.tellg();
    size_t fsize = (size_t) ss;
    num = (unsigned) (fsize / (dim + 1) / 4);
    raw_data.resize(num);
    for (unsigned i = 0; i < num; i++) raw_data[i].resize(dim);
    in.seekg(0, std::ios::beg);
    for (size_t i = 0; i < num; i++) {
        in.seekg(4, std::ios::cur);
        in.read((char *) raw_data[i].data(), dim * 4);

    }
    in.close();

}


vector<vector<int>> buildBruteForceIndex(vector<vector<float>> raw_data, int k) {
    vector<vector<int>> result;
    vector<vector<pair<int, float>>> distancePairs;


    for (int i = 0; i < raw_data.size(); i++) {
        vector<pair<int, float>> temp_pairs;
        for (int j = 0; j < raw_data.size(); j++) {
            if (i != j) {
                pair<int, float> temp_pair;
                temp_pair.first = j;
                temp_pair.second = calcEulideanDistance(raw_data[i], raw_data[j]);
                temp_pairs.push_back(temp_pair);
            }
        }
        sort(temp_pairs.begin(), temp_pairs.end(), distanceComparison);
        distancePairs.push_back(temp_pairs);
    }

    for (int x = 0; x < distancePairs.size(); x++) {
        vector<int> temp_res;
        for (int y = 0; y < k; y++) {
            temp_res.push_back(distancePairs[x][y].first);
        }
        result.push_back(temp_res);
    }
    return result;
}


bool distanceComparison(const pair<int, float> &a, const pair<int, float> &b) {
    return a.second < b.second;
}


float calcEulideanDistance(vector<float> v1, vector<float> v2) {
    float result = 0.0;
    for (int i = 0; i < v1.size(); i++) {
        result += pow((v1[i] - v2[i]), 2);
    }
    result = sqrt(result);
    return result;
}

