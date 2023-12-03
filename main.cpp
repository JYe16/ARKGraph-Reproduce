//
// Created by Jackie on 10/4/2023.

#include <vector>
#include "iostream"
#include "string"
#include "fstream"
#include <cmath>
#include <algorithm>
#include <assert.h>

using namespace std;

struct compactGraph {
    int bl = 0;
    int br = 0;
    int el = 0;
    int er = 0;
    vector<int> C;
};

struct naiveGraph {
    int b = 0;
    int e = 0;
    vector<int> C;
};

struct partialRangeGraph {
    int bl;
    int br;
    vector<int> C;
};

struct deltaCompressionGraph {
    int lMin;
    int j;
    vector<int> C;
};

void loadData(const string &filename, vector<vector<float> > &raw_data, unsigned &num, unsigned &dim);

float calcEuclideanDistance(vector<float> v1, vector<float> v2);

bool distanceComparison(const pair<int, float> &a, const pair<int, float> &b);

vector<vector<naiveGraph>> buildBruteForceIndex(vector<vector<float>> raw_data, int k);

vector<vector<compactGraph>> buildCompactGraph(vector<vector<float>> raw_data, int k);

vector<pair<int, float>> calcAndSortDistance(vector<vector<float>> raw_data, vector<float> vi, int index);

vector<vector<float>> cutData(vector<vector<float>> originalData, int x, int y);

long countGraphSize(vector<vector<compactGraph>> G);

vector<int> mergeLR(vector<int> L, vector<int> R, vector<int> LR);

long countBruteForceSize(vector<vector<naiveGraph>> G);

void buildPartialRangeGraph(vector<vector<float>> raw_data, int k, vector<vector<partialRangeGraph>> &Gb,
                            vector<vector<partialRangeGraph>> &Ge);

vector<int> buildKNNPartialRange(vector<vector<float>> rawData, vector<vector<partialRangeGraph>> Gb,
                                 vector<vector<partialRangeGraph>> Ge, int x, int y, int i, int k);

vector<int> buildKNNCompactList(vector<vector<compactGraph>> G, int x, int y, int i);

vector<int> buildKNNBruteForce(vector<vector<naiveGraph>> G, int x, int y, int i);

vector<pair<int, float>>
calcAndSortDistanceBF(vector<vector<float>> raw_data, vector<float> vi, int index, int x, int y);

void buildDeltaCompressionIndex(vector<vector<float>> raw_data, int k, vector<vector<deltaCompressionGraph>> &dGb,
                                vector<vector<deltaCompressionGraph>> &dGe);

vector<int> buildKNNDeltaCompression(vector<vector<float>> rawData, vector<vector<deltaCompressionGraph>> dGb,
                                     vector<vector<deltaCompressionGraph>> dGe, int x, int y, int i, int k);

/*
 *
 * The following two functions are borrowed from the official implementation of the paper
 * GitHub Link: https://github.com/rutgers-db/ARKGraph
 * These functions will be used to read BigGraph Data
 *
 */

void ReadMatFromTsv(const string &path, vector<vector<float>> &data, const int length_limit);

void Split(std::string &s, std::string &delim, std::vector<std::string> *ret);

int main() {
    char dataset;
    int k, N;
    string data_path;
    vector<vector<float> > raw_data;
    unsigned dim, num;

    cout << "Program Started... Please enter the dataset you would like to index: 'b' for BigGraph, 'd' for DEEP10M: " << endl;
    cin >> dataset;
    cout << "Please enter the value of N (total number of data): " << endl;
    cin >> N;
    cout << "Please enter the value of k (number of nearest neighbors): " << endl;
    cin >> k;

    if (dataset == 'b') {
        data_path = "data/biggraph.tsv";
        ReadMatFromTsv(data_path, raw_data, N);
    } else if (dataset == 'd') {
        data_path = "data/deep10M.fvecs";
        loadData(data_path, raw_data, num, dim);
    } else {
        cout << "Invalid argument: dataset. Please try again." << endl;
        return -1;
    }

    if (N > 200) {
        char cont;
        cout << "Warning: N is greater than 200. The indexing time of brute force approach could be prohibitively long..." << endl;
        cout << "Continue? [Y/N]: " << endl;

        cin >> cont;

        if (cont == 'N') {
            cout << "Program Aborted..." << endl;
            return 1;
        }
    }

    raw_data.resize(N);

    cout << "Program Begin. Data size: " << raw_data.size() << " K: " << k << endl;
    vector<vector<naiveGraph>> bruteForceIndex = buildBruteForceIndex(raw_data, k);
    vector<vector<compactGraph>> compactIndex = buildCompactGraph(raw_data, k);

    vector<vector<partialRangeGraph>> Gb;
    vector<vector<partialRangeGraph>> Ge;

    buildPartialRangeGraph(raw_data, k, Gb, Ge);
    int partialRangeNumber = 0;

    for (int i = 0; i < Gb.size(); i++) {
        partialRangeNumber += Gb[i].size();
    }

    for (int i = 0; i < Ge.size(); i++) {
        partialRangeNumber += Ge[i].size();
    }

    vector<vector<deltaCompressionGraph>> dGb;
    vector<vector<deltaCompressionGraph>> dGe;

    buildDeltaCompressionIndex(raw_data, k, dGb, dGe);
    int deltaCompressionNumber = 0;

    for (int i = 0; i < Gb.size(); i++) {
        deltaCompressionNumber += dGb[i].size();
    }

    for (int i = 0; i < Ge.size(); i++) {
        deltaCompressionNumber += dGe[i].size();
    }

    cout << "Program Finished... Calculating Stats..." << endl;

    int bruteForceNumber = countBruteForceSize(bruteForceIndex);
    long compactGraphNumber = countGraphSize(compactIndex);

    cout << "Brute Force Index Number: " << bruteForceNumber << "\t\tTotal Size: "
         << (bruteForceNumber * (8 + 4 * k)) / (1024.0 * 1024.0) << " MiB" << endl;
    cout << "Compact Graph Index Number: " << compactGraphNumber << "\t\tTotal Size: "
         << (compactGraphNumber * (16 + 4 * k)) / (1024.0 * 1024.0) << " MiB" << endl;
    cout << "Partial Range Index Number: " << partialRangeNumber << "\t\tTotal Size: "
         << (partialRangeNumber * (8 + 4 * k)) / (1024.0 * 1024.0) << " MiB" << endl;
    cout << "Delta Compression Index Number: " << deltaCompressionNumber << "\t\tTotal Size: "
         << (deltaCompressionNumber * 8 / (1024.0 * 1024.0)) << " MiB" << endl;

    int indexForKNN;
    int rangeStart;
    int rangeEnd;
    cout << "Please enter the vector you would like to generate KNN (from 0 to " << N - 1 << "): " << endl;
    cin >> indexForKNN;
    cout << "Please enter the starting index you would like to generate KNN (from 0 to " << N - 1 << "): " << endl;
    cin >> rangeStart;
    cout << "Please enter the ending index you would like to generate KNN (from 0 to " << N - 1 << "): " << endl;
    cin >> rangeEnd;

    if (rangeEnd - rangeStart < k) {
        cout << "Error: range is less than k... Program Aborted..." << endl;
        return -1;
    }

    vector<int> bfknn = buildKNNBruteForce(bruteForceIndex, rangeStart, rangeEnd, indexForKNN);
    vector<int> cgknn = buildKNNCompactList(compactIndex, rangeStart, rangeEnd, indexForKNN);
    vector<int> prknn = buildKNNPartialRange(raw_data, Gb, Ge, rangeStart, rangeEnd, indexForKNN, k);
    vector<int> dcknn = buildKNNDeltaCompression(raw_data, dGb, dGe, rangeStart, rangeEnd, indexForKNN, k);

    sort(bfknn.begin(), bfknn.end());
    sort(cgknn.begin(), cgknn.end());
    sort(prknn.begin(), prknn.end());
    sort(dcknn.begin(), dcknn.end());

    cout << "BF\tCG\tPR\tDC" << endl;
    for (int i = 0; i < k; i++) {
        cout << bfknn[i] << "\t" << cgknn[i] << "\t" << prknn[i] << "\t" << dcknn[i] << endl;
    }

//    cout << "CG\tPR" << endl;
//    for (int i = 0; i < k; i ++) {
//        cout << cgknn[i] << "\t" << prknn[i] << endl;
//    }
    return 0;
}


void ReadMatFromTsv(const string &path, vector<vector<float>> &data, const int length_limit = -1) {
    ifstream infile;
    string bline;
    string delim = "\t";
    int numCols = 0;
    infile.open(path, ios::in);
    getline(infile, bline, '\n');
    if (getline(infile, bline, '\n')) {
        vector<string> ret;
        Split(bline, delim, &ret);
        numCols = ret.size();
    }
    infile.close();

    int counter = 0;
    if (length_limit == -1) counter = -9999999;
    infile.open(path, ios::in);
    // skip the first line
    getline(infile, bline, '\n');
    while (getline(infile, bline, '\n')) {
        if (counter >= length_limit) break;
        counter++;

        vector<string> ret;
        Split(bline, delim, &ret);
        vector<float> arow(numCols - 1);
        assert(ret.size() == numCols);
        for (int i = 0; i < ret.size() - 1; i++) {
            arow[i] = static_cast<float>(stod(ret[i + 1]));
        }
        data.emplace_back(arow);
    }
    infile.close();
}

void Split(std::string &s, std::string &delim, std::vector<std::string> *ret) {
    size_t last = 0;
    size_t index = s.find_first_of(delim, last);
    while (index != std::string::npos) {
        ret->push_back(s.substr(last, index - last));
        last = index + 1;
        index = s.find_first_of(delim, last);
    }
    if (index - last > 0) {
        ret->push_back(s.substr(last, index - last));
    }
}

int YT8M2Int(const string id) {
    int res = 0;
    for (size_t i = 0; i < 4; i++) {
        res *= 100;
        res += (int) id[i] - 38;
    }
    return res;
}


void loadData(const string &filename, vector<vector<float> > &raw_data, unsigned &num, unsigned &dim) {
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


vector<vector<compactGraph>> buildCompactGraph(vector<vector<float>> raw_data, int k) {

    cout << "Indexing using Compact Graph Approach..." << endl;
    int indicator = raw_data.size() / 10;
    vector<vector<compactGraph>> G;
    int window;

    for (int i = 0; i < raw_data.size(); i++) {
        if (i % indicator == 0) {
            cout << "#";
        }

        vector<pair<int, float>> sortedDistancePairs = calcAndSortDistance(raw_data, raw_data[i], i);
        // declare and initialize L and R
        vector<int> L;
        vector<int> R;
        vector<int> LR;
        vector<compactGraph> tempG;

        int prevLMin = 0;
        int LMin = 0;
        int prevRMax = raw_data.size() + 1;
        int RMax = raw_data.size() + 1;

        for (int j = 0; j < sortedDistancePairs.size(); j++) {
            int j_val = sortedDistancePairs[j].first;
            if (j_val < i) {
                if (j_val > LMin) {
                    // insert j into L
                    L.push_back(j_val);
                    sort(L.begin(), L.end());
                    if (L.size() == k + 1) {
                        prevLMin = L[0];
                        L.erase(L.begin());
                        LMin = L[0];
                    }
                    // merge L and R
                    if (L.size() + R.size() >= k && L.size() != 0 && R.size() != 0) {
                        LR = mergeLR(L, R, LR);
                        if (LR.size() == k) {
                            compactGraph graph = *new compactGraph();

                            graph.bl = prevLMin;
                            graph.br = LR[0];
                            graph.el = LR[k - 1];
                            graph.er = prevRMax;
                            for (int z = 0; z < k; z++) {
                                graph.C.push_back(LR[z]);
                            }
                            tempG.push_back(graph);
                        } else {
                            if (L.size() < k) {
                                compactGraph graph = *new compactGraph();

                                graph.bl = prevLMin;
                                graph.br = LR[0];
                                graph.el = LR[k - 1];
                                graph.er = LR[k];
                                for (int z = 0; z < k; z++) {
                                    graph.C.push_back(LR[z]);
                                }
                                tempG.push_back(graph);

                                window = 1;

                                while (window + k < LR.size()) {
                                    compactGraph graph = *new compactGraph();

                                    graph.bl = LR[window - 1];
                                    graph.br = LR[window];
                                    graph.el = LR[window + k - 1];
                                    graph.er = LR[window + k];
                                    for (int z = 0; z < k; z++) {
                                        graph.C.push_back(LR[window + z]);
                                    }
                                    tempG.push_back(graph);
                                    window++;
                                }

                                graph = *new compactGraph();

                                graph.bl = LR[window - 1];
                                graph.br = LR[window];
                                graph.el = LR[window + k - 1];
                                graph.er = prevRMax;

                                for (int z = 0; z < k; z++) {
                                    graph.C.push_back(LR[window + z]);
                                }
                                tempG.push_back(graph);

                            } else if (L.size() == k) {
                                compactGraph graph = *new compactGraph();;

                                graph.bl = prevLMin;
                                graph.br = LR[0];
                                graph.el = i + 1;
                                //R-min: The 0 element in R, which is k
                                graph.er = LR[k] + 1;

                                for (int z = 0; z < k; z++) {
                                    graph.C.push_back(LR[window + z]);
                                }
                                tempG.push_back(graph);

                                window = 1;
                                while (window + k < LR.size()) {
                                    graph = *new compactGraph();;

                                    graph.bl = LR[window - 1];
                                    graph.br = LR[window];
                                    graph.el = LR[window + k - 1];
                                    graph.er = LR[window + k];

                                    for (int z = 0; z < k; z++) {
                                        graph.C.push_back(LR[window + z]);
                                    }
                                    tempG.push_back(graph);
                                    window++;
                                }

                                graph = *new compactGraph();;

                                graph.bl = LR[window - 1];
                                graph.br = LR[window];
                                graph.el = LR[window + k - 1];
                                graph.er = prevRMax;

                                for (int z = 0; z < k; z++) {
                                    graph.C.push_back(LR[window + z]);
                                }
                                tempG.push_back(graph);
                            }
                        }
                    }
                }

            } else if (j_val > i) {
                if (j_val < RMax) {
                    // insert j into R
                    R.push_back(j_val);
                    sort(R.begin(), R.end());
                    if (R.size() == k + 1) {
                        prevRMax = R[R.size() - 1];
                        R.erase(R.end());
                        RMax = R[R.size() - 1];
                    }

                    if (L.size() + R.size() >= k && L.size() != 0 && R.size() != 0) {
                        LR = mergeLR(L, R, LR);
                        if (LR.size() == k) {
                            compactGraph graph = *new compactGraph();

                            graph.bl = prevLMin;
                            graph.br = LR[0];
                            graph.el = LR[k - 1];
                            graph.er = prevRMax;
                            for (int z = 0; z < k; z++) {
                                graph.C.push_back(LR[z]);
                            }
                            tempG.push_back(graph);
                        } else {
                            if (R.size() < k) {
                                compactGraph graph = *new compactGraph();

                                graph.bl = prevLMin;
                                graph.br = LR[0];
                                graph.el = LR[k - 1];
                                graph.er = LR[k];
                                for (int z = 0; z < k; z++) {
                                    graph.C.push_back(LR[z]);
                                }
                                tempG.push_back(graph);

                                window = 1;
                                while (window + k < LR.size()) {
                                    compactGraph graph = *new compactGraph();

                                    graph.bl = LR[window - 1];
                                    graph.br = LR[window];
                                    graph.el = LR[window + k - 1];
                                    graph.er = LR[window + k];
                                    for (int z = 0; z < k; z++) {
                                        graph.C.push_back(LR[window + z]);
                                    }
                                    tempG.push_back(graph);
                                    window++;
                                }

                                graph = *new compactGraph();;

                                graph.bl = LR[window - 1];
                                graph.br = LR[window];
                                graph.el = LR[window + k - 1];
                                graph.er = prevRMax;

                                for (int z = 0; z < k; z++) {
                                    graph.C.push_back(LR[window + z]);
                                }
                                tempG.push_back(graph);
                            } else if (R.size() == k) {
//                            if (L.size() == 0)
                                compactGraph graph = *new compactGraph();

                                graph.bl = prevLMin;
                                graph.br = LR[0];
                                graph.el = LR[k - 1];
                                graph.er = LR[k];

                                for (int z = 0; z < k; z++) {
                                    graph.C.push_back(LR[window + z]);
                                }
                                tempG.push_back(graph);

                                window = 1;
                                while (window + k < LR.size()) {
                                    graph = *new compactGraph();;

                                    graph.bl = LR[window - 1];
                                    graph.br = LR[window];
                                    graph.el = LR[window + k - 1];
                                    graph.er = LR[window + k];

                                    for (int z = 0; z < k; z++) {
                                        graph.C.push_back(LR[window + z]);
                                    }
                                    tempG.push_back(graph);
                                    window++;
                                }

                                graph = *new compactGraph();;

                                //Should be LMax, which is L.size - 1
                                graph.bl = LR[L.size() - 1] - 1;
                                graph.br = i - 1;
                                //Should be RMax
                                graph.el = LR[LR.size() - 1];
                                graph.er = prevRMax;

                                for (int z = 0; z < k; z++) {
                                    graph.C.push_back(LR[window + z]);
                                }
                                tempG.push_back(graph);
                            }
                        }
                    }
                }
            }
        }
        G.push_back(tempG);
    }
    cout << "\nFinished..." << endl;
    return G;
}


vector<pair<int, float>> calcAndSortDistance(vector<vector<float>> raw_data, vector<float> vi, int index) {
    vector<pair<int, float>> result;
    for (int i = 0; i < raw_data.size(); i++) {
        if (i != index) {
            pair<int, float> temp_pair;
            temp_pair.first = i;
            temp_pair.second = calcEuclideanDistance(raw_data[i], vi);
            result.push_back(temp_pair);
        }
    }
    sort(result.begin(), result.end(), distanceComparison);
    return result;
}

vector<pair<int, float>>
calcAndSortDistanceBF(vector<vector<float>> raw_data, vector<float> vi, int index, int x, int y) {
    vector<pair<int, float>> result;
    for (int i = x; i <= y; i++) {
        if (i != index) {
            pair<int, float> temp_pair;
            temp_pair.first = i;
            temp_pair.second = calcEuclideanDistance(raw_data[i], vi);
            result.push_back(temp_pair);
        }
    }
    sort(result.begin(), result.end(), distanceComparison);
    return result;
}


void buildPartialRangeGraph(vector<vector<float>> raw_data, int k, vector<vector<partialRangeGraph>> &Gb,
                            vector<vector<partialRangeGraph>> &Ge) {

    cout << "Indexing using Partial Range Approach..." << endl;
    int indicator = raw_data.size() / 10;

    for (int i = 0; i < raw_data.size(); i++) {
        if (i % indicator == 0) {
            cout << "#";
        }
        vector<pair<int, float>> sortedDistancePairs = calcAndSortDistance(raw_data, raw_data[i], i);
        // declare and initialize L and R
        vector<int> L;
        vector<int> R;

        int prevLMin = 0;
        int LMin = 0;
        int prevRMax = raw_data.size() + 1;
        int RMax = raw_data.size() + 1;
        L.push_back(LMin);
        R.push_back(RMax);

        vector<partialRangeGraph> tempGb;
        vector<partialRangeGraph> tempGe;

        for (int j = 0; j < sortedDistancePairs.size(); j++) {
            int j_val = sortedDistancePairs[j].first;
            if (j_val < i) {
                if (j_val > LMin) {
                    // insert j into L
                    L.push_back(j_val);
                    sort(L.begin(), L.end());
                    if (L.size() == k + 1) {
                        prevLMin = L[0];
                        L.erase(L.begin());
                        LMin = L[0];

                        partialRangeGraph graph = *new partialRangeGraph();

                        graph.bl = prevLMin;
                        graph.br = L[0];
                        for (int z = 0; z < k; z++) {
                            graph.C.push_back(L[z]);
                        }
                        tempGb.push_back(graph);
                    }
                }

            } else if (j_val > i) {
                if (j_val < RMax) {
                    // insert j into R
                    R.push_back(j_val);
                    sort(R.begin(), R.end());
                    if (R.size() == k + 1) {
                        prevRMax = R[R.size() - 1];
                        R.erase(R.end());
                        RMax = R[R.size() - 1];

                        partialRangeGraph graph = *new partialRangeGraph();

                        graph.bl = R[R.size() - 1];
                        graph.br = prevRMax;
                        for (int z = 0; z < k; z++) {
                            graph.C.push_back(R[z]);
                        }
                        tempGe.push_back(graph);
                    }

                }
            }
        }
        Gb.push_back(tempGb);
        Ge.push_back(tempGe);
    }
    cout << "\nFinished..." << endl;
}


void buildDeltaCompressionIndex(vector<vector<float>> raw_data, int k, vector<vector<deltaCompressionGraph>> &dGb,
                                vector<vector<deltaCompressionGraph>> &dGe) {

    cout << "Indexing using Delta Compression Approach..." << endl;
    int indicator = raw_data.size() / 10;

    for (int i = 0; i < raw_data.size(); i++) {
        if (i % indicator == 0) {
            cout << "#";
        }
        vector<pair<int, float>> sortedDistancePairs = calcAndSortDistance(raw_data, raw_data[i], i);
        // declare and initialize L and R
        vector<int> L;
        vector<int> R;

        int prevLMin = 0;
        int LMin = 0;
        int prevRMax = raw_data.size() + 1;
        int RMax = raw_data.size() + 1;
        L.push_back(LMin);
        R.push_back(RMax);

        vector<deltaCompressionGraph> tempGb;
        vector<deltaCompressionGraph> tempGe;

        for (int j = 0; j < sortedDistancePairs.size(); j++) {
            int j_val = sortedDistancePairs[j].first;
            if (j_val < i) {
                if (j_val > LMin) {
                    // insert j into L
                    L.push_back(j_val);
                    sort(L.begin(), L.end());
                    if (L.size() == k + 1) {
                        prevLMin = L[0];
                        L.erase(L.begin());
                        LMin = L[0];

                        deltaCompressionGraph graph = *new deltaCompressionGraph();

                        graph.lMin = prevLMin;
                        graph.j = L[0];
                        for (int z = 0; z < k; z++) {
                            graph.C.push_back(L[z]);
                        }
                        tempGb.push_back(graph);
                    }
                }

            } else if (j_val > i) {
                if (j_val < RMax) {
                    // insert j into R
                    R.push_back(j_val);
                    sort(R.begin(), R.end());
                    if (R.size() == k + 1) {
                        prevRMax = R[R.size() - 1];
                        R.erase(R.end());
                        RMax = R[R.size() - 1];

                        deltaCompressionGraph graph = *new deltaCompressionGraph();

                        graph.lMin = R[R.size() - 1];
                        graph.j = prevRMax;
                        for (int z = 0; z < k; z++) {
                            graph.C.push_back(R[z]);
                        }
                        tempGe.push_back(graph);
                    }

                }
            }
        }
        dGb.push_back(tempGb);
        dGe.push_back(tempGe);
    }
    cout << "\nFinished..." << endl;
}


vector<vector<naiveGraph>> buildBruteForceIndex(vector<vector<float>> raw_data, int k) {
    cout << "Indexing using Brute Force Approach..." << endl;
    int indicator = raw_data.size() / 10;
    vector<vector<naiveGraph>> result;
    naiveGraph graph;

    for (int i = 0; i < raw_data.size(); i++) {
        if (i % indicator == 0) {
            cout << "#";
        }
        vector<float> vi = raw_data[i];
        vector<naiveGraph> tempGraph;
        for (int x = 0; x < i; x++) {
            if (i + 1 < raw_data.size()) {
                for (int y = i + 1; y < raw_data.size(); y++) {
                    if (y - x >= k) {
                        vector<pair<int, float>> distancePairs = calcAndSortDistanceBF(raw_data, vi, i, x, y);
                        graph = *new naiveGraph();
                        graph.b = x;
                        graph.e = y;
                        for (int j = 0; j < k; j++) {
                            graph.C.push_back(distancePairs[j].first);
                        }
                        tempGraph.push_back(graph);
                    }
                }
            }
        }
        result.push_back(tempGraph);
    }
    cout << "\nFinished..." << endl;
    return result;
}


vector<vector<float>> cutData(vector<vector<float>> originalData, int x, int y) {
    vector<vector<float>> result;

    for (int i = x; i <= y; i++) {
        result.push_back(originalData[i]);
    }

    return result;
}


// this function is for sorting the key, value pair
bool distanceComparison(const pair<int, float> &a, const pair<int, float> &b) {
    return a.second < b.second;
}


float calcEuclideanDistance(vector<float> v1, vector<float> v2) {
    float result = 0.0;
    for (int i = 0; i < v1.size(); i++) {
        result += pow((v1[i] - v2[i]), 2);
    }
    result = sqrt(result);
    return result;
}

vector<int> mergeLR(vector<int> L, vector<int> R, vector<int> LR) {
    LR.clear();

    for (int i = 0; i < L.size(); i++) {
        LR.push_back(L[i]);
    }

    for (int j = 0; j < R.size(); j++) {
        LR.push_back(R[j]);
    }

    return LR;
}

long countGraphSize(vector<vector<compactGraph>> G) {
    long result = 0;
    for (int i = 0; i < G.size(); i++) {
        result += G[i].size();
    }
    return result;
}

long countBruteForceSize(vector<vector<naiveGraph>> G) {
    long result = 0;
    for (int i = 0; i < G.size(); i++) {
        result += G[i].size();
    }
    return result;
}

vector<int> buildKNNPartialRange(vector<vector<float>> rawData, vector<vector<partialRangeGraph>> Gb,
                                 vector<vector<partialRangeGraph>> Ge, int x, int y, int i, int k) {
    vector<int> result;
    vector<int> tempNN;
    vector<float> vi = rawData[i];
    vector<pair<int, float>> distancePairs;

    if (Gb[i].size() != 0) {
        for (int j = 0; j < Gb[i].size(); j++) {
            if (x > Gb[i][j].bl && x <= Gb[i][j].br) {
                for (int z = 0; z < Gb[i][j].C.size(); z++) {
                    tempNN.push_back(Gb[i][j].C[z]);
                }
                break;
            }
        }
    }

    if (Ge.size() != 0) {
        for (int j = 0; j < Ge[i].size(); j++) {
            if (y >= Ge[i][j].bl && y < Ge[i][j].br) {
                for (int z = 0; z < Ge[i][j].C.size(); z++) {
                    tempNN.push_back(Ge[i][j].C[z]);
                }
                break;
            }
        }
    }

    for (int j = 0; j < tempNN.size(); j++) {
        pair<int, float> tempPair;
        tempPair.first = tempNN[j];
        tempPair.second = calcEuclideanDistance(vi, rawData[tempNN[j]]);
        distancePairs.push_back(tempPair);
    }
    sort(distancePairs.begin(), distancePairs.end(), distanceComparison);

    for (int j = 0; j < k; j++) {
        result.push_back(distancePairs[j].first);
    }
    return result;
}

vector<int> buildKNNCompactList(vector<vector<compactGraph>> G, int x, int y, int i) {
    vector<int> result;
    if (G[i].size() == 0) {
        return result;
    }
    for (int j = 0; j < G[i].size(); j++) {
        if (x > G[i][j].bl && x <= G[i][j].br && y >= G[i][j].el && y < G[i][j].er) {
            return G[i][j].C;
        }
    }
    return result;
}

vector<int> buildKNNBruteForce(vector<vector<naiveGraph>> G, int x, int y, int i) {
    vector<int> result;
    if (G[i].size() == 0) {
        return result;
    }
    for (int j = 0; j < G[i].size(); j++) {
        if (x == G[i][j].b && y == G[i][j].e) {
            return G[i][j].C;
        }
    }
    return result;
}

vector<int> buildKNNDeltaCompression(vector<vector<float>> rawData, vector<vector<deltaCompressionGraph>> dGb,
                                 vector<vector<deltaCompressionGraph>> dGe, int x, int y, int i, int k) {
    vector<int> result;
    vector<int> tempNN;
    vector<float> vi = rawData[i];
    vector<pair<int, float>> distancePairs;

    if (dGb[i].size() != 0) {
        for (int j = 0; j < dGb[i].size(); j++) {
            if (x > dGb[i][j].lMin && x <= dGb[i][j].j) {
                for (int z = 0; z < dGb[i][j].C.size(); z++) {
                    tempNN.push_back(dGb[i][j].C[z]);
                }
                break;
            }
        }
    }

    if (dGe.size() != 0) {
        for (int j = 0; j < dGe[i].size(); j++) {
            if (y >= dGe[i][j].lMin && y < dGe[i][j].j) {
                for (int z = 0; z < dGe[i][j].C.size(); z++) {
                    tempNN.push_back(dGe[i][j].C[z]);
                }
                break;
            }
        }
    }

    for (int j = 0; j < tempNN.size(); j++) {
        pair<int, float> tempPair;
        tempPair.first = tempNN[j];
        tempPair.second = calcEuclideanDistance(vi, rawData[tempNN[j]]);
        distancePairs.push_back(tempPair);
    }
    sort(distancePairs.begin(), distancePairs.end(), distanceComparison);

    for (int j = 0; j < k; j++) {
        result.push_back(distancePairs[j].first);
    }
    return result;
}
