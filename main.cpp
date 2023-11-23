//
// Created by Jackie on 10/4/2023.
// 数据读取: vector<vector<float>>
// ifstream很可能不行

#include <vector>
#include "iostream"
#include "string"
#include "fstream"
#include <cmath>
#include <algorithm>

using namespace std;

struct indexGraph {
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

void loadData(const string &filename, vector<vector<float> > &raw_data, unsigned &num, unsigned &dim);

float calcEuclideanDistance(vector<float> v1, vector<float> v2);

bool distanceComparison(const pair<int, float> &a, const pair<int, float> &b);

vector<vector<naiveGraph>> buildBruteForceIndex(vector<vector<float>> raw_data, int k);

vector<vector<indexGraph>> buildCompactGraph(vector<vector<float>> raw_data, int k);

vector<pair<int, float>> calcAndSortDistance(vector<vector<float>> raw_data, vector<float> vi, int index);

vector<vector<float>> cutData(vector<vector<float>> originalData, int x, int y);

long countGraphSize(vector<vector<indexGraph>> G);

vector<int> mergeLR(vector<int> L, vector<int> R, vector<int> LR);

long countBruteForceSize(vector<vector<naiveGraph>> G);

void buildPartialRangeGraph(vector<vector<float>> raw_data, int k, vector<vector<partialRangeGraph>> &Gb,
                            vector<vector<partialRangeGraph>> &Ge);

vector<int> buildKNNPartialRange(vector<vector<float>> rawData, vector<vector<partialRangeGraph>> Gb,
                                 vector<vector<partialRangeGraph>> Ge, int x, int y, int i, int k);

vector<int> buildKNNCompactList(vector<vector<indexGraph>> G, int x, int y, int i);

vector<int> buildKNNBruteForce(vector<vector<naiveGraph>> G, int x, int y, int i);

int main() {
    int k = 10;
    string data_path = "data/siftsmall_base.fvecs";
    vector<vector<float> > raw_data;
    unsigned dim, num;
    loadData(data_path, raw_data, num, dim);
    raw_data.resize(100);

    cout << "Program Begin. Data size: " << raw_data.size() << " K: " << k << endl;
    vector<vector<naiveGraph>> bruteForceIndex = buildBruteForceIndex(raw_data, k);
    vector<vector<indexGraph>> compactIndex = buildCompactGraph(raw_data, k);

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

    cout << "Program Finished... Calculating Stats..." << endl;

    int bruteForceNumber = countBruteForceSize(bruteForceIndex);
    long compactGraphNumber = countGraphSize(compactIndex);

    cout << "Brute Force Index Number: " << bruteForceNumber << "\tTotal Size: "
         << (bruteForceNumber * (8 + 4 * k)) / (1024.0 * 1024.0) << " MiB" << endl;
    cout << "Compact Graph Index Number: " << compactGraphNumber << "\tTotal Size: "
         << (compactGraphNumber * (16 + 4 * k)) / (1024.0 * 1024.0) << " MiB" << endl;
    cout << "Partial Range Index Number: " << partialRangeNumber << "\tTotal Size: "
         << (partialRangeNumber * (8 + 4 * k)) / (1024.0 * 1024.0) << " MiB" << endl;


    vector<int> bfknn = buildKNNBruteForce(bruteForceIndex, 12, 66, 30);
    vector<int> cgknn = buildKNNCompactList(compactIndex, 12, 66, 30);
    vector<int> prknn = buildKNNPartialRange(raw_data, Gb, Ge, 12, 66, 30, k);

    cout << "BF\tCG\tPR" << endl;
    for (int i = 0; i < k; i++) {
        cout << bfknn[i] << "\t" << cgknn[i] << "\t" << prknn[i] << endl;
    }

//    cout << "CG\tPR" << endl;
//    for (int i = 0; i < k; i ++) {
//        cout << cgknn[i] << "\t" << prknn[i] << endl;
//    }
    return 0;
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


vector<vector<indexGraph>> buildCompactGraph(vector<vector<float>> raw_data, int k) {

    cout << "Indexing using Compact Graph Approach..." << endl;
    int indicator = raw_data.size() / 100;
    vector<vector<indexGraph>> G;
    int window;

    for (int i = 0; i < raw_data.size(); i++) {
        if (i % indicator == 0) {
            cout << "#";
        }

        if (i == 5) {
            cout << "sdf";
        }


        vector<pair<int, float>> sortedDistancePairs = calcAndSortDistance(raw_data, raw_data[i], i);
        // declare and initialize L and R
        vector<int> L;
        vector<int> R;
        vector<int> LR;
        vector<indexGraph> tempG;

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
                    if (L.size() == k + 1) {
                        prevLMin = L[0];
                        L.erase(L.begin());
                        LMin = L[0];
                    }
                    sort(L.begin(), L.end());
                    // merge L and R
                    if (L.size() + R.size() >= k && L.size() != 0 && R.size() != 0) {
                        LR = mergeLR(L, R, LR);
                        //cout << LR.size() << endl;

                        if (L.size() < k) {
                            indexGraph graph = *new indexGraph();

                            graph.bl = prevLMin;
                            graph.br = LR[0];
                            graph.el = LR[k - 1];
                            graph.er = LR[k];
                            for (int z = 0; z < k; z++) {
                                graph.C.push_back(LR[z]);
                            }
                            tempG.push_back(graph);

                            window = 0;

                            while (window + k < LR.size()) {
                                indexGraph graph = *new indexGraph();

                                graph.bl = LR[window];
                                graph.br = LR[window + 1];
                                graph.el = LR[window + k - 1];
                                graph.er = LR[window + k];
                                for (int z = 0; z < k; z++) {
                                    graph.C.push_back(LR[window + z]);
                                }
                                tempG.push_back(graph);
                                window++;
                            }

                            graph = *new indexGraph();;

                            graph.bl = LR[window];
                            graph.br = LR[window + 1];
                            graph.el = LR[window + k - 1];
                            graph.er = prevRMax;

                            for (int z = 0; z < k; z++) {
                                graph.C.push_back(LR[window + z]);
                            }
                            tempG.push_back(graph);

                        } else if (L.size() == k) {
                            indexGraph graph = *new indexGraph();;

                            graph.bl = prevLMin;
                            graph.br = LR[0];
                            graph.el = i - 1;
                            graph.er = LR[k];

                            for (int z = 0; z < k; z++) {
                                graph.C.push_back(LR[window + z]);
                            }
                            tempG.push_back(graph);

                            window = 0;
                            while (window + k < LR.size()) {
                                graph = *new indexGraph();;

                                graph.bl = LR[window];
                                graph.br = LR[window + 1];
                                graph.el = LR[window + k - 1];
                                graph.er = LR[window + k];

                                for (int z = 0; z < k; z++) {
                                    graph.C.push_back(LR[window + z]);
                                }
                                tempG.push_back(graph);
                                window++;
                            }

                            graph = *new indexGraph();;

                            graph.bl = LR[window];
                            graph.br = LR[window + 1];
                            graph.el = LR[window + k - 1];
                            graph.er = prevRMax;

                            for (int z = 0; z < k; z++) {
                                graph.C.push_back(LR[window + z]);
                            }
                            tempG.push_back(graph);
                        }
                    }
                }

            } else if (j_val > i) {
                if (j_val < RMax) {
                    // insert j into R
                    R.push_back(j_val);
                    if (R.size() == k + 1) {
                        prevRMax = R[R.size() - 1];
                        R.erase(R.end());
                        RMax = R[R.size() - 1];
                    }
                    sort(R.begin(), R.end());

                    if (L.size() + R.size() >= k && L.size() != 0 && R.size() != 0) {
                        LR = mergeLR(L, R, LR);
//                        cout << LR.size() << endl;
                        if (R.size() < k) {
                            indexGraph graph = *new indexGraph();

                            graph.bl = prevLMin;
                            graph.br = LR[0];
                            graph.el = LR[k - 1];
                            graph.er = LR[k];
                            for (int z = 0; z < k; z++) {
                                graph.C.push_back(LR[z]);
                            }
                            tempG.push_back(graph);

                            window = 0;
                            while (window + k + 1 < LR.size()) {
                                indexGraph graph = *new indexGraph();

                                graph.bl = LR[window];
                                graph.br = LR[window + 1];
                                graph.el = LR[window + k - 1];
                                graph.er = LR[window + k];
                                for (int z = 0; z < k; z++) {
                                    graph.C.push_back(LR[window + z]);
                                }
                                tempG.push_back(graph);
                                window++;
                            }

                            graph = *new indexGraph();;

                            graph.bl = LR[window];
                            graph.br = LR[window + 1];
                            graph.el = LR[window + k - 1];
                            graph.er = prevRMax;

                            for (int z = 0; z < k; z++) {
                                graph.C.push_back(LR[window + z]);
                            }
                            tempG.push_back(graph);
                        } else if (R.size() == k) {
//                            if (L.size() == 0)
                            indexGraph graph = *new indexGraph();

                            graph.bl = prevLMin;
                            graph.br = LR[0];
                            graph.el = LR[k - 1];
                            graph.er = LR[k];

                            for (int z = 0; z < k; z++) {
                                graph.C.push_back(LR[window + z]);
                            }
                            tempG.push_back(graph);

                            window = 0;
                            while (window + k < LR.size()) {
                                graph = *new indexGraph();;

                                graph.bl = LR[window];
                                graph.br = LR[window + 1];
                                graph.el = LR[window + k - 1];
                                graph.er = LR[window + k];

                                for (int z = 0; z < k; z++) {
                                    graph.C.push_back(LR[window + z]);
                                }
                                tempG.push_back(graph);
                                window++;
                            }

                            graph = *new indexGraph();;

                            graph.bl = LR[LR.size() - k - 1];
                            graph.br = i + 1;
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
        }
        if (tempG.size() != 0) {
            G.push_back(tempG);
        }
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


void buildPartialRangeGraph(vector<vector<float>> raw_data, int k, vector<vector<partialRangeGraph>> &Gb,
                            vector<vector<partialRangeGraph>> &Ge) {

    cout << "Indexing using Partial Range Approach..." << endl;
    int indicator = raw_data.size() / 100;

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


//vector<vector<int>> buildBruteForceIndex(vector<vector<float>> raw_data, int k) {
//    vector<vector<int>> result;
//
//    for (int i = 0; i < raw_data.size(); i++) {
//        vector<pair<int, float>> temp_pairs;
//        for (int j = 0; j < raw_data.size(); j++) {
//            if (i != j) {
//                pair<int, float> temp_pair;
//                temp_pair.first = j;
//                temp_pair.second = calcEuclideanDistance(raw_data[i], raw_data[j]);
//                temp_pairs.push_back(temp_pair);
//            }
//        }
//        // sort the key-value pair by distance
//        sort(temp_pairs.begin(), temp_pairs.end(), distanceComparison);
//        vector<int> temp_res;
//        // push k indexes into the result
//        for (int x = 0; x < k; x++) {
//            // if k is greater than the total search key, push k values
//            if (x < temp_pairs.size()) {
//                temp_res.push_back(temp_pairs[x].first);
//            }
//        }
//        result.push_back(temp_res);
//    }
//
//    return result;
//}

vector<vector<naiveGraph>> buildBruteForceIndex(vector<vector<float>> raw_data, int k) {
    cout << "Indexing using Brute Force Approach..." << endl;
    int indicator = raw_data.size() / 100;
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
                        vector<pair<int, float>> distancePairs = calcAndSortDistance(cutData(raw_data, x, y), vi,
                                                                                     i);
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

long countGraphSize(vector<vector<indexGraph>> G) {
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
            if (x >= Gb[i][j].bl && x <= Gb[i][j].br) {
                for (int z = 0; z < Gb[i][j].C.size(); z++) {
                    tempNN.push_back(Gb[i][j].C[z]);
                }
                break;
            }
        }
    }

    if (Ge.size() != 0) {
        for (int j = 0; j < Ge[i].size(); j++) {
            if (y >= Ge[i][j].bl && y <= Ge[i][j].br) {
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

vector<int> buildKNNCompactList(vector<vector<indexGraph>> G, int x, int y, int i) {
    vector<int> result;
    if (G[i].size() == 0) {
        return result;
    }
    for (int j = 0; j < G[i].size(); j++) {
        if (x >= G[i][j].bl && x <= G[i][j].br && y >= G[i][j].el && y <= G[i][j].er) {
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
