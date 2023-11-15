//
// Created by Jackie on 10/4/2023.
// 数据读取: vector<vector<float>>
// ifstream很可能不行

#include <vector>
#include "iostream"
#include "string"
#include "fstream"

using namespace std;

void ReadDataWrapper(vector<vector<float>> &raw_data, vector<int> &search_keys,
                     const string &dataset, string &dataset_path,
                     const int item_num);

int main() {
    ifstream in1;
    string data_path = "../data/siftsmall_base.fvecs";
    in1.open(data_path);
    int count = 0;
    char ch;

    while (count < 5) {
        in1.get(ch);
        cout << (ch - 48) << endl;
        count++;
    }

    return 0;
}

//void ReadDataWrapper(vector<vector<float>> &raw_data, vector<int> &search_keys,
//                     const string &dataset, string &dataset_path,
//                     const int item_num) {
//    raw_data.clear();
////    if (dataset == "glove") {
////        ReadMatFromTxtTwitter(dataset_path, raw_data, item_num);
////    } else if (dataset == "ml25m") {
////        ReadMatFromTxt(dataset_path, raw_data, item_num);
////    } else if (dataset == "sift") {
////        raw_data = pqdescent::ReadTopN(dataset_path, "bvecs", item_num);
////    } else if (dataset == "biggraph") {
////        ReadMatFromTsv(dataset_path, raw_data, item_num);
////    } else if (dataset == "local") {
////        if (dataset_path == "")
//    dataset_path = "../data/siftsmall_base.fvecs";
//    cout << dataset_path << endl;
//    raw_data = pqdescent::ReadTopN(dataset_path, "fvecs", item_num);
////    } else if (dataset == "deep1b") {
////        raw_data = pqdescent::ReadTopN(dataset_path, "fvecs", item_num);
////    } else if (dataset == "deep10m") {
////        raw_data = pqdescent::ReadTopN(dataset_path, "fvecs", item_num);
////    } else if (dataset == "yt8m") {
////        ReadMatFromTsvYT8M(dataset_path, raw_data, search_keys, item_num);
////    } else {
////        std::cerr << "Wrong Datset!" << endl;
////        assert(false);
////    }
//}
//
//std::vector<std::vector<float>> ReadTopN(std::string filename, std::string ext, int top_n) {
//    std::vector<std::vector<float>> vecs;
//    if (top_n != -1) {
//        vecs.reserve(top_n);
//    }
//    ItrReader reader(filename, ext);
//    while (!reader.IsEnd()) {
//        if (top_n != -1 && top_n <= (int) vecs.size()) {
//            break;
//        }
//        vecs.push_back(reader.Next());
//    }
//    return vecs;
//}

