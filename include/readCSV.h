#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

class dataPreprocessing
{
public:
    vector<vector<double>> A;
    vector<vector<double>> B;
    vector<long> labels;
    vector<vector<long>> newLabels;

    dataPreprocessing(string in_addr, double coefficient)
    {
        addr = in_addr;
        readCSV(addr);
        num_customers = datasets.size();
        normalize(0);
        partitionFeatures(coefficient);
    }

    void normalize(int column);

private:
    int num_customers;  // 数据集大小
    int featureNum;     // 总特征数量
    string addr;
    vector<vector<double>> datasets;
    void readCSV(string addr);
    void partitionFeatures(double coefficient);
};