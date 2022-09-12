#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <assert.h>
#include "readCSV.h"

using namespace std;

void dataPreprocessing::readCSV(string addr)
{
    ifstream data(addr, ios::in);
    assert(data.is_open() == true);

    // 跳过标题行
    string line;
    getline(data, line);

    while (getline(data, line))
    {
        vector<double> tmp;
        istringstream customer(line);
        string num;
        // 跳过行号(第0列)
        getline(customer, num, ',');
        // 跳过样本ID(第1列)
        getline(customer, num, ',');
        // 读取标签
        getline(customer, num, ',');
        labels.push_back(static_cast<long>(stod(num)));
        while (getline(customer, num, ','))
        {
            tmp.push_back(stod(num));
        }
        datasets.push_back(tmp);
    }

    data.close();

    num_customers = datasets.size();
    featureNum = datasets[0].size();
}

void dataPreprocessing::partitionFeatures(double coefficient)
{
    unsigned A_features = static_cast<unsigned>(featureNum * coefficient);
    for (int i = 0; i < num_customers; ++i)
    {
        vector<double> tempA;
        vector<double> tempB;
        for (int j = 0; j < A_features; ++j)
        {
            tempA.push_back(datasets[i][j]);
        }
        for (int j = A_features; j < featureNum; ++j)
        {
            tempB.push_back(datasets[i][j]);
        }
        A.push_back(tempA);
        B.push_back(tempB);
    }
}

void dataPreprocessing::normalize(int column)
{
    double max = 0;
    for (int row = 0; row < num_customers; ++row)
    {
        if (datasets[row][column] > max)
            max = datasets[row][column];
    }
    for (int row = 0; row < num_customers; ++row)
    {
        datasets[row][column] /= max;
    }
}
