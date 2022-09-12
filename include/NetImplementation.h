#pragma once

#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

namespace NetImpementation
{
    // 线性层
    void myLinear(vector<vector<double>> &fc_out, vector<vector<double>> &x, vector<double> &weight, vector<double> &bias);
    // Softmax层
    void mySigmoid(vector<vector<double>> &out, vector<vector<double>> &fc_out);
    // 预测结果聚合
    void resultIntegration(vector<vector<double>> &res, vector<vector<double>> &A, vector<vector<double>> &B, double coef = 1);
}