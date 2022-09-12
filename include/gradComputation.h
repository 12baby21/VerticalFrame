/**
 * @file gradComputation.h
 * @author Michael Wu (18917131713@163.com)
 * @brief Compute plaintext gradient for one-level LR model.
 * @version 0.1
 * @date 2022-09-11
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#pragma once

#include <iostream>
#include <gmpxx.h>
#include <vector>
#include <cmath>

using namespace std;


namespace gradComputation
{
    double myCrossEntropyLoss(vector<vector<double>> &prediction, vector<vector<long>> &label);
    void computeDiff(vector<vector<double>> &diff, vector<vector<long>> &labels, vector<vector<double>> &prediction);
    void CalGradWeight(vector<vector<double>> &gradWeight, vector<vector<double>> &diff, vector<vector<double>> &input);
    void CalGradBias(vector<double> &gradBias, vector<vector<double>> &diff);
    void mySGDUpdateWeight(vector<double> &weight, vector<vector<double>> &grad, double lr);
    void mySGDUpdateBias(vector<double> &bias, vector<double> &grad, double lr);
}