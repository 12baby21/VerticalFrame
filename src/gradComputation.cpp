#include <iostream>
#include <gmpxx.h>
#include <vector>
#include "util.h"
#include <cmath>
#include <assert.h>

namespace gradComputation
{
    double myCrossEntropyLoss(vector<vector<double>> &prediction, vector<vector<long>> &labels)
    {
        double loss = 0;
        int batchSize = prediction.size();
        for (int i = 0; i < batchSize; ++i)
        {
            loss = labels[i][1] * log(prediction[i][1]) + labels[i][0] * log(prediction[i][0]);
        }
        loss = -1 * loss;

        return loss;
    }

    void computeDiff(vector<vector<double>> &diff, vector<vector<long>> &labels, vector<vector<double>> &prediction)
    {
        assert(labels.size() == prediction.size());
        assert(labels[0].size() == prediction[0].size());
        int batchSize = diff.size();
        int l = labels[0].size();
        for (int i = 0; i < batchSize; ++i)
        {
            for (int j = 0; j < l; ++j)
            {
                diff[i][j] = labels[i][j] - prediction[i][j];
                diff[i][j] = -1 * diff[i][j];
            }
        }
    }

    void CalGradWeight(vector<vector<double>> &gradWeight, vector<vector<double>> &diff, vector<vector<double>> &input)
    {
        // 梯度一定要初始化为0
        int batchSize = diff.size();
        int featureNum = input[0].size();
        int outDimension = diff[0].size();
        for (int i = 0; i < batchSize; ++i)
        {
            for (int j = 0; j < outDimension; ++j)
            {
                for (int k = 0; k < featureNum; ++k)
                {
                    gradWeight[j][k] += diff[i][j] * input[i][k];
                }
            }
        }
        for (int i = 0; i < outDimension; ++i)
        {
            for (int j = 0; j < featureNum; ++j)
            {
                gradWeight[i][j] /= batchSize;
            }
        }
    }

    void CalGradBias(vector<double> &gradBias, vector<vector<double>> &diff)
    {
        int batchSize = diff.size();
        int outDimension = diff[0].size();
        for (int i = 0; i < batchSize; ++i)
        {
            for (int j = 0; j < outDimension; ++j)
            {
                gradBias[j] += diff[i][j];
            }
        }
        for (int i = 0; i < outDimension; ++i)
        {
            gradBias[i] /= batchSize;
        }
    }

    void mySGDUpdateWeight(vector<double> &weight, vector<vector<double>> &grad, double lr)
    {
        int featureNum = grad[0].size();
        for (int i = 0; i < grad.size(); ++i)
        {
            for (int j = 0; j < grad[0].size(); ++j)
            {
                weight[i * featureNum + j] = weight[i * featureNum + j] - lr * grad[i][j];
            }
        }
    }

    void mySGDUpdateBias(vector<double> &bias, vector<double> &grad, double lr)
    {
        for (int i = 0; i < bias.size(); ++i)
        {
            bias[i] = bias[i] - lr * grad[i];
        }
    }

}

