#include <iostream>
#include <vector>
#include <cmath>
#include <modelCreation.h>
using namespace std;

class Net {
private:
    int batchSize;
    int featureNum;
    int outDimension;
    vector<double> weight;
    vector<double> bias;

public:
    Net(int _batchSize, int _featureNum, int _outDimension) {
        batchSize = batchSize;
        featureNum = _featureNum;
        outDimension = _outDimension;
        weight = vector<double>(2 * featureNum, 0);
        bias = vector<double>(2, 0);
        createModel<double>(weight);
    }

    void linear(vector<vector<double>>& fc_out, vector<vector<double>>& x) {
        for (int i = 0; i < batchSize; ++i)
        {
            for (int j = 0; j < outDimension; ++j)
            {
                double tmp = 0;
                for (int k = 0; k < featureNum; ++k)
                {
                    tmp += weight[j * featureNum + k] * x[i][k];
                }
                tmp += bias[j];
                fc_out[i][j] = tmp;
            }
        }
    }

    void sigmoid(vector<vector<double>> &out, vector<vector<double>> &fc_out)
    {
        int batchSize = fc_out.size();
        int outDimension = fc_out[0].size();
        for (int i = 0; i < batchSize; ++i)
        {
            for (int j = 0; j < outDimension; ++j)
            {
                out[i][j] = 1 / (1 + exp(-fc_out[i][j]));
            }
        }
    }

    void forward(vector<vector<double>>& out, vector<vector<double>>& x) {
        linear(out, x);
        sigmoid(out, out);
    }
};

namespace NetImpementation
{
    void myLinear(vector<vector<double>> &fc_out, vector<vector<double>> &x, vector<double> &weight, vector<double> &bias)
    {
        int batchSize = x.size();
        int featureNum = x[0].size(); // 特征数量
        int outDimension = fc_out[0].size();
        for (int i = 0; i < batchSize; ++i)
        {
            for (int j = 0; j < outDimension; ++j)
            {
                double tmp = 0;
                for (int k = 0; k < featureNum; ++k)
                {
                    tmp += weight[j * featureNum + k] * x[i][k];
                }
                tmp += bias[j];
                fc_out[i][j] = tmp;
            }
        }
    }

    void mySigmoid(vector<vector<double>> &out, vector<vector<double>> &fc_out)
    {
        int batchSize = fc_out.size();
        int outDimension = fc_out[0].size();
        for (int i = 0; i < batchSize; ++i)
        {
            for (int j = 0; j < outDimension; ++j)
            {
                out[i][j] = 1 / (1 + exp(-fc_out[i][j]));
            }
        }
    }

    void resultIntegration(vector<vector<double>> &res, vector<vector<double>> &A, vector<vector<double>> &B, double coef = 1)
    {
        int n = res.size();
        int m = 2;
        for (int i = 0; i < n; ++i)
        {
            double sum = A[i][0] + B[i][0] + A[i][1] + B[i][1];

            res[i][0] = (A[i][0] + B[i][0]) / sum;
            res[i][1] = (A[i][1] + B[i][1]) / sum;
        }
    }
}
