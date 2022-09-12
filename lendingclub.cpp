/**
 * Function Name
 * @ baseline: Encryption_QEXP
 * @ SBR: Encryption_SBR
 * @ pre: Encryption_pre
 * @ pre+SBR: Encryption_pre_SBR
 */
#include <iostream>
#include <string>
#include <fstream>
#include <gmpxx.h>
#include "util.h"
#include <ctime>
#include <chrono>
#include <thread>
#include "readCSV.h"
#include "modelCreation.h"
#include "NetImplementation.h"
#include "gradComputation.h"
// 加密算法
#include "paillier_base.h"
#include "paillier_opt.h"
#include "paillier_djn.h"
#include "paillier_ant.h"

using namespace std;
using namespace std::chrono;

// DEBUG MACRO
// #define DEBUG_MODEL_CREATION
// #define DEBUG_PROCESS
// #define MNIST
// #define TEST_SIZE
// #define TEST_TIME

// global variables
mpz_class R;
mpz_class zero(0);
mpz_class enZero;

int main()
{
    // 系统初始化，bits是p和q的位宽，n_length用于ANT
    mp_bitcnt_t bits = 1024;
    mp_bitcnt_t n_length = 2 * bits;

    // Encryption Protocol
    // BASE Proto(bits);
    // OPT Proto(bits);
    // DJN Proto(bits);
    ANT Proto(n_length);

    // double testA = -0.4957;
    // mpz_class res;
    // Proto.Encode(res, testA, 1e6);
    // gmp_printf("%Zd\n", res);
    // Proto.Encryption_QEXP(res, res);
    // Proto.Decryption(res, res);
    // cout << "testA = " << testA << endl;

    // return 0;

    // Please Use Encryption_QEXP as The Baseline
    R = myGenRand::GenRandomPrime(256);
    mpz_class mask;
    Proto.Encryption_QEXP(mask, R);
    Proto.Encryption_QEXP(enZero, zero);

    double lr = 0.1;
    std::ofstream sfile("../trainLog.txt", ios::out);
    srand((int)time(0));

    // 加载并预处理图片和标签
    const string addr = "../Financial_data/lendingclub_11w_minmax_test.csv";
    dataPreprocessing mydata(addr, 0.5);

    // 训练超参数
    int batchNum = mydata.A.size() / batchSize;
    int featureNumA = mydata.A[0].size();
    int featureNumB = mydata.B[0].size();

    // 模型初始化
    vector<double> weightA(2 * featureNumA);
    vector<double> biasA(2, 0);
    vector<double> weightB(2 * featureNumB);
    vector<double> biasB(2, 0);
    createModel<double>(weightA);
    createModel<double>(weightB);

    

    // 开始训练
    for (size_t epoch = 1; epoch <= 1; ++epoch)
    {
        double lossSum = 0;
        // Iterate the data loader to yield batches from the dataset.
        for (int batchIndex = 0; batchIndex < batchNum; ++batchIndex)
        {
            auto iter_st = clock();
            // 准备batch数据
            vector<vector<double>> batchA;
            vector<vector<double>> batchB;
            vector<vector<long>> labelA;
            for (int i = 0; i < batchSize; ++i)
            {
                batchA.push_back(mydata.A[batchIndex * batchSize + i]);
                batchB.push_back(mydata.B[batchIndex * batchSize + i]);
                vector<long> tmp(2);
                if (mydata.labels[batchIndex * batchSize + i] == 0)
                {
                    tmp[0] = 1;
                    tmp[1] = 0;
                }
                else
                {
                    tmp[0] = 0;
                    tmp[1] = 1;
                }
                labelA.push_back(tmp);
            }

            // Execute the model on the input data.
            vector<vector<double>> fc_outA(batchSize, vector<double>(2));
            vector<vector<double>> fc_outB(batchSize, vector<double>(2));
            vector<vector<double>> prediction(batchSize, vector<double>(2));
            NetImpementation::myLinear(fc_outA, batchA, weightA, biasA);
            NetImpementation::myLinear(fc_outB, batchB, weightB, biasB);
            NetImpementation::mySigmoid(fc_outA, fc_outA);
            NetImpementation::mySigmoid(fc_outB, fc_outB);
            NetImpementation::resultIntegration(prediction, fc_outA, fc_outB, 0.5);
#ifdef DEBUG_PROCESS
            cout << "前向传播完成" << endl;
#endif

            // 计算loss
            double loss = gradComputation::myCrossEntropyLoss(prediction, labelA);
            lossSum += loss;
#ifdef DEBUG_PROCESS
            cout << "loss计算完成" << endl;
            cout << "loss = " << loss << endl;
#endif

            // 计算diff
            vector<vector<double>> diff(batchSize, vector<double>(2));
            gradComputation::computeDiff(diff, labelA, prediction);
#ifdef DEBUG_PROCESS
            cout << "diff计算完成" << endl;
#endif
#ifdef TEST_SIZE
            cout << diff.size() << endl;
            cout << diff[0].size() << endl;
#endif

            // 计算参与方A的梯度
            vector<vector<double>> gradWeightA(2, vector<double>(featureNumA, 0));
            vector<double> gradBiasA(2);
            gradComputation::CalGradWeight(gradWeightA, diff, batchA);
            gradComputation::CalGradBias(gradBiasA, diff);
            // 这里应该不需要除batchsize
            for (int i = 0; i < 2; ++i)
            {
                gradBiasA[i] /= batchSize;
                for (int j = 0; j < featureNumA; ++j)
                {
                    gradWeightA[i][j] /= batchSize;
                }
            }
#ifdef DEBUG_PROCESS
            cout << "服务发起方的梯度计算完成" << endl;
#endif

            // PartyA把diff加密后传给B
            Proto.resetTime();
            vector<vector<mpz_class>> enDiff(batchSize, vector<mpz_class>(2));
            // 编码后的第0,1维数据
            for (int i = 0; i < batchSize; ++i)
            {
                for (int j = 0; j < 2; ++j)
                {
                    Proto.Encode(enDiff[i][j], diff[i][j], 1e6);
                    Proto.Encryption_SBR(enDiff[i][j], enDiff[i][j]);
                }
            }
#ifdef DEBUG_PROCESS
            cout << "编码和加密diff完成" << endl;
#endif
            printf("加密diff花费%lfs.\n", Proto.encryptionTime);
            
            // 输入编码
            int row = 2, column = featureNumB;
            vector<vector<mpz_class>> mpz_input(batchSize, vector<mpz_class>(column));
            for (int i = 0; i < batchSize; ++i)
            {
                for (int j = 0; j < batchB[0].size(); ++j)
                {
                    Proto.Encode(mpz_input[i][j], batchB[i][j], 1e6);
                }
            }
#ifdef DEBUG_PROCESS
            cout << "编码input完成" << endl;
#endif

            // [weight * diff + R]
            double gradComputationTime = 0;
            auto gradComputationStart = system_clock::now();
            std::vector<mpz_class> encrypt_gradWeight_B(row * column, enZero);
            for (int i = 0; i < row; ++i)
            {
                for (int j = 0; j < column; ++j)
                {
                    for (int k = 0; k < batchSize; ++k)
                    {
                        mpz_class tmp;
                        Proto.EncryptMul(tmp, enDiff[k][i], mpz_input[k][j]);
                        Proto.EncryptAdd(encrypt_gradWeight_B[i * column + j], encrypt_gradWeight_B[i * column + j], tmp);
                    }
                }
            }
            for (int i = 0; i < row; ++i)
            {
                Proto.EncryptAdd(encrypt_gradWeight_B[i], encrypt_gradWeight_B[i], mask);
            }
#ifdef DEBUG_PROCESS
            cout << "[weight * diff + R]完成" << endl;
#endif

            // [bias * diff + R]
            std::vector<mpz_class> encrypt_gradBias_B(row, enZero);
            for (int i = 0; i < row; ++i)
            {
                for (int j = 0; j < batchSize; ++j)
                {
                    Proto.EncryptAdd(encrypt_gradBias_B[i], encrypt_gradBias_B[i], enDiff[j][i]);
                }
            }
            for (int i = 0; i < row; ++i)
            {
                Proto.EncryptAdd(encrypt_gradBias_B[i], encrypt_gradBias_B[i], mask);
            }
            auto gradComputationEnd = system_clock::now();
            auto duration = double(duration_cast<microseconds>(gradComputationEnd - gradComputationStart).count());
            gradComputationTime = duration * microseconds::period::num / microseconds::period::den;
            printf("B在密文域上计算梯度花费了%lfs.\n", gradComputationTime);

            // 2. decrypt weight and bias
            Proto.resetTime();
            std::vector<mpz_class> decrypt_gradWeight_B(row * column);
            std::vector<mpz_class> decrypt_gradBias_B(row);
            for (int i = 0; i < row * column; ++i)
            {
                Proto.Decryption(decrypt_gradWeight_B[i], encrypt_gradWeight_B[i]);
            }
            
            for (int i = 0; i < row; ++i)
            {
                Proto.Decryption(decrypt_gradBias_B[i], encrypt_gradBias_B[i]);
            }
            printf("A解密梯度花费%lfs.\n", Proto.decryptionTime);

            // 3. weight - R || bias - R
            std::vector<vector<double>> gradWeightB(row, vector<double>(column));
            std::vector<double> gradBiasB(row, 0);
            for (int i = 0; i < row; ++i)
            {
                for (int j = 0; j < column; ++j)
                {
                    decrypt_gradWeight_B[i] = decrypt_gradWeight_B[i] - R;
                    Proto.Decode(gradWeightB[i][j], decrypt_gradWeight_B[i * column + j], true, 1e6);
                    gradWeightB[i][j] /= batchSize;
                }
            }
            for (int i = 0; i < row; ++i)
            {
                decrypt_gradBias_B[i] = decrypt_gradBias_B[i] - R;
                Proto.Decode(gradBiasB[i], decrypt_gradBias_B[i], false, 1e6);
                gradBiasB[i] /= batchSize;
            }

            // 更新模型A
            gradComputation::mySGDUpdateWeight(weightA, gradWeightA, lr);
            gradComputation::mySGDUpdateBias(biasA, gradBiasA, lr);

            // 更新模型B
            gradComputation::mySGDUpdateWeight(weightB, gradWeightB, lr);
            gradComputation::mySGDUpdateBias(biasB, gradBiasB, lr);

            // Output the loss and checkpoint every 100 batches.
            if ((batchIndex + 1) % 2 == 0)
            {
                sfile << "Epoch: " << epoch << " | Batch: " << batchIndex + 1
                      << " | Average Loss: " << lossSum / (batchIndex + 1) << std::endl;
            }
            auto iter_et = clock();
            cout << "iteration time = " << 1.0 * (iter_et - iter_st) / CLOCKS_PER_SEC << "s" << endl;
        }
    }
    sfile.close();
    return 0;
}
