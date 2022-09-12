#pragma once
#include <gmp.h>
#include <gmpxx.h>
#include <iostream>
#include <vector>
#include <ctime>

using namespace std;

// Configs
using limbType = int;

// 训练配置
const int batchSize = 512;

// Generate Random Number
namespace myGenRand
{
    mpz_class GenRandomPrime(mp_bitcnt_t bits);

    mpz_class GenRandomNumber(mp_bitcnt_t bits);

    // 生成位宽范围在[lower, upper]范围内的随机数
    mpz_class GenRandom_Range(mp_bitcnt_t lower, mp_bitcnt_t upper);

    // 生成大小不超过N的位宽为bits的随机数
    mpz_class GenRandomNumber_N(mp_bitcnt_t bits, mpz_class &N);

    unsigned genRandUnsinged();

    // 生成位宽为bits的随机奇数
    mpz_class GenRandomOdd(mp_bitcnt_t bits);
}

// Useful Functions
namespace myFun
{
    mpz_class L_function(mpz_class &x, mpz_class &n);

    mpz_class H_function(mpz_class &g, mpz_class &p, mpz_class &psquare);

    void quickExp(mpz_class &res, mpz_class &base, mpz_class &exponent, mpz_class &n);

}

void preCompute(vector<mpz_class> &v_pre, mpz_class &g, mpz_class &n);