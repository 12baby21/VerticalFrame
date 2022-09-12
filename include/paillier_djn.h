#pragma once
#include <iostream>
#include <gmpxx.h>
#include "util.h"
using namespace std;

// 可以继承一个有编码和解码的父类
class DJN
{
private:
    mp_bitcnt_t bits;
    mpz_class p, q, psquare, qsquare, p_1, q_1, pInvertq, gInvert;
    mpz_class n, g, nsquare, gn, gnInvert, h, h_s, h_sInvert;
    mpz_class x, xsquare, _xsquare; // Exponent
    mpz_class h_p, h_q;             // CRT

    // Precomputation
    vector<mpz_class> v_h_s;
    vector<mpz_class> v_h_sInvert;

    // SBR
    limbType *oriLimb, *newLimb;

    size_t n_bits, x_bits, a_bits;
    void GenKey();

public:
    // For Test
    double preProcessTime, encryptionTime, decryptionTime;
    int oriHamWeight, newHamWeight, minusOne;

    DJN(mp_bitcnt_t _bits);
    void Encryption(mpz_class &c, mpz_class &m);
    void Encryption_QEXP(mpz_class &c, mpz_class &m);
    void Encryption_SBR(mpz_class &c, mpz_class &m);
    void Encryption_pre(mpz_class &c, mpz_class &m);
    void Encryption_pre_SBR(mpz_class &c, mpz_class &m);
    void Decryption(mpz_class &m, mpz_class &c);
    void Decryption_QEXP(mpz_class &m, mpz_class &c);
    void Encode(mpz_class &res, double scalar, const unsigned scale);
    void Decode(double &res, mpz_class plain, bool isMul, int scale_factor);
    void EncryptMul(mpz_class &res, mpz_class &c, mpz_class &m);
    void EncryptAdd(mpz_class &res, mpz_class &c1, mpz_class &c2);
    void resetTime();
};