#pragma once
#include <iostream>
#include <algorithm>
#include <ctime>
#include <gmpxx.h>
#include <vector>
#include "util.h"

using namespace std;


class BASE
{
private:
    mp_bitcnt_t bits;
    mpz_class  n, g, lambda, mu, nsquare, PrimeP, PrimeQ, gInvert, gn, gnInvert;
    mpz_class h_p, h_q, p, q, psquare, qsquare, p_1, q_1, pInvertq;
    mpz_class r, rInvert;

    // SBR
    limbType *oriLimb, *newLimb;
    size_t n_bits;

    void GenKey();

public:
    // For Test
    double preProcessTime, encryptionTime, decryptionTime;
    int oriHamWeight, newHamWeight, minusOne;

    BASE(mp_bitcnt_t _bits);
    void Encryption_QEXP(mpz_class &c, mpz_class &m);
    void Encryption(mpz_class &c, mpz_class &m);
    void Encryption_SBR(mpz_class &c, mpz_class &m);
    void Decryption(mpz_class &m, mpz_class &c);
    void Encode(mpz_class &res, double scalar, const unsigned scale);
    void Decode(double &res, mpz_class plain, bool isMul, int scale_factor);
    void EncryptMul(mpz_class &res, mpz_class &c, mpz_class &m);
    void EncryptAdd(mpz_class &res, mpz_class &c1, mpz_class &c2);

    void resetTime();
};