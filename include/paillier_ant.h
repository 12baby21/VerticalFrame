#pragma once
#include <iostream>
#include <gmpxx.h>
#include "util.h"

using namespace std;

class ANT
{
private:
    mp_bitcnt_t n_length;
    mpz_class p, q, P, Q, p1, q1;
    mpz_class alpha, b, h, hn, hnInvert;
    mpz_class nsquare, n1, a2, a2Invert;
    unsigned ll;

    // Precomputation Table
    vector<mpz_class> v_hn;
    vector<mpz_class> v_hnInvert;

    void GenKey();

public:
    mpz_class n;
    // For Test
    double preProcessTime, encryptionTime, decryptionTime;
    int oriHamWeight, newHamWeight, minusOne;

    ANT(mp_bitcnt_t _n_length);
    void Encryption_QEXP(mpz_class & c, mpz_class & m);
    void Encryption_SBR(mpz_class & c, mpz_class & m);
    void Encryption_pre(mpz_class & c, mpz_class & m);
    void Encryption_pre_SBR(mpz_class & c, mpz_class & m);
    void Decryption(mpz_class & m, mpz_class & c);
    void Encode(mpz_class &res, double scalar, const unsigned scale);
    void Decode(double &res, mpz_class plain, bool isMul, int scale_factor);
    void EncryptMul(mpz_class &res, mpz_class &c, mpz_class &m);
    void EncryptAdd(mpz_class &res, mpz_class &c1, mpz_class &c2);
    void resetTime();
};
