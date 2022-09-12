#include <iostream>
#include <algorithm>
#include <ctime>
#include <gmpxx.h>
#include <vector>
#include "util.h"
#include "readLimb.h"
#include "SBR_MOF.h"
#include "paillier_base.h"

using namespace std;

BASE::BASE(mp_bitcnt_t _bits)
{
    bits = _bits;
    preProcessTime = 0;
    encryptionTime = 0;
    decryptionTime = 0;
    oriHamWeight = newHamWeight = minusOne = 0;
    GenKey();
}

void BASE::GenKey()
{
    p = myGenRand::GenRandomPrime(bits);
    q = myGenRand::GenRandomPrime(bits);
    psquare = p * p;
    qsquare = q * q;
    p_1 = p - 1;
    q_1 = q - 1;
    mpz_invert(pInvertq.get_mpz_t(), p.get_mpz_t(), q.get_mpz_t());

    n = p * q;
    nsquare = n * n;
    n_bits = mpz_sizeinbase(n.get_mpz_t(), 2);
    
    g = n + 1;
    lambda = p_1 * q_1;
    mu = lambda - 1;
    mpz_powm(gn.get_mpz_t(), g.get_mpz_t(), n.get_mpz_t(), nsquare.get_mpz_t());
    mpz_invert(gnInvert.get_mpz_t(), gn.get_mpz_t(), nsquare.get_mpz_t());

    readLimbStr(oriLimb, n_bits, n);
    SBR_MOF::replace(newLimb, oriLimb, n_bits);
    mpz_invert(gInvert.get_mpz_t(), g.get_mpz_t(), nsquare.get_mpz_t());

    // preComputation for CRT
    h_p = myFun::H_function(g, p, psquare);
    h_q = myFun::H_function(g, q, qsquare);

    // Obfrustration
    r = myGenRand::GenRandomNumber_N(n_bits, n);
    rInvert;
    mpz_invert(rInvert.get_mpz_t(), r.get_mpz_t(), nsquare.get_mpz_t());
}

void BASE::Encryption_QEXP(mpz_class &c, mpz_class &m)
{
    // Random Key
    auto st = clock();
    mpz_class gm;
    mpz_class rn;
    myFun::quickExp(gm, g, m, nsquare);
    myFun::quickExp(rn, r, n, nsquare);

    c = gm * rn % nsquare;
    auto et = clock();
    encryptionTime += 1.0 * (et - st) / CLOCKS_PER_SEC;
}

void BASE::Encryption(mpz_class &c, mpz_class &m)
{
    auto st = clock();
    // Random Key
    mpz_class r = myGenRand::GenRandomNumber_N(n_bits, n);

    // gm = g^m mod n^2
    // rn = r^n mod n^2
    mpz_class gm;
    mpz_class rn;

    mpz_powm(gm.get_mpz_t(), g.get_mpz_t(), m.get_mpz_t(), nsquare.get_mpz_t());
    mpz_powm(rn.get_mpz_t(), r.get_mpz_t(), n.get_mpz_t(), nsquare.get_mpz_t());

    c = gm * rn % nsquare;
    auto et = clock();
    encryptionTime += 1.0 * (et - st) / CLOCKS_PER_SEC;
}

void BASE::Encryption_SBR(mpz_class &c, mpz_class &m)
{
#ifdef TEST_HAMMING_WEIGHT
    oriHamWeight += mpz_popcount(n.get_mpz_t());
    for (int i = 0; i <= n_bits; ++i)
    {
        if (newLimb[i] != 0)
        {
            ++newHamWeight;
        }
        if (newLimb[i] == -1)
        {
            ++minusOne;
        }
    }
#endif

    // gm = g^m mod n^2
    // rn = r^n mod n^2
    auto st = clock();
    mpz_class gm;
    mpz_class rn;

    myFun::quickExp(gm, g, m, nsquare);
    SBR_MOF::_powm_sbr(rn, r, rInvert, newLimb, n_bits, nsquare);

    c = gm * rn % nsquare;

    auto et = clock();
    encryptionTime += 1.0 * (et - st) / CLOCKS_PER_SEC;
}

void BASE::Decryption(mpz_class &m, mpz_class &c)
{
    auto st = clock();
    mpz_class a, b, u;
    mpz_class tmp_p, tmp_q;
    mpz_powm(tmp_p.get_mpz_t(), c.get_mpz_t(), p_1.get_mpz_t(), psquare.get_mpz_t());
    mpz_powm(tmp_q.get_mpz_t(), c.get_mpz_t(), q_1.get_mpz_t(), qsquare.get_mpz_t());
    a = myFun::L_function(tmp_p, p);
    b = myFun::L_function(tmp_q, q);
    a = a * h_p % p;
    b = b * h_q % q;
    u = (a - b) * pInvertq % q;
    m = b + u * p;
    auto et = clock();
    decryptionTime += 1.0 * (et - st) / CLOCKS_PER_SEC;
    // mpz_class L;
    // mpz_powm(L.get_mpz_t(), c.get_mpz_t(), lambda.get_mpz_t(), nsquare.get_mpz_t());
    // L = (L - 1) / n;

    // mpz_class lambdainvert;
    // mpz_invert(lambdainvert.get_mpz_t(), lambda.get_mpz_t(), n.get_mpz_t());
    // L = L * lambdainvert;
    // L = L % n;
    // // L = L * mu
    // m = L;
}

void BASE::Encode(mpz_class &res, double scalar, const unsigned scale = 1e6)
{
    bool flag = (scalar < 0);
    if (flag)
        scalar = -scalar;
    unsigned after_scale = static_cast<unsigned>(scalar * scale);
    if (flag)
    {
        res = n - after_scale;
    }
    else
    {
        res = after_scale;
    }
}

void BASE::Decode(double &res, mpz_class plain, bool isMul, int scale_factor = 1e6)
{
    long ret;
    mpz_class max_int = n / 3;           // n/3
    mpz_class forNegative = max_int * 2; // 2n/3
    int isPositive = cmp(max_int, plain) > 0;
    int isNegative = cmp(plain, forNegative) > 0;

    if (!isMul)
    {
        if (isNegative == 1)
        {
            mpz_class tmp = n - plain;
            ret = tmp.get_ui();
            ret = -ret;
        }
        else if (isPositive == 1)
        {
            ret = plain.get_si();
        }
        else
        {
            cout << "There is a possible OVERFLOW!\n";
            exit(-1);
        }
    }
    else
    {
        if (isNegative == 1)
            plain = n - plain;
        plain = plain / scale_factor;
        ret = plain.get_si();
        if (isNegative == 1)
            ret = -ret;
    }

    res = static_cast<double>(ret) / scale_factor;
}

void BASE::EncryptMul(mpz_class &res, mpz_class &c, mpz_class &m)
{
    mpz_class max_int = n / 3;           // n/3
    mpz_class forNegative = max_int * 2; // 2n/3
    if (cmp(m, forNegative) == 1)
    {
        mpz_class neg_c, neg_scalar;
        mpz_invert(neg_c.get_mpz_t(), c.get_mpz_t(), nsquare.get_mpz_t());
        neg_scalar = n - m;
        mpz_powm(res.get_mpz_t(), neg_c.get_mpz_t(), neg_scalar.get_mpz_t(), nsquare.get_mpz_t());
        return;
    }

    mpz_powm(res.get_mpz_t(), c.get_mpz_t(), m.get_mpz_t(), nsquare.get_mpz_t());
}

void BASE::EncryptAdd(mpz_class &res, mpz_class &c1, mpz_class &c2)
{
    res = (c1 * c2) % nsquare;
}

void BASE::resetTime()
{
    encryptionTime = decryptionTime = preProcessTime = 0;
}