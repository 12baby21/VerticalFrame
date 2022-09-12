#include <iostream>
#include <gmpxx.h>
#include "util.h"
#include "readLimb.h"
#include "SBR_MOF.h"
#include "paillier_djn.h"

using namespace std;

DJN::DJN(mp_bitcnt_t _bits)
{
    bits = _bits;
    preProcessTime = 0;
    encryptionTime = 0;
    decryptionTime = 0;
    oriHamWeight = newHamWeight = minusOne = 0;
    GenKey();
    v_h_s = vector<mpz_class>(a_bits + 1);
    v_h_sInvert = vector<mpz_class>(a_bits);
    preCompute(v_h_s, h_s, nsquare);
    preCompute(v_h_sInvert, h_sInvert, nsquare);
}

void DJN::GenKey()
{
    p = myGenRand::GenRandomPrime(bits);
    q = myGenRand::GenRandomPrime(bits);
    psquare = p * p;
    qsquare = q * q;
    p_1 = p - 1;
    q_1 = q - 1;
    mpz_invert(pInvertq.get_mpz_t(), p.get_mpz_t(), q.get_mpz_t());

    n = p * q;
    g = n + 1;
    nsquare = n * n;
    n_bits = mpz_sizeinbase(n.get_mpz_t(), 2);
    x_bits = n_bits;
    a_bits = (n_bits + 1) / 2;

    mpz_class x;
    x = myGenRand::GenRandomNumber_N(x_bits, n);
    xsquare = x * x;
    _xsquare = -1 * xsquare;
    mpz_mod(h.get_mpz_t(), _xsquare.get_mpz_t(), nsquare.get_mpz_t());
    mpz_powm(h_s.get_mpz_t(), h.get_mpz_t(), n.get_mpz_t(), nsquare.get_mpz_t());
    mpz_invert(h_sInvert.get_mpz_t(), h_s.get_mpz_t(), nsquare.get_mpz_t());

    // preComputation for CRT
    h_p = myFun::H_function(g, p, psquare);
    h_q = myFun::H_function(g, q, qsquare);
}

void DJN::Encryption(mpz_class &c, mpz_class &m)
{
    mpz_class a = myGenRand::GenRandomNumber(a_bits);
    mpz_class gm, rn;
    gm = 1 + m * n;
    mpz_powm(rn.get_mpz_t(), h_s.get_mpz_t(), a.get_mpz_t(), nsquare.get_mpz_t());
    c = gm * rn % nsquare;
}

void DJN::Encryption_QEXP(mpz_class &c, mpz_class &m)
{
    auto st = clock();
    mpz_class a = myGenRand::GenRandomNumber(a_bits);
    mpz_class gm, rn;
    gm = 1 + m * n;
    myFun::quickExp(rn, h_s, a, nsquare);
    c = gm * rn % nsquare;
    auto et = clock();
    encryptionTime += 1.0 * (et - st) / CLOCKS_PER_SEC;
}

void DJN::Encryption_SBR(mpz_class &c, mpz_class &m)
{
    // SBR Generation
    auto st = clock();
    mpz_class a = myGenRand::GenRandomNumber(a_bits);
    readLimbStr(oriLimb, a_bits, a);
    SBR_MOF::replace(newLimb, oriLimb, a_bits);

#ifdef TEST_HAMMING_WEIGHT
    oriHamWeight += mpz_popcount(a.get_mpz_t());
    for (int i = 0; i <= a_bits; ++i)
    {
        if (newLimb[i] != 0)
        {
            ++hamming_new;
        }
        if (newLimb[i] == -1)
        {
            ++minusone;
        }
    }
#endif

    auto et = clock();
    preProcessTime += 1.0 * (et - st) / CLOCKS_PER_SEC;

    // 加密
    st = clock();
    mpz_class gm, rn;
    gm = (1 + m * n) % nsquare;
    SBR_MOF::_powm_sbr(rn, h_s, h_sInvert, newLimb, a_bits, nsquare);
    c = gm * rn % nsquare;
    et = clock();
    encryptionTime += 1.0 * (et - st) / CLOCKS_PER_SEC;
}

void DJN::Encryption_pre(mpz_class &c, mpz_class &m)
{
    // 预处理
    auto st = clock();
    mpz_class a = myGenRand::GenRandomNumber(a_bits);
    readLimbStr(oriLimb, a_bits, a);
    auto et = clock();
    preProcessTime += 1.0 * (et - st) / CLOCKS_PER_SEC;

    // 加密
    st = clock();
    mpz_class gm, rn;
    gm = (1 + m * n) % nsquare;
    SBR_MOF::_powm_pre(rn, v_h_s, oriLimb, a_bits, nsquare);
    c = gm * rn % nsquare;
    et = clock();
    encryptionTime += 1.0 * (et - st) / CLOCKS_PER_SEC;
}

void DJN::Encryption_pre_SBR(mpz_class &c, mpz_class &m)
{
    // 预处理
    auto st = clock();
    mpz_class a = myGenRand::GenRandomNumber(a_bits);
    readLimbStr(oriLimb, a_bits, a);
    SBR_MOF::replace(newLimb, oriLimb, a_bits);
    auto et = clock();
    preProcessTime += 1.0 * (et - st) / CLOCKS_PER_SEC;

    // 加密
    st = clock();
    mpz_class gm, rn;
    gm = 1 + m * n;
    SBR_MOF::_powm_pre_sbr(rn, v_h_s, v_h_sInvert, newLimb, a_bits, nsquare);
    c = gm * rn % nsquare;
    et = clock();
    encryptionTime += 1.0 * (et - st) / CLOCKS_PER_SEC;
}

void DJN::Decryption(mpz_class &m, mpz_class &c)
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
}

void DJN::Encode(mpz_class &res, double scalar, const unsigned scale = 1e6)
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

void DJN::Decode(double &res, mpz_class plain, bool isMul, int scale_factor = 1e6)
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

void DJN::EncryptMul(mpz_class &res, mpz_class &c, mpz_class &m) 
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

void DJN::EncryptAdd(mpz_class &res, mpz_class &c1, mpz_class &c2)
{
    res = (c1 * c2) % nsquare;
}

void DJN::resetTime()
{
    encryptionTime = decryptionTime = preProcessTime = 0;
}