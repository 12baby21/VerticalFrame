#include <iostream>
#include <gmpxx.h>
#include "util.h"
#include "readLimb.h"
#include "SBR_MOF.h"
#include "paillier_opt.h"

using namespace std;

OPT::OPT(mp_bitcnt_t _bits)
{
    bits = _bits;
    preProcessTime = 0;
    encryptionTime = 0;
    decryptionTime = 0;
    oriHamWeight = newHamWeight = minusOne = 0;
    v_g = vector<mpz_class>(3073);
    v_gn = vector<mpz_class>(3073);
    v_gInvert = vector<mpz_class>(3073);
    v_gnInvert = vector<mpz_class>(3073);
    GenKey();
}

void OPT::GenKey()
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
    // g = myGenRand::GenRandom_Range(1, n_bits);
    g = n + 1;
    mpz_invert(gInvert.get_mpz_t(), g.get_mpz_t(), nsquare.get_mpz_t());
    mpz_powm(gn.get_mpz_t(), g.get_mpz_t(), n.get_mpz_t(), nsquare.get_mpz_t());
    mpz_invert(gnInvert.get_mpz_t(), gn.get_mpz_t(), nsquare.get_mpz_t());

    // preComputation Table
    preCompute(v_g, g, nsquare);
    preCompute(v_gInvert, gInvert, nsquare);
    preCompute(v_gn, gn, nsquare);
    preCompute(v_gnInvert, gnInvert, nsquare);

    // preComputation for CRT
    h_p = myFun::H_function(g, p, psquare);
    h_q = myFun::H_function(g, q, qsquare);
}

void OPT::Encryption_QEXP(mpz_class &c, mpz_class &m)
{
    auto st = clock();
    mpz_class r = myGenRand::GenRandomNumber(448);
    mpz_class e = r;

    mpz_class gm, rn;
    myFun::quickExp(gm, g, m, nsquare);
    myFun::quickExp(rn, gn, e, nsquare);
    c = gm * rn % nsquare;
    // mpz_powm(c.get_mpz_t(), g.get_mpz_t(), e.get_mpz_t(), nsquare.get_mpz_t());
    auto et = clock();
    encryptionTime += 1.0 * (et - st) / CLOCKS_PER_SEC;
}

void OPT::Encryption_SBR(mpz_class &c, mpz_class &m)
{
    mpz_class r = myGenRand::GenRandomNumber(448);
    mpz_class e = r;

    // 预处理
    auto st = clock();
    auto size = mpz_sizeinbase(e.get_mpz_t(), 2);
    limbType *oriLimb, *newLimb;
    readLimbStr(oriLimb, size, e);
    SBR_MOF::replace(newLimb, oriLimb, size);

#ifdef TEST_HAMMING_WEIGHT
    oriHamWeight += mpz_popcount(e.get_mpz_t());
    for (int i = 0; i <= size; ++i)
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

    auto et = clock();
    preProcessTime += 1.0 * (et - st) / CLOCKS_PER_SEC;

    // 加密
    st = clock();
    mpz_class gm, rn;
    myFun::quickExp(gm, g, m, nsquare);
    SBR_MOF::_powm_sbr(rn, gn, gnInvert, newLimb, size, nsquare);
    c = gm * rn % nsquare;
    et = clock();
    encryptionTime += 1.0 * (et - st) / CLOCKS_PER_SEC;
}

void OPT::Encryption_pre(mpz_class &c, mpz_class &m)
{
    mpz_class r = myGenRand::GenRandomNumber(448);
    mpz_class e = r;

    // 预处理
    auto st = clock();
    auto sizeE = mpz_sizeinbase(e.get_mpz_t(), 2);
    auto sizeM = mpz_sizeinbase(m.get_mpz_t(), 2);
    limbType *eLimb, *mLimb;
    readLimbStr(eLimb, sizeE, e);
    readLimbStr(mLimb, sizeM, m);
    auto et = clock();
    preProcessTime += 1.0 * (et - st) / CLOCKS_PER_SEC;

    // 加密
    st = clock();
    mpz_class gm, rn;
    mpz_powm(gm.get_mpz_t(), g.get_mpz_t(), m.get_mpz_t(), nsquare.get_mpz_t());
    // SBR_MOF::_powm_pre(gm, v_g, mLimb, sizeM, nsquare);
    SBR_MOF::_powm_pre(rn, v_gn, eLimb, sizeE, nsquare);
    c = gm * rn % nsquare;
    et = clock();
    encryptionTime += 1.0 * (et - st) / CLOCKS_PER_SEC;
}

void OPT::Encryption_pre_SBR(mpz_class &c, mpz_class &m)
{
    mpz_class r = myGenRand::GenRandomNumber(448);
    mpz_class e = r;

    // 预处理
    auto st = clock();
    auto sizeM = mpz_sizeinbase(m.get_mpz_t(), 2);
    auto sizeE = mpz_sizeinbase(e.get_mpz_t(), 2);

    limbType *oriMLimb, *newMLimb;
    limbType *oriELimb, *newELimb;
    readLimbStr(oriMLimb, sizeM, m);
    readLimbStr(oriELimb, sizeE, e);
    SBR_MOF::replace(newMLimb, oriMLimb, sizeM);
    SBR_MOF::replace(newELimb, oriELimb, sizeE);
    auto et = clock();
    preProcessTime += 1.0 * (et - st) / CLOCKS_PER_SEC;

    // 加密
    st = clock();
    mpz_class gm, rn;
    mpz_powm(gm.get_mpz_t(), g.get_mpz_t(), m.get_mpz_t(), nsquare.get_mpz_t());
    // SBR_MOF::_powm_pre_sbr(gm, v_g, v_gInvert, newMLimb, sizeM, nsquare);
    SBR_MOF::_powm_pre_sbr(rn, v_gn, v_gnInvert, newELimb, sizeE, nsquare);
    c = gm * rn % nsquare;
    et = clock();
    encryptionTime += 1.0 * (et - st) / CLOCKS_PER_SEC;
}

void OPT::Decryption(mpz_class &m, mpz_class &c)
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

void OPT::Encode(mpz_class &res, double scalar, const unsigned scale = 1e6)
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

void OPT::Decode(double &res, mpz_class plain, bool isMul, int scale_factor = 1e6)
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

void OPT::EncryptMul(mpz_class &res, mpz_class &c, mpz_class &m)
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

void OPT::EncryptAdd(mpz_class &res, mpz_class &c1, mpz_class &c2)
{
    res = (c1 * c2) % nsquare;
}


void OPT::resetTime()
{
    encryptionTime = decryptionTime = preProcessTime = 0;
}