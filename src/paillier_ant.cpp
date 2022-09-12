#include <iostream>
#include <gmpxx.h>
#include <ctime>
#include "util.h"
#include "readLimb.h"
#include "SBR_MOF.h"
#include "paillier_ant.h"

using namespace std;

ANT::ANT(mp_bitcnt_t _n_length)
{
    n_length = _n_length;
    preProcessTime = 0;
    encryptionTime = 0;
    decryptionTime = 0;
    oriHamWeight = newHamWeight = minusOne = 0;
    GenKey();

    // Precomputation
    v_hn = vector<mpz_class>(ll + 1);
    v_hnInvert = vector<mpz_class>(ll + 1);
    preCompute(v_hn, hn, nsquare);
    preCompute(v_hnInvert, hn, nsquare);
}

void ANT::GenKey()
{
    switch (n_length)
    {
    case 2048:
        ll = 448;
        break;
    case 3072:
        ll = 512;
        break;
    case 7680:
        ll = 768;
        break;
    default:
    {
        if (n_length < 2048)
            ll = 448;
        else if (n_length > 7680)
            ll = 768;
        break;
    }
    }
    p = 4;
    q = 4;
    n = 4;
    P = 4;
    Q = 4;
    p1 = 4;
    q1 = 4;

    unsigned n_len = 0;
    int i = 0;
    bool prim = false;
    mpz_class gcd1, gcd2, gcd3, gcd4;
    mpz_gcd(gcd1.get_mpz_t(), p.get_mpz_t(), p1.get_mpz_t());
    mpz_gcd(gcd2.get_mpz_t(), p.get_mpz_t(), q1.get_mpz_t());
    mpz_gcd(gcd3.get_mpz_t(), q.get_mpz_t(), p1.get_mpz_t());
    mpz_gcd(gcd4.get_mpz_t(), q.get_mpz_t(), q1.get_mpz_t());

    while (!(prim && gcd1 == 1 && gcd2 == 1 && gcd3 == 1 && gcd4 == 1))
    {
        p = myGenRand::GenRandomPrime(ll / 2);
        q = p;
        while (q == p)
        {
            q = myGenRand::GenRandomPrime(ll / 2);
        }
        p1 = myGenRand::GenRandomOdd((n_length - ll) / 2 - 1);
        q1 = p1;
        while (q1 == p1)
        {
            q1 = myGenRand::GenRandomOdd((n_length - ll) / 2 - 1);
        }
        P = 2 * p * p1 + 1;
        Q = 2 * q * q1 + 1;
        i = 0;

        while (mpz_probab_prime_p(P.get_mpz_t(), 64) == 0 ||
               mpz_sizeinbase(P.get_mpz_t(), 2) != n_length / 2)
        {
            if (i > 1000)
            {
                i = 0;
                p = myGenRand::GenRandomOdd(ll / 2);
            }
            p1 = myGenRand::GenRandomOdd((n_length - ll) / 2 - 1);
            P = 2 * p * p1 + 1;
            i += 1;
        }
        i = 0;
        n = 0;

        while (mpz_probab_prime_p(Q.get_mpz_t(), 64) == 0 ||
               mpz_sizeinbase(Q.get_mpz_t(), 2) != n_length / 2)
        {
            if (i > 1000)
            {
                i = 0;
                q = myGenRand::GenRandomOdd(ll / 2);
                while (q == p)
                {
                    q = myGenRand::GenRandomOdd(ll / 2);
                }
            }
            q1 = myGenRand::GenRandomOdd((n_length - ll) / 2 - 1);
            while (q1 == p1)
            {
                q1 = myGenRand::GenRandomOdd((n_length - ll) / 2 - 1);
            }
            Q = 2 * q * q1 + 1;
            i += 1;
        }

        n = P * Q;
        n_len = mpz_sizeinbase(n.get_mpz_t(), 2);
        int isPrimeP = mpz_probab_prime_p(P.get_mpz_t(), 64);
        int isPrimeQ = mpz_probab_prime_p(Q.get_mpz_t(), 64);
        prim = (isPrimeP >= 1 && isPrimeQ >= 1);

        mpz_gcd(gcd1.get_mpz_t(), p.get_mpz_t(), p1.get_mpz_t());
        mpz_gcd(gcd2.get_mpz_t(), p.get_mpz_t(), q1.get_mpz_t());
        mpz_gcd(gcd3.get_mpz_t(), q.get_mpz_t(), p1.get_mpz_t());
        mpz_gcd(gcd4.get_mpz_t(), q.get_mpz_t(), q1.get_mpz_t());
    }
    alpha = p * q;
    b = (P - 1) * (Q - 1) / (4 * q * p);
    mpz_class y = myGenRand::GenRandomNumber(n_length - 1);
    mpz_class b2 = 2 * b;
    mpz_powm(h.get_mpz_t(), y.get_mpz_t(), b2.get_mpz_t(), n.get_mpz_t());
    h = -h;
    mpz_mod(h.get_mpz_t(), h.get_mpz_t(), n.get_mpz_t());

    // 一些其他初始化
    nsquare = n * n;
    n1 = n + 1;
    a2 = alpha * 2;
    mpz_invert(a2Invert.get_mpz_t(), a2.get_mpz_t(), n.get_mpz_t());
    mpz_powm(hn.get_mpz_t(), h.get_mpz_t(), n.get_mpz_t(), nsquare.get_mpz_t()); // h^n mod n^2
    mpz_invert(hnInvert.get_mpz_t(), hn.get_mpz_t(), nsquare.get_mpz_t());
}

void ANT::Encryption_QEXP(mpz_class &c, mpz_class &m)
{
    auto st = clock();
    mpz_class gm, rn;
    // mpz_class r = myGenRand::GenRandom_Range(1, ll);
    mpz_class r = myGenRand::GenRandomNumber(ll);
    gm = (1 + m * n) % nsquare;
    // mpz_powm(rn.get_mpz_t(), hn.get_mpz_t(), r.get_mpz_t(), nsquare.get_mpz_t());
    myFun::quickExp(rn, hn, r, nsquare);
    c = gm * rn % nsquare;
    auto et = clock();
    encryptionTime += 1.0 * (et - st) / CLOCKS_PER_SEC;
}

void ANT::Encryption_SBR(mpz_class &c, mpz_class &m)
{
    // mpz_class r = myGenRand::GenRandom_Range(1, ll);
    mpz_class r = myGenRand::GenRandomNumber(ll);

    // 预处理
    auto st = clock();
    auto size = mpz_sizeinbase(r.get_mpz_t(), 2);
    limbType *oriLimb, *newLimb;
    readLimbStr(oriLimb, size, r);
    SBR_MOF::replace(newLimb, oriLimb, size);

#ifdef TEST_HAMMING_WEIGHT
    oriHamWeight += mpz_popcount(r.get_mpz_t());
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
    gm = (1 + m * n) % nsquare;
    SBR_MOF::_powm_sbr(rn, hn, hnInvert, newLimb, size, nsquare);
    c = gm * rn % nsquare;
    et = clock();
    encryptionTime += 1.0 * (et - st) / CLOCKS_PER_SEC;
}

void ANT::Encryption_pre(mpz_class &c, mpz_class &m)
{
    // mpz_class r = myGenRand::GenRandom_Range(1, ll);
    mpz_class r = myGenRand::GenRandomNumber(ll);

    // preProcessing
    auto st = clock();
    auto size = mpz_sizeinbase(r.get_mpz_t(), 2);
    limbType *oriLimb, *newLimb;
    readLimbStr(oriLimb, size, r);
    auto et = clock();
    preProcessTime += 1.0 * (et - st) / CLOCKS_PER_SEC;

    // Encryption
    st = clock();
    mpz_class gm, rn;
    gm = (1 + m * n) % nsquare;
    SBR_MOF::_powm_pre(rn, v_hn, oriLimb, size, nsquare);
    c = gm * rn % nsquare;
    et = clock();
    encryptionTime += 1.0 * (et - st) / CLOCKS_PER_SEC;
}

void ANT::Encryption_pre_SBR(mpz_class &c, mpz_class &m)
{
    mpz_class r = myGenRand::GenRandomNumber(ll);

    // 预处理
    auto st = clock();
    auto size = mpz_sizeinbase(r.get_mpz_t(), 2);
    limbType *oriLimb, *newLimb;
    readLimbStr(oriLimb, size, r);
    SBR_MOF::replace(newLimb, oriLimb, size);
    auto et = clock();
    preProcessTime += 1.0 * (et - st) / CLOCKS_PER_SEC;

    // 加密
    st = clock();
    mpz_class gm, rn;
    gm = (1 + m * n) % nsquare;
    SBR_MOF::_powm_pre_sbr(rn, v_hn, v_hnInvert, newLimb, size, nsquare);
    c = gm * rn % nsquare;
    et = clock();
    encryptionTime += 1.0 * (et - st) / CLOCKS_PER_SEC;
}

void ANT::Decryption(mpz_class &m, mpz_class &c)
{
    auto st = clock();
    mpz_class ca;
    mpz_class L;

    mpz_powm(ca.get_mpz_t(), c.get_mpz_t(), a2.get_mpz_t(), nsquare.get_mpz_t());
    L = (ca - 1) / n;
    m = L * a2Invert % n;
    auto et = clock();
    decryptionTime += 1.0 * (et - st) / CLOCKS_PER_SEC;
}

void ANT::Encode(mpz_class &res, double scalar, const unsigned scale = 1e6)
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

void ANT::Decode(double &res, mpz_class plain, bool isMul, int scale_factor = 1e6)
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

void ANT::EncryptMul(mpz_class &res, mpz_class &c, mpz_class &m) 
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

void ANT::EncryptAdd(mpz_class &res, mpz_class &c1, mpz_class &c2)
{
    res = (c1 * c2) % nsquare;
}

void ANT::resetTime()
{
    encryptionTime = decryptionTime = preProcessTime = 0;
}