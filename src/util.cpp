#include <iostream>
#include <gmpxx.h>
#include <vector>
#include <string>
#include <ctime>
using namespace std;

namespace myGenRand
{
    mpz_class GenRandomPrime(mp_bitcnt_t bits)
    {
        clock_t time = clock();
        mpz_class rand_num;
        gmp_randclass rc(gmp_randinit_default);
        rc.seed(time);
        rand_num = rc.get_z_bits(bits);
        mpz_setbit(rand_num.get_mpz_t(), bits - 1);
        mpz_nextprime(rand_num.get_mpz_t(), rand_num.get_mpz_t());
        return rand_num;
    }

    mpz_class GenRandomNumber(mp_bitcnt_t bits)
    {
        clock_t time = clock();
        mpz_class rand_num;
        gmp_randclass rc(gmp_randinit_default);
        rc.seed(time);
        rand_num = rc.get_z_bits(bits);
        mpz_setbit(rand_num.get_mpz_t(), bits - 1);
        return rand_num;
    }

    unsigned genRandUnsinged()
    {
        auto time = clock();
        srand(time);
        return rand();
    }

    // 生成一个小于N的随机数
    mpz_class GenRandomNumber_N(mp_bitcnt_t bits, mpz_class &N)
    {
        mpz_class res;
        do
        {
            res = GenRandomNumber(bits);
        } while (res >= N);
        return res;
    }

    mpz_class GenRandomOdd(mp_bitcnt_t bits)
    {
        mpz_class res = GenRandomNumber(bits);
        mpz_setbit(res.get_mpz_t(), 0);
        return res;
    }

    mpz_class GenRandom_Range(mp_bitcnt_t lower, mp_bitcnt_t upper)
    {
        auto time = clock();
        srand(time);
        auto bitWidth = (rand() % (upper - lower + 1)) + lower;
        return GenRandomNumber(bitWidth);
    }
}

namespace myFun
{
    mpz_class L_function(mpz_class &x, mpz_class &n)
    {
        mpz_class res;
        res = (x - 1) / n;
        return res;
    }

    mpz_class H_function(mpz_class &g, mpz_class &p, mpz_class &psquare)
    {
        mpz_class res;
        mpz_class x;
        mpz_class p_1 = p - 1;
        mpz_powm(x.get_mpz_t(), g.get_mpz_t(), p_1.get_mpz_t(), psquare.get_mpz_t());
        mpz_class L = L_function(x, p);
        mpz_invert(res.get_mpz_t(), L.get_mpz_t(), p.get_mpz_t());
        return res;
    }

    void quickExp(mpz_class &res, mpz_class &base, mpz_class &exponent, mpz_class &n)
    {
        res = 1;
        string str = exponent.get_str(2);
        for (int i = 0; i < str.size(); ++i)
        {
            mpz_powm_ui(res.get_mpz_t(), res.get_mpz_t(), 2, n.get_mpz_t());
            if (str[i] == '1')
            {
                res = res * base % n;
            }
        }
    }

}

void preCompute(vector<mpz_class> &v_pre, mpz_class &g, mpz_class &n)
{
    auto size = v_pre.size();
    mpz_class x = g;
    for (int i = 0; i < size; ++i)
    {
        v_pre[i] = x;
        x = x * x % n;
    }
}

