#include <iostream>
#include <vector>
#include <time.h>
#include "util.h"
#include "SBR_MOF.h"
using namespace std;

namespace SBR_MOF
{
    void replace(limbType *&newLimb, limbType *oriLimb, const int size)
    {
        // size of newLimb = size of oriLimb + 1
        newLimb = new limbType[size + 1];
        newLimb[size] = oriLimb[size - 1];
        for (int i = size - 1; i >= 1; --i)
        {
            newLimb[i] = oriLimb[i - 1] - oriLimb[i];
        }
        newLimb[0] = -1 * oriLimb[0];
        // substitutePattern(newLimb, size + 1);
        substitutePattern3(newLimb, size+1);
    }

    void substitutePattern(limbType *limb, const int size)
    {

        int i = size - 1;
        while (i >= 1)
        {
            // 1 | -1
            if (limb[i] == 1 && limb[i - 1] == -1)
            {
                limb[i] = 0;
                limb[i - 1] = 1;
            }
            // -1 | 1
            // else if (limb[i] == -1 && limb[i - 1] == 1)
            // {
            //     if (i == 1 || limb[i - 2] == 0)
            //     {
            //         limb[i] = 0;
            //         limb[i - 1] = -1;
            //     }
            // }
            else if (limb[i] == -1 && limb[i - 1] == 1)
            {
                limb[i] = 0;
                limb[i - 1] = -1;
            }
            --i;
        }
    }

    void substitutePattern3(limbType *limb, const int size)
    {

        int i = size - 1;
        while (i >= 2)
        {
            // 0 1 -1
            if (limb[i] == 0 && limb[i - 1] == 1 && limb[i - 2] == -1)
            {
                limb[i] = 0, limb[i - 1] = 0, limb[i - 2] = 1;
            }
            // 0 -1 1
            if (limb[i] == 0 && limb[i - 1] == -1 && limb[i - 2] == 1)
            {
                if (i - 3 >= 0 && limb[i - 3] == 0)
                {
                    limb[i] = 0, limb[i - 1] = 0, limb[i - 2] = -1;
                }
            }
            // 1 0 -1 | 1 -1 1
            if (limb[i] == 1 && limb[i - 1] == 0 && limb[i - 2] == -1)
            {
                if (i - 3 >= 0 && limb[i - 3] == 0)
                {
                    limb[i] = 0, limb[i - 1] = 1, limb[i - 2] = 1;
                }
            }
            if (limb[i] == 1 && limb[i - 1] == -1 && limb[i - 2] == 1)
            {
                limb[i] = 0, limb[i - 1] = 1, limb[i - 2] = 1;
            }
            // -1 1 -1 | 1 1 -1 额外情况
            if (limb[i] == -1 && limb[i - 1] == 1 && limb[i - 2] == -1)
            {
                limb[i] = -1, limb[i - 1] = 0, limb[i - 2] = 1;
            }
            if (limb[i] == 1 && limb[i - 1] == 1 && limb[i - 2] == -1)
            {
                limb[i] = 1, limb[i - 1] = 0, limb[i - 2] = 1;
            }
            --i;
        }
    }

    void _powm_sbr(mpz_class &res, mpz_class &base, mpz_class &baseInvert, limbType *exponent, int size, mpz_class &mod)
    {
        res = 1;
        for (int i = size; i >= 0; --i)
        {
            res = res * res % mod;

            if (exponent[i] == 1)
            {
                res = res * base % mod; // degradation的原因在于，先乘再模而不是蒙哥马利模乘
            }
            else if (exponent[i] == -1)
            {
                res = res * baseInvert % mod;
            }
        }
    }

    void _powm_pre(mpz_class &res, vector<mpz_class> &base, limbType *exponent, int size, mpz_class &mod)
    {
        mpz_class x = 1;
        for (int i = 0; i < size; ++i)
        {
            if (exponent[i] == 1)
            {
                x = x * base[i] % mod;
            }
        }
        res = x;
    }

    void _powm_pre_sbr(mpz_class &res, vector<mpz_class> &base, vector<mpz_class> &baseInvert, limbType *exponent, int size, mpz_class &mod)
    {
        mpz_class x = 1;
        for (int i = 0; i <= size; ++i)
        {
            if (exponent[i] == 1)
            {
                x = x * base[i] % mod;
            }
            else if (exponent[i] == -1)
            {
                x = x * baseInvert[i] % mod;
            }
        }
        res = x;
    }
};