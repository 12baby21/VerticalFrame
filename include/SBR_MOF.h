#pragma once
#include <iostream>
#include <vector>


using namespace std;

namespace SBR_MOF
{
    void replace(limbType *&newLimb, limbType *oriLimb, const int size);

    void substitutePattern(limbType *limb, const int size);

    void substitutePattern3(limbType *limb, const int size);

    void _powm_sbr(mpz_class& res, mpz_class& base, mpz_class& baseInvert, limbType *exponent, int size, mpz_class& mod);

    void _powm_pre(mpz_class &res, vector<mpz_class> &base, limbType *exponent, int size, mpz_class &mod);

    void _powm_pre_sbr(mpz_class& res, vector<mpz_class>& base, vector<mpz_class>& baseInvert, limbType *exponent, int size, mpz_class& mod);
};