#include <iostream>
#include <gmpxx.h>
#include "util.h"

using namespace std;

extern int ham_ori, ham_sbr;

void readLimbStr(limbType *&res, size_t size, mpz_class &num)
{
    string str;
    str = num.get_str(2);
    res = new limbType[size];
    for (int i = 0; i < size; ++i)
    {
        res[size - i - 1] = str[i] - '0';
    }
}
