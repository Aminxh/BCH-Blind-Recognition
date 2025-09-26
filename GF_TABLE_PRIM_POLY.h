#ifndef GF_TABLE_PRIM_POLY_H
#define GF_TABLE_PRIM_POLY_H

#include <armadillo>
#include <iostream>

using namespace std;
using namespace arma;

struct GF_PRIM_POLY
{
    rowvec r1 = {};
    rowvec r2 = {7};
    rowvec r3 = {11, 13};
    rowvec r4 = {19, 25};
    rowvec r5 = {37};
    rowvec r6 = {67};
    rowvec r7 = {137};
    rowvec r8 = {285};
    rowvec r9 = {529};
    rowvec r10 = {1033};
    rowvec r11 = {2053};
    rowvec r12 = {4179};
    rowvec r13 = {8219};
    rowvec r14 = {17475};
    rowvec r15 = {32771};
    rowvec r16 = {69643};


    vector<rowvec> prim_poly = {r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15, r16};
};

#endif // GF_TABLE_PRIM_POLY_H
