#ifndef BCH_BLIND_H
#define BCH_BLIND_H

#include <iostream>
#include <tuple>
#include <armadillo>
#include "GF_TABLE_PRIM_POLY.h"
#include "GF_TABLE1.h"
#include "GF_TABLE2.h"

using namespace arma;
using namespace std;

struct dataParameters
{
    int extended;
    int n_et;
    int shortend_et;
    int k_et;
    int t;
    int primPoly_et;
    string status;
};

struct object
{
    rowvec p_vec;
    int64_t m;
    int64_t prim_poly;
    int8_t is_prime;
    int64_t q;
    umat x;
    rowvec twos;
};


class BCH_Blind
{
public:
    BCH_Blind(const uchar_rowvec& X, size_t n, double threshold);
    dataParameters run_BCH_Blind();
    vec allPrimPoly(int m);
    object myGaliosField(uchar_mat& x, int64_t m, int64_t prim_poly);
    object myMPower(object A, double pwr);
    void getTables(object A);
    umat myPower(object x, rowvec y);
    int myGFMultiply(uint32_t a, uint32_t b, int64_t irr, rowvec twos, int64_t m, int64_t q, umat E2P, umat P2E);
    umat myGFMtimes(umat& x, umat& y, int64_t irr, rowvec twos, int64_t m, int64_t q, umat E2P, umat P2E);
    umat myGFtimes(uvec& x, uvec& y, int64_t irr, rowvec twos, int64_t m, int64_t q, umat E2P, umat P2E);
    umat myMTimes(uchar_mat& x, object y);
    mat myBCHNumerr(int n);

private:
    uchar_rowvec m_X_rowvec;                    // This is the encoded data which we take as input
    size_t m_n;                                 // This is code word length which we take as input
    double m_threshold;                         // This is the threshold data which we take as input
    uchar_mat m_X_mat;                          // We will reshape m_X_rowvec so we need this variable
    int64_t m_m;                                // Order of galios field
    int64_t m_n0;                               // Number of members in galios field
    int64_t GF_TABLE_M;                         // Galios field order that we will get from table
    int64_t GF_TABLE_PRIM_POLY;                 // Galios field primPoly that we will get from table
    umat GF_TABLE1;                             // Galios field 1 which we get from table
    umat GF_TABLE2;                             // Galios field 2 which we get from table
    struct GF_TABLE1 GF_MAIN_TABLE1;
    struct GF_TABLE2 GF_MAIN_TABLE2;
    struct GF_PRIM_POLY GF_MAIN_PRIM_POLY;
};

#endif // BCH_BLIND_H
