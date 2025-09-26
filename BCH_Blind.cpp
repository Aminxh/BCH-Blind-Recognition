#include "BCH_Blind.h"


BCH_Blind::BCH_Blind(const uchar_rowvec &X, size_t n, double threshold)
{
    m_X_rowvec = X;
    m_n = n;
    m_threshold = threshold;
}

dataParameters BCH_Blind::run_BCH_Blind()
{
    dataParameters block;
    int shorten;
    vec k;
    double mlcr;
    int primPoly;

    if(m_n <= 0)        // Check to see n is a posetive number or not
    {
        cout << "invalid codeword length" << endl;

        block.extended = 0;         //
        block.shortend_et = -1;     //
        block.k_et = -1;            //
        block.n_et = -1;            //  Just set non-sense numbers
        block.primPoly_et = -1;     //
        block.t = -1;               //
        block.status = "failed";    //

        return block;
    }

    if(m_X_rowvec.size() <= 0)         // Check size of data to not be empty
    {
        cout << "input data is empty" << endl;

        block.extended = 0;         //
        block.shortend_et = -1;     //
        block.k_et = -1;            //
        block.n_et = -1;            //  Just set non-sense numbers
        block.primPoly_et = -1;     //
        block.t = -1;               //
        block.status = "failed";    //

        return block;
    }

    if(m_threshold > 1 || m_threshold < 0)      // Check threshold so not be invalid
    {
        cout << "threshold is invalid" << endl;

        block.extended = 0;         //
        block.shortend_et = -1;     //
        block.k_et = -1;            //
        block.n_et = -1;            //  Just set non-sense numbers
        block.primPoly_et = -1;     //
        block.t = -1;               //
        block.status = "failed";    //

        return block;
    }

    if(m_X_rowvec.size() % m_n != 0)       // Check to see if size of data is dividable by n or not
        m_X_rowvec = m_X_rowvec.cols(0, (m_X_rowvec.size() - 1) - (m_X_rowvec.size() % m_n));

    //___________________________ now we start main algorithm ___________________________

    m_X_rowvec = m_X_rowvec.cols(0, (floor((double)(m_X_rowvec.size() / m_n)) * m_n) - 1);

    if(log2(m_n) - (double)floor(log2(m_n)) < 10e-5)        // Check to see is floor(log2(m_n)) == log2(m_n)
    {
        block.extended = 1;
        m_n = m_n - 1;
        m_X_mat = reshape(m_X_rowvec, m_n + 1, m_X_rowvec.n_elem / ( m_n + 1)).t();
        m_X_mat = m_X_mat.cols(0, m_n - 1);
        m_m = ceil(log2(m_n + 1));
        m_n0 = (1 << m_m) - 1;      // 2 ^ m_m - 1
    }
    else
    {
        m_m = ceil(log2(m_n + 1));
        m_n0 = (1 << m_m) - 1;      // 2 ^ m_m - 1
        shorten = (1 << m_m) - 1 - m_n;     // 2 ^ m_m - 1 - n

        m_X_mat = reshape(m_X_rowvec, m_n, m_X_rowvec.n_elem / m_n).t();
        m_X_mat = join_rows(zeros<uchar_mat>(m_X_mat.n_rows, shorten), m_X_mat);
    }

    vec prim_poly_list = allPrimPoly(m_m);

    if(m_X_mat.n_rows > 1000)
        m_X_mat = m_X_mat.rows(0, 999);         // Just take the first 1000 rows for calculations

    int M = m_X_mat.n_rows;         // M is the number of messages (frame length)

    // Main Algorithm
    for (int i = 0; i < prim_poly_list.n_elem; i++)
    {
        primPoly = prim_poly_list(i);
        block.status = "success";

        uchar_mat x(1, 1);     x(0, 0) = 2;
        object alpha = myGaliosField(x, m_m, primPoly);

        x = zeros<uchar_mat>(m_n0, m_n0);
        object F = myGaliosField(x, m_m, primPoly);

        umat a;
        rowvec rvec;
        rvec.resize(m_n0);
        for (int j = 0; j < m_n0; j++)
            rvec(j) = j;

        for (int j = 0; j < m_n0; j++)
        {
            umat col = myPower(myMPower(alpha, j), rvec);
            a = join_rows(a, col.row(0).t());
        }

        F.x = a;

        uchar_mat temp_mat = fliplr(m_X_mat);
        umat r = myMTimes(temp_mat, F);

        mat left = ones<mat>(M, 1);
        mat right = ones<mat>(M, 1);
        mat A = join_rows(join_rows(left, conv_to<mat>::from(r)) , right);
        umat B = (A == 0);
        mat Y = conv_to<mat>::from(B.cols(1, B.n_cols - 1)) - conv_to<mat>::from(B.cols(0, B.n_cols - 2));

        umat sz_nz;

        for (int j = 0; j < Y.n_rows; j++)
        {
            rowvec row = Y.row(j);
            uvec idx_neg = find(row == -1);
            uvec idx_pos = find(row == +1);
            long min_len = min(idx_neg.n_elem, idx_pos.n_elem);
            ivec zero_lengths = conv_to<ivec>::from(idx_neg.head(min_len)) - conv_to<ivec>::from(idx_pos.head(min_len));

            if(zero_lengths.is_empty() || zero_lengths.max() == 1)
                continue;

            uword lz_idx = -1;
            int max_lz = zero_lengths.max(lz_idx);

            uvec tmp = find(row == 1);

            if(max_lz % 2 == 1)
            {
                if(tmp(lz_idx)+1 + 1 <= m_n0)
                {
                    uvec new_row = {tmp(lz_idx)+1 + 1, static_cast<unsigned long long>(max_lz - 1)};
                    sz_nz = join_vert(sz_nz, new_row.t());
                }

                if(tmp(lz_idx)+1 - 1 >= 1)
                {
                    uvec new_row = {tmp(lz_idx)+1 - 1, static_cast<unsigned long long>(max_lz - 1)};
                    sz_nz = join_vert(sz_nz, new_row.t());
                }
            }

            else
            {
                uvec new_row = {tmp(lz_idx)+1, static_cast<unsigned long long>(max_lz)};
                sz_nz = join_vert(sz_nz, new_row.t());
            }
        }

        if(!sz_nz.is_empty())
        {
            vec col1 = conv_to<vec>::from(sz_nz.col(0));
            uvec idx = find_finite(col1);
            vec sz_tmp = col1.elem(idx);

            vec sz_set = unique(sz_tmp);

            if(sz_set.is_empty())
                continue;

            vec counts = zeros<vec>(sz_set.n_rows, sz_set.n_cols);

            for (int ii = 0; ii < sz_set.n_elem; ii++)
                counts(ii) = accu(sz_tmp == sz_set(ii));

            uword tmp_idx;
            double max_tmp = counts.max(tmp_idx);

            if((max_tmp < m_threshold * Y.n_rows))
            {
                block.status = "failed";
                continue;
            }

            double svcr = sz_set(tmp_idx);

            if(svcr != 2)
            {
                block.status = "failed";
                continue;
            }

            uvec t_idx = find(sz_nz.col(0) == svcr);
            umat tt = sz_nz.rows(idx);
            vec lz = conv_to<vec>::from(tt.col(1));

            vec lz_set = unique(lz);

            if(lz_set.is_empty())
            {
                block.status = "failed";
                continue;
            }

            counts = zeros<vec>(lz_set.n_elem);
            for (int ii = 0; ii < lz_set.n_elem; ii++)
                counts(ii) = accu(lz == lz_set(ii));

            max_tmp = counts.max(tmp_idx);

            if((max_tmp < m_threshold * Y.n_rows))
            {
                block.status = "failed";
                continue;
            }

            mlcr = lz_set(tmp_idx);

            if(block.status == "success")
            {
                mat T = myBCHNumerr(m_n0);
                uvec idx = find(T.col(2) == mlcr / 2);

                if(idx.n_elem > 0)
                {
                    mat temp = T.rows(idx);
                    k = temp.col(1);
                }

                else
                    k.reset();

                break;
            }
        }
    }

    block.n_et = m_n0;
    block.shortend_et = shorten;

    if(!k.is_empty())
    {
        block.k_et = k(0);
        block.t = mlcr/2;
        block.primPoly_et = primPoly;
    }

    return block;
}

vec BCH_Blind::allPrimPoly(int m)
{
    switch (m)
    {
    case 2: return vec{7};
    case 3: return vec{11, 13};
    case 4: return vec{19, 25};
    case 5: return vec{37, 41, 47, 55, 59, 61};
    case 6: return vec{67, 91, 97, 103, 109, 115};
    case 7: return vec{131, 137, 143, 145, 157, 167, 171, 185, 191, 193, 203, 211, 213, 229, 239, 241, 247, 253};
    case 8: return vec{285, 299, 301, 333, 351, 355, 357, 361, 369, 391, 397, 425, 451, 463, 487, 501};

    default: cout << "m is not a supported value" << endl; return vec{};
    }
}

object BCH_Blind::myGaliosField(uchar_mat& x, int64_t m, int64_t prim_poly)
{
    object obj;
    obj.p_vec = {3, 7, 11, 19, 37, 67, 137, 285, 529, 1033, 2053, 4179, 8219, 17475, 32771, 69643};
    obj.m = m;
    obj.prim_poly = prim_poly;
    obj.is_prime = 1;
    obj.q = (1 << obj.m) - 1;       // 2^obj.m-1;

    if(x.min() >= 0 && x.max() <= obj.q)
        obj.x = conv_to<umat>::from(x);

    rowvec tmp(2 * obj.m - 1);
    for (int i = 0; i < tmp.n_elem ; i++)
        tmp(i) = (1 << i);

    obj.twos = tmp;

    return obj;
}

object BCH_Blind::myMPower(object A, double pwr)
{
    pwr = round(pwr);
    umat z;

    if(A.m != GF_TABLE_M || A.prim_poly != GF_TABLE_PRIM_POLY)
        getTables(A);

    z = myPower(A, rowvec{pwr});

    object out = A;
    out.x = z;

    return out;
}

void BCH_Blind::getTables(object A)
{
    uvec index;

    if(GF_MAIN_PRIM_POLY.prim_poly[A.m].is_empty())
        index = {};
    else
        index = find(GF_MAIN_PRIM_POLY.prim_poly[A.m - 1] == A.prim_poly);

    if(!index.is_empty())       // Checking index so must not be empty
    {
        GF_TABLE1 = conv_to<umat>::from(GF_MAIN_TABLE1.table1[A.m - 1].t()).cols(index);
        GF_TABLE2 = conv_to<umat>::from(GF_MAIN_TABLE2.table2[A.m - 1].t()).cols(index);
    }

    GF_TABLE_M = A.m;
    GF_TABLE_PRIM_POLY = A.prim_poly;
}

umat BCH_Blind::myPower(object x, rowvec y)
{
    getTables(x);

    x.x = umat(y.n_rows, y.n_cols, arma::fill::value(x.x(0, 0)));

    if(x.m != GF_TABLE_M || x.prim_poly != GF_TABLE_PRIM_POLY)
        getTables(x);

    umat E2P = GF_TABLE1;
    umat P2E = GF_TABLE2;

    object w = x;
    umat v = w.x;
    y = round(y);

    umat z(v.n_rows, v.n_cols, arma::fill::zeros);

    if(!E2P.is_empty() && !P2E.is_empty())
    {
        umat r = z;
        umat nzs = (v != 0);

        uvec idx = find(nzs);       // Get indices where nzs is true

        uvec y_sel = conv_to<uvec>::from(y.elem(idx));      // Select y and convert to uint32

        uvec v_sel = P2E(conv_to<uvec>::from(v.elem(idx) - 1));      // Apply P2E to selected v

        uvec temp = (y_sel % v_sel);
        temp.transform([x](uint32_t val){return val % x.q;});       // All of these line implement this line :
        r.elem(find(nzs)) = conv_to<umat>::from(temp);

        z(find((r == 0) % nzs)).fill(1);       // This line is equivalent to this matlab line : z(r==0 & nzs) = uint32(1);
        z.elem(find(r != 0)) = E2P(r(find(r != 0)) - 1);        // This line is equivalent to this matlab line : z(r~=0) = E2P(r(r~=0));
        z.elem(find(y==0)).fill(1);       // This line is equivalent to this matlab line : z(y==0) = uint32(1);
    }

    else
    {
        for (int idx = 0; idx < v.n_elem; idx++)
        {
            if(y(idx) * v(idx) == 0)
                z(idx) = y(idx) == 0;

            else
            {
                z(idx) = v(idx);
                for (int count = 1; count < y(idx); count++)
                    z(idx) = myGFMultiply(z(idx), v(idx), x.prim_poly, x.twos, x.m, x.q, E2P, P2E);
            }
        }
    }

    w.x = z;
    umat out = w.x;

    return out;
}

int BCH_Blind::myGFMultiply(uint32_t a, uint32_t b, int64_t irr, rowvec twos, int64_t m, int64_t q, umat E2P, umat P2E)
{
    int64_t y = -1;
    if(irr <= 3)
        return a * b;

    if(a * b == 0)
        return 0;

    if(!E2P.is_empty() && !P2E.is_empty())
    {
        int r = (P2E(a) + P2E(b)) % q;
        y = r == 0 ? 1 : E2P(r);
    }
    else
    {
        y = 0;
        for (int i = 0; i < m; i++)
            y = y ^ static_cast<int>(((a >> i) & 1) * (b) * twos(i));

        int64_t degY = floor(log2(y));
        int64_t degQ = degY - m;

        if(degQ < 0)
            return y;

        else if (degQ == 0)
        {
            y ^= irr;
            return y;
        }

        else
        {
            for (int idx = degQ; idx >= 0; idx--)
            {
                uint32_t bit = ((m + idx + 1) >> 0) & 1;    // Least significant bit
                uint32_t xor_term = bit * irr * twos(idx);  // Multiply by irr and
                y ^= xor_term;
            }
        }
    }

    return y;
}

umat BCH_Blind::myGFMtimes(umat &x, umat &y, int64_t irr, rowvec twos, int64_t m, int64_t q, umat E2P, umat P2E)
{
    umat z(x.n_rows, y.n_cols, fill::zeros);
    umat x_t = x.t();

    for (int idx = 0; idx < x_t.n_cols; idx++)
    {
        for (int jdx = 0; jdx < y.n_cols; jdx++)
        {
            uvec temp_1 = x_t.col(idx);
            uvec temp_2 = y.col(jdx);
            umat prod = myGFtimes(temp_1, temp_2, irr, twos, m, q, E2P, P2E);

            for (int k = 0; k < x_t.n_rows; k++)
                z(idx, jdx) = static_cast<int>(z(idx, jdx)) ^ static_cast<int>(prod(k));
        }
    }

    return z;
}

umat BCH_Blind::myGFtimes(uvec &x, uvec &y, int64_t irr, rowvec twos, int64_t m, int64_t q, umat E2P, umat P2E)
{
    umat z;
    if(!E2P.is_empty() && !P2E.is_empty())
    {
        z = conv_to<umat>::from(zeros(x.n_rows, x.n_cols));
        umat nzs = conv_to<umat>::from((x != 0) % (y != 0));
        umat r = z;

        uvec idx = find(nzs);
        uvec x_sel = x(idx);
        uvec y_sel = y(idx);
        uvec vals = P2E.elem(x_sel - 1) + P2E.elem(y_sel - 1);
        vals.transform([q](uint32_t val){return val % q;});
        r.elem(find(nzs)) = conv_to<umat>::from(vals);

        uvec idx_temp = find((r == 0) % nzs);
        z.elem(idx_temp).fill(1);

        z.elem(find(r != 0)) = E2P.elem(r.elem(find(r != 0)) - 1);
    }

    else
    {
        if(irr <= 3)
            z = x % y;
        else
        {
            z = conv_to<umat>::from(zeros(x.n_rows, x.n_cols));

            for (int idx = 0; idx < y.n_elem; idx++)
                z(idx) = myGFMultiply(x(idx), y(idx), irr, twos, m, q, E2P, P2E);
        }
    }

    return z;
}

umat BCH_Blind::myMTimes(uchar_mat &x, object y)
{
    object xx = myGaliosField(x, y.m, y.prim_poly);

    if(xx.m != GF_TABLE_M || xx.prim_poly != GF_TABLE_PRIM_POLY)
        getTables(xx);

    if(xx.x.n_cols != y.x.n_rows || xx.x.n_cols <= 1 || xx.x.n_rows <= 1 || y.x.n_rows <= 1 || y.x.n_cols <= 1)
    {
        cout << "comm:gf_mtimes:MismatchSize" << endl;
        return umat(1, 1);
    }

    umat final_out = myGFMtimes(xx.x, y.x, xx.prim_poly, xx.twos, xx.m, xx.q, GF_TABLE1, GF_TABLE2);

    return final_out;
}

mat BCH_Blind::myBCHNumerr(int n)
{
    switch (n)
    {
        case 7:   return mat{{7, 4, 1}};
        case 15:  return mat{{15, 11, 1}, {15, 7, 2}, {15, 5, 3}};
        case 31:  return mat{{31, 26, 1}, {31, 21, 2}, {31, 16, 3}, {31, 11, 5}, {31, 6, 7}};
        case 63:  return mat{{63, 57, 1}, {63, 51, 2}, {63, 45, 3}, {63, 39, 4}, {63, 36, 5}, {63, 30, 6}, {63, 24, 7}, {63, 18, 10}, {63, 16, 11}, {63, 10, 13}, {63, 7, 15}};
        case 127: return mat{{127, 120, 1}, {127, 113, 2}, {127, 106, 3}, {127, 99, 4}, {127, 92, 5}, {127, 85, 6}, {127, 78, 7}, {127, 71, 9}, {127, 64, 10}, {127, 57, 11}, {127, 50, 13}, {127, 43, 14}, {127, 36, 15}, {127, 29, 21}, {127, 22, 2}, {127, 15, 27}, {127, 8, 31}};
        case 255: return mat{{255, 247, 1}, {255, 239, 2}, {255, 231, 3}, {255, 223, 4}, {255, 215, 5}, {255, 207, 6}, {255, 199, 7}, {255, 191, 8}, {255, 187, 9}, {255, 179, 10}, {255, 171, 11}, {255, 163, 12}, {255, 155, 13}, {255, 147, 14}, {255, 139, 15}, {255, 131, 18}, {255, 123, 19}, {255, 115, 21}, {255, 107, 22}, {255, 99, 23}, {255, 91, 25}, {255, 87, 26}, {255, 79, 27}, {255, 71, 29}, {255, 63, 30}, {255, 55, 31}, {255, 47, 42}, {255, 45, 43}, {255, 37, 45}, {255, 29, 47}, {255, 21, 55}, {255, 13, 59}, {255, 9, 63}};
    }

    return mat{};
}
