#include <armadillo>
#include <iostream>
#include <BCH_Blind.h>

using namespace std;
using namespace arma;

int main(int argc, char *argv[])
{
    // Initial parameters
    int extend = 0;
    int shortend = 0;
    int N = 8;
    double threshold = 0.2;

    mat tmp;
    tmp.load("./r.bin", arma::raw_binary);      // This is the encoded data

    uchar_rowvec r = conv_to<uchar_rowvec>::from(tmp.t());

    BCH_Blind blind_test(r, N - shortend + extend, threshold);
    dataParameters out = blind_test.run_BCH_Blind();

    cout << "status : " << out.status << endl;
    cout << "extended : " << out.extended << endl;
    cout << "k_et : " << out.k_et << endl;
    cout << "n_et : " << out.n_et << endl;
    cout << "shortend_et : " << out.shortend_et << endl;
    cout << "primPoly_et : " << out.primPoly_et << endl;
    cout << "t : " << out.t << endl;
}
