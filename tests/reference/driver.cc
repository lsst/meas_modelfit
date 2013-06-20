// This is driver code that calls the FORTRAN routine BVND, defined in bvn.f, and writes
// output for a number of input values to be used to test the C++ implementation in
// integrals.h.  Because this relies on some platform-dependent assumptions about how
// to link C++ and FORTRAN, we do not expect to build this code often; instead we'll just
// commit the output file in git.

#include <iostream>

extern "C" {
    double phid_(double * z);
    double bvnd_(double * dh, double * dk, double *r);
}

int main() {
    int const n = 10;
    double const hkMin = -3.0;
    double const hkMax = 3.0;
    double const rMin = -0.98;
    double const rMax = 0.98;
    std::cout.precision(15);
    for (int i = 0; i < n; ++i) {
        double h = hkMin + i * (hkMax - hkMin) / (n-1);
        for (int j = 0; j < n; ++j) {
            double k = hkMin + j * (hkMax - hkMin) / (n-1);
            for (int p = 0; p < n; ++p) {
                double r = rMin + p * (rMax - rMin ) / (n-1);
                std::cout << h << ", " << k << ", " << r << ", " << bvnd_(&h, &k, &r) << std::endl;
            }
        }
    }
    return 0;
}
