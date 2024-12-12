#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cassert>
#include <chrono> // For high-resolution timing
#include <omp.h>  // For parallelization

typedef unsigned long ui32;
typedef unsigned long long ui64;

// High-Resolution Timer Using std::chrono
double dml_micros() {
    auto now = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
                        now.time_since_epoch());
    return duration.count();
}

// Even function for sign flipping
inline int even(int n) {
    return (n % 2 == 0) ? 1 : -1;
}

// Efficient theta function using polynomial expansion
double theta(double t) {
    const double pi = 3.1415926535897932385;
    return (t / 2.0 * log(t / (2.0 * pi)) - t / 2.0 - pi / 8.0 +
            1.0 / (48.0 * t) + 7.0 / (5760.0 * pow(t, 3.0)) +
            31.0 / (80640.0 * pow(t, 5.0)) + 127.0 / (1216512.0 * pow(t, 7.0)));
}

// Horner's method for polynomial evaluation
double horner_eval(double x, const double coeffs[], int degree) {
    double result = coeffs[degree];
    for (int i = degree - 1; i >= 0; --i) {
        result = result * x + coeffs[i];
    }
    return result;
}

// Precompute constants for efficiency using Horner's method
double C(int n, double z) {
    static const double coeffs[5][24] = {
        {0.38268343236508977173, 0.43724046807752044936, 0.13237657548034352332, -0.01360502604767418865, /* ... coefficients for n = 0 */},
        {-0.02682510262837534703, 0.01378477342635185305, 0.03849125048223508223, 0.00987106629906207647, /* ... coefficients for n = 1 */},
        {0.00518854283029316849, 0.00030946583880634746, -0.01133594107822937338, 0.00223304574195814477, /* ... coefficients for n = 2 */},
        {-0.00133971609071945690, 0.00374421513637939370, -0.00133031789193214681, -0.00226546607654717871, /* ... coefficients for n = 3 */},
        {0.00046483389361763382, -0.00100566073653404708, 0.00024044856573725793, 0.00102830861497023219, /* ... coefficients for n = 4 */}
    };

    if (n >= 0 && n < 5) {
        int degree = sizeof(coeffs[n]) / sizeof(coeffs[n][0]) - 1;
        return horner_eval(z, coeffs[n], degree);
    } else {
        return 0.0; // Default for out-of-range n
    }
}

// Optimized Z function with precomputation and parallelization
double Z(double t, int n) {
    const double pi = 3.1415926535897932385;
    int N = sqrt(t / (2.0 * pi));
    double p = sqrt(t / (2.0 * pi)) - N;
    double tt = theta(t);
    double ZZ = 0.0;

    // Parallelize loop for efficiency
    #pragma omp parallel for reduction(+:ZZ)
    for (int j = 1; j <= N; ++j) {
        ZZ += 1.0 / sqrt((double)j) * cos(fmod(tt - t * log((double)j), 2.0 * pi));
    }
    ZZ *= 2.0;

    double R = 0.0;
    for (int k = 0; k <= n; ++k) {
        R += C(k, 2.0 * p - 1.0) * pow(2.0 * pi / t, ((double)k) * 0.5);
    }
    R *= even(N - 1) * pow(2.0 * pi / t, 0.25);

    return ZZ + R;
}

// Count zeros efficiently
void count_zeros(double LOWER, double UPPER, double SAMP) {
    const double pi = 3.1415926535897932385;
    double STEP = 1.0 / SAMP;
    ui64 NUMSAMPLES = floor((UPPER - LOWER) * SAMP + 1.0);
    double prev = 0.0, count = 0.0;

    double t1 = dml_micros();
    #pragma omp parallel for reduction(+:count)
    for (ui64 i = 0; i < NUMSAMPLES; ++i) {
        double t = LOWER + i * STEP;
        double zout = Z(t, 4);
        if (i > 0) {
            if (((zout < 0.0) && (prev > 0.0)) || ((zout > 0.0) && (prev < 0.0))) {
                count++;
            }
        }
        prev = zout;
    }
    double t2 = dml_micros();
    printf("I found %.0lf Zeros in %.3lf seconds\n", count, (t2 - t1) / 1000000.0);
}

int main(int argc, char **argv) {
    if (argc < 4) {
        std::cerr << argv[0] << " START END SAMPLING" << std::endl;
        return -1;
    }

    double LOWER = atof(argv[1]);
    double UPPER = atof(argv[2]);
    double SAMP = atof(argv[3]);

    double estimate_zeros = theta(UPPER) / 3.1415926535897932385;
    printf("I estimate I will find %.3lf zeros\n", estimate_zeros);

    count_zeros(LOWER, UPPER, SAMP);

    return 0;
}
