#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

/*
 *        b0 + b1*z^-1 + b2*z^-2
 * H(z) = ----------------------
 *        1 - a1*z^-1 - a2*z^-2
 */
struct biquad_coeffs {
    double b0, b1, b2;
    double a1, a2;
};

static const struct biquad_coeffs coeffs[] = {
    { 3.4054e-04, 6.8373e-04, 3.4320e-04, -1.0321e+00, 2.7571e-01 },
    { 1.0000e+00, 2.0000e+00, 1.0000e+00, -1.1430e+00, 4.1280e-01 },
    { 1.0000e+00, 1.9922e+00, 9.9224e-01, -1.4044e+00, 7.3592e-01 },
};

struct biquad_state {
    double x1, x2;
    double y1, y2;
};

static double filter(double input, struct biquad_state* state, const struct biquad_coeffs* coeffs) {
    int i;
    double output;

    for (i = 0; i != 3; ++i) {
        output = 0;
        output += coeffs[i].b0 * input;
        output += coeffs[i].b1 * state[i].x1;
        output += coeffs[i].b2 * state[i].x2;
        output -= coeffs[i].a1 * state[i].y1;
        output -= coeffs[i].a2 * state[i].y2;

        state[i].x2 = state[i].x1;
        state[i].x1 = input;
        state[i].y2 = state[i].y1;
        state[i].y1 = output;

        input = output;
    }

    return output;
}

static void stopband_test(void) {
    struct biquad_state state[3] = {0};
    double t = 0.0;
    const double fs = 1000;
    const double f = 200;
    while (t < 0.3) {
        double x = cos(2*M_PI*f*t);
        double y = filter(x, state, coeffs);
        printf("%f,%f\n", x, y);
        t += 1/fs;
    }
}

static void impulse_response_test(void) {
    struct biquad_state state[3] = {0};
    double t = 0.0;
    const double fs = 1000;
    while (t < 0.3) {
        double x = t == 0.0 ? 1.0 : 0.0;
        double y = filter(x, state, coeffs);
        printf("%f,%f\n", x, y);
        t += 1/fs;
    }
}

static void passband_test(void) {
    struct biquad_state state[3] = {0};
    double t = 0.0;
    const double fs = 1000;
    const double f = 5;
    while (t < 0.3) {
        double x = cos(2*M_PI*f*t);
        double y = filter(x, state, coeffs);
        printf("%f,%f\n", x, y);
        t += 1/fs;
    }
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <test index>\n", argv[0]);
        return EXIT_FAILURE;
    }

    switch (atoi(argv[1])) {
        case 1: stopband_test(); break;
        case 2: impulse_response_test(); break;
        case 3: passband_test(); break;
        default: fprintf(stderr, "Unknown test index: %s\n", argv[1]); return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
