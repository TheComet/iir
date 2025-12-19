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
    int16_t b0, b1, b2;
    int16_t a1, a2;
};

const int Q = 14;
const int K = 1 << (Q-1);
static const struct biquad_coeffs coeffs[] = {
   {     6,    11,     6, -16909,  4517 },
   { 16384, 32767, 16384, -18727,  6763 },
   { 16384, 32641, 16257, -23009, 12057 },
};

struct biquad_state {
    int16_t x1, x2;
    int16_t y1, y2;
    int16_t saved_fraction;
};

static int16_t filter(int16_t input, struct biquad_state* state, const struct biquad_coeffs* coeffs) {
    int i;
    int16_t output, frac;

    for (i = 0; i != 3; ++i) {
        int32_t acc = state[i].saved_fraction;

        acc += coeffs[i].b0 * input;
        acc += coeffs[i].b1 * state[i].x1;
        acc += coeffs[i].b2 * state[i].x2;
        acc -= coeffs[i].a1 * state[i].y1;
        acc -= coeffs[i].a2 * state[i].y2;

        output = (int16_t)(acc >= 0 ? (acc+K)>>Q : (acc-K)>>Q);
        frac = (int16_t)(acc - (output << Q));

        state[i].x2 = state[i].x1;
        state[i].x1 = input;
        state[i].y2 = state[i].y1;
        state[i].y1 = output;
        state[i].saved_fraction = frac;

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
        int16_t x = (int16_t)(cos(2*M_PI*f*t) * (1<<Q));
        int16_t y = filter(x, state, coeffs);
        printf("%d,%d\n", x, y);
        t += 1/fs;
    }
}

static void impulse_response_test(void) {
    struct biquad_state state[3] = {0};
    double t = 0.0;
    const double fs = 1000;
    while (t < 0.3) {
        int16_t x = t == 0.0 ? (int16_t)(1<<Q) : 0;
        int16_t y = filter(x, state, coeffs);
        printf("%d,%d\n", x, y);
        t += 1/fs;
    }
}

static void passband_test(void) {
    struct biquad_state state[3] = {0};
    double t = 0.0;
    const double fs = 1000;
    const double f = 5;
    while (t < 0.3) {
        int16_t x = (int16_t)(cos(2*M_PI*f*t) * (1<<Q));
        int16_t y = filter(x, state, coeffs);
        printf("%d,%d\n", x, y);
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
