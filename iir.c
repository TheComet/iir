#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#define N 3
const int Q = 14;
const int K = 1 << (Q-1);

const int input_Q = 14;
const int INPUT_SCALE = (1<<input_Q)-1;

/*
 *        b0 + b1*z^-1 + b2*z^-2
 * H(z) = ----------------------
 *        1 - a1*z^-1 - a2*z^-2
 */
struct biquad_coeffs {
    int16_t b0, b1, b2;
    int16_t a1, a2;
};
struct biquad_state {
    int16_t history[N*2+2];
    int16_t saved_fractions[N];
};

static const struct biquad_coeffs coeffs[N] = {
#if N == 3
   //{     6,    12,     6, -16909,  4518 },
   //{ 16385, 32767, 16385, -18726,  6764 },
   //{ 16385, 32641, 16257, -23009, 12058 },
   { 1144, 2297, 1153, -16909,  4517 },
   { 1144, 2288, 1144, -18727,  6763 },
   { 1144, 2279, 1135, -23009, 12057 },
#elif N == 2
   {  1265,  2529,  1265, -17180,  4852 },
   {  1024,  2048,  1024, -21642, 10367 },
#elif N == 1
    { 1105, 2210, 1105, -18727, 6763 },
#endif
};

static int16_t filter(int16_t input, struct biquad_state* state, const struct biquad_coeffs* coeffs) {
    int i;
    int16_t output, frac;

    for (i = 0; i != N; ++i) {
        int32_t acc = state->saved_fractions[i];

        acc += coeffs[i].b0 * input;
        acc += coeffs[i].b1 * state->history[i*2+0];
        acc += coeffs[i].b2 * state->history[i*2+1];
        acc -= coeffs[i].a1 * state->history[i*2+2];
        acc -= coeffs[i].a2 * state->history[i*2+3];

        output = (int16_t)(acc >> Q);
        frac = (int16_t)(acc - (output << Q));

        state->history[i*2+1] = state->history[i*2+0];
        state->history[i*2+0] = input;
        state->saved_fractions[i] = frac;

        input = output;
    }

    state->history[i*2+1] = state->history[i*2+0];
    state->history[i*2+0] = input;

    return output;
}

static void stopband_test(void) {
    struct biquad_state state[N] = {0};
    double t = 0.0;
    const double fs = 1000;
    const double f = 200;
    while (t < 0.3) {
        int16_t x = (int16_t)(cos(2*M_PI*f*t) * INPUT_SCALE);
        int16_t y = filter(x, state, coeffs);
        printf("%d,%d\n", x, y);
        t += 1/fs;
    }
}

static void impulse_response_test(void) {
    struct biquad_state state[N] = {0};
    double t = 0.0;
    const double fs = 1000;
    while (t < 0.3) {
        int16_t x = t == 0.0 ? (int16_t)INPUT_SCALE : 0;
        int16_t y = filter(x, state, coeffs);
        printf("%d,%d\n", x, y);
        t += 1/fs;
    }
}

static void passband_test(void) {
    struct biquad_state state[N] = {0};
    double t = 0.0;
    const double fs = 1000;
    const double f = 5;
    while (t < 0.3) {
        int16_t x = (int16_t)(cos(2*M_PI*f*t) * INPUT_SCALE);
        int16_t y = filter(x, state, coeffs);
        printf("%d,%d\n", x, y);
        t += 1/fs;
    }
}

static void sweepsine_test(void) {
    struct biquad_state state[N] = {0};
    const double fs = 1000;
    double f = 0;
    double t = 0;
    while (f < fs/4) {
        int16_t x = (int16_t)(sin(2*M_PI*f*t) * INPUT_SCALE);
        int16_t y = filter(x, state, coeffs);
        printf("%d,%d\n", x, y);
        t += 1/fs;
        f += 0.2;
    }
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <test index>\n", argv[0]);
        fprintf(stderr, "  Test 1: Cosine f=200Hz\n");
        fprintf(stderr, "  Test 2: Impulse reponse\n");
        fprintf(stderr, "  Test 3: Cosine f=5Hz\n");
        fprintf(stderr, "  Test 4: Sine Sweep f=0..250Hz\n");
        return EXIT_FAILURE;
    }

    switch (atoi(argv[1])) {
        case 1: stopband_test(); break;
        case 2: impulse_response_test(); break;
        case 3: passband_test(); break;
        case 4: sweepsine_test(); break;
        default: fprintf(stderr, "Unknown test index: %s\n", argv[1]); return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
