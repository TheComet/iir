#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#define N 4
const int Q = 14;
const int K = 1 << (Q-1);

const int input_Q = 12;
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
#if N == 16
    {    1149,    4242,    3972,  -12215,    2298 },
    {    1149,    4036,    3744,  -12353,    2514 },
    {    1149,    3712,    3391,  -12790,    3009 },
    {    1149,    3292,    2916,  -13575,    3889 },
    {    1149,    2844,    2394,  -14301,    4975 },
    {    1149,    2431,    1901,  -15036,    4288 },
    {    1149,    2081,    1477,  -15460,    6419 },
    {    1149,    1802,    1132,  -17137,    8294 },
    {    1149,    1590,     865,  -19371,   10596 },
    {    1149,    1439,     667,  -22110,   13181 },
    {    1149,    1333,     525,  -23654,   10607 },
    {    1149,    1260,     425,  -25063,   15636 },
    {    1149,    1216,     360,  -27377,   16707 },
    {    1149,    1196,     324,  -29339,   16182 },
    {    1149,    1196,     312,  -31206,   15709 },
    {    1149,    3095,    1757,  -23522,    7348 },
#elif N == 12
    {    1149,    3576,    2792,  -15642,    3788 },
    {    1149,    3413,    2614,  -15744,    4043 },
    {    1149,    3132,    2304,  -16255,    4666 },
    {    1149,    2805,    1941,  -17150,    5657 },
    {    1149,    2501,    1603,  -18421,    7009 },
    {    1149,    2234,    1313,  -20196,    8702 },
    {    1149,    1993,    1055,  -21178,    7633 },
    {    1149,    1794,     840,  -21285,    9660 },
    {    1149,    1646,     677,  -22233,   11132 },
    {    1149,    1544,     564,  -23785,   13019 },
    {    1149,    1481,     494,  -25529,   15171 },
    {    1149,    1451,     460,  -18426,    5074 },
#elif N == 6
    {    1148,    2497,    1359,  -16749,    4319 },
    {    1148,    2433,    1294,  -17180,    4852 },
    {    1148,    2331,    1192,  -18079,    5963 },
    {    1148,    2229,    1089,  -19524,    7749 },
    {    1148,    2157,    1016,  -21642,   10367 },
    {    1148,    2126,     985,  -24621,   14049 },
#elif N == 4
    {    1146,    2328,    1182,  -16816,    4401 },
    {    1146,    2292,    1146,  -17807,    5627 },
    {    1146,    2257,    1111,  -19984,    8318 },
    {    1146,    2293,    1146,  -23783,   13013 },
#elif N == 3
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
    const double f = 150;
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

static void mixed_signal_test(void) {
    struct biquad_state state[N] = {0};
    double t = 0.0;
    const double fs = 1000;
    const double f1 = 80;
    const double f2 = 150;
    while (t < 0.3) {
        double xf = cos(2*M_PI*f1*t) + cos(2*M_PI*f2*t);
        int16_t x = (int16_t)(xf * INPUT_SCALE);
        int16_t y = filter(x, state, coeffs);
        printf("%d,%d\n", x, y);
        t += 1/fs;
    }
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <test index>\n", argv[0]);
        fprintf(stderr, "  Test 1: Cosine f=150Hz\n");
        fprintf(stderr, "  Test 2: Impulse reponse\n");
        fprintf(stderr, "  Test 3: Cosine f=5Hz\n");
        fprintf(stderr, "  Test 4: Sine Sweep f=0..250Hz\n");
        fprintf(stderr, "  Test 5: Two cosines, f=80Hz and f=150Hz\n");
        return EXIT_FAILURE;
    }

    switch (atoi(argv[1])) {
        case 1: stopband_test(); break;
        case 2: impulse_response_test(); break;
        case 3: passband_test(); break;
        case 4: sweepsine_test(); break;
        case 5: mixed_signal_test(); break;
        default: fprintf(stderr, "Unknown test index: %s\n", argv[1]); return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
