clear all, close all;
pkg load signal;

fs = 1000; fg = 100;
N=4 % Number of biquads
[B, A] = butter(N*2, fg/(fs/2));
[sos, G] = tf2sos(B, A);

freqz(B, A, 1024, fs);

% Redistribute gain evenly among numerator coefficients
num_stages = length(sos(:,1));
gain_per_stage = G^(1/num_stages);
sos(:,1:3) *= gain_per_stage;

% Determine maximum Q format to fit into a signed 16-bit integer (2^15)
largest_coeff = max(max(abs(sos)));
Q = floor(log2(2^15/largest_coeff))
% Print output as C code
sos_q = round(sos*2^Q);
for n = 1:length(sos_q(:,1))
    printf("    { %7d, %7d, %7d, %7d, %7d },\n",
        sos_q(n,1), sos_q(n,2), sos_q(n,3), sos_q(n,5), sos_q(n,6));
endfor

% Example of filtering a cosine
%t = linspace(0, 0.3, fs*0.3);
%f = 5;
%x = cos(2*pi*t*f);
%y = sosfilt(sos, x);
%plot(t, x, '-o', t, y, '-o');
