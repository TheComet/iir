clear all, close all;
pkg load signal;

fs = 1000;
fg = 100;
[b, a] = butter(6, fg/(fs/2));
sos = tf2sos(b, a)

round(sos * 2^14)

t = linspace(0, 0.3, fs*0.3);
f = 5;
x = cos(2*pi*t*f);
y = sosfilt(sos, x);
plot(t, x, '-o', t, y, '-o');
