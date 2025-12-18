clear all, close all;
pkg load signal;

fs = 1000;
fg = 100;
[b, a] = butter(6, fg/(fs/2));
sos = tf2sos(b, a);

round(sos * 2^14)

