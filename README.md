Quick demo on how a simple IIR filter is implemented using cascaded biquads

The executable has three tests built into it.

Test  1 will generate a cosine at f=200 Hz and pass it through the filter.  The
filter is a Butterworth filter  with  fg=100  Hz,  meaning,  the  output of the
filter will be an attenuated cosine.

```
make && ./iir 1 | ./plot.py
```

Test 2 measures the impulse response of the filter:

```
make && ./iir 2 | ./plot.py
```

Test 3 will generate a  cosine  at  f=5  Hz,  which  is  in the passband of the
filter.  The   output   should   be  an  identical  cosine,  with  some  delay.

```
make && ./iir 3 | ./plot.py
```

Additionally, there is an identical implementation ```iir_f.c``` using floating
point numbers instead of integers, for reference.

You will need Python 3 with matplotlib and numpy installed.
