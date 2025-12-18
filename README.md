Quick demo on how a simple IIR filter is implemented using cascaded biquads

The executable has two tests built into it.  Test  1  will generate a cosine at
f=200  Hz  and pass it through the filter. The filter is a  Butterworth  filter
with fg=100 Hz, meaning, the output of the filter will be an attenuated cosine.

```
make && ./iir 1 | ./plot.py
```

Test 2 measures the impulse response of the filter:

```
make && ./iir 2 | ./plot.py
```

You will need Python 3 with matplotlib and numpy installed.
