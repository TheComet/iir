Demo on how an IIR  filter can be implemented using cascaded biquads and shared
state optimization.

Compile the executable with:

```
make
```

The executable has six tests built into it. You can specify which test you want
to run as the first argument:

```
./iir <test number>
```

The  program  will  write  the  input  and output signals to stdout. There is a
```plot.py``` script which  can read the output and generate a visual. You will
need Python 3 with matplotlib and numpy installed:

```
./iir 1 | ./plot.py
```

```make plot``` is a shorthand for ```make &&  ./iir  1 | ./plot.py```. You can
set the test number using the environment variable ```CASE```:

```
make plot CASE=1
```

Additionally, there is an identical implementation ```iir_f.c``` using floating
point numbers instead of integers, for reference.
