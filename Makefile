.PHONY: all clean plot
all: iir iir_f
clean:
	rm iir iir_f

CASE ?= 1
plot: all
	./iir $(CASE) | ./plot.py

iir: iir.c iir_f.c
	cc iir.c -o iir -lm
iir_f: iir_f.c
	cc iir_f.c -o iir_f -lm

