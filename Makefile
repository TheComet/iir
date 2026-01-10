CC = cc
LD = cc

CFLAGS = -g -W -Wall -Wextra -Wpedantic -Werror
LDFLAGS = -lm

all: iir iir_f
clean:
	rm iir iir_f

CASE ?= 1
plot: all
	./iir $(CASE) | ./plot.py

iir: iir.c iir_f.c Makefile
	$(CC) $(CFLAGS) iir.c -o iir $(LDFLAGS)
iir_f: iir_f.c
	$(CC) $(CFLAGS) iir_f.c -o iir_f $(LDFLAGS)

.PHONY: all clean plot
