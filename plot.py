#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import sys

pairs = (
    tuple(map(float, pair.split(",")))
    for pair in sys.stdin.read().split("\n")
    if len(pair) > 0
)
x, y = map(np.array, zip(*pairs))
t = np.linspace(0, 0.3, len(x))

plt.plot(t, x)
plt.scatter(t, x)
plt.plot(t, y)
plt.scatter(t, y)
plt.show()
