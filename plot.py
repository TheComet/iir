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

color_table = (
    ("#FF6060", "#FF8080", "#C02020"),
    ("#60FFFF", "#A0FFFF", "#20C0C0"),
    ("#FF8000", "#FFA000", "#C07010"),
    ("#FF60FF", "#FFA0FF", "#C020C0"),
    ("#6060FF", "#A0A0FF", "#2020C0"),
    ("#60FF60", "#A0FFA0", "#20C020"),
)

plt.plot(t, x, color=color_table[0][0])
plt.scatter(t, x, color=color_table[0][1])
plt.plot(t, y, color=color_table[1][0])
plt.scatter(t, y, color=color_table[1][1])
plt.gca().set_facecolor("#011627")
plt.grid(visible=True)
plt.show()
