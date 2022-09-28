# All rights reserved (c) 2022, Vladimir Botka <vbotka@gmail.com>
# Simplified BSD License, https://opensource.org/licenses/BSD-2-Clause

from hamming_digital_filters import gradient as hdfg
from hamming_digital_filters import lowpass as hdfl

import numpy as np
import matplotlib.pyplot as plt
import math
import random


dx = 0.01
x = np.array([i * dx for i in range(2 * 314)])
y = np.array([math.sin(xi) for xi in x])
z = np.array([math.cos(xi) for xi in x])

r = np.array([yi + (random.random() - 0.5) / 100 for yi in y])
drdx_7 = hdfg.gradient_nr(r, dx, 0.1, 7)

s1 = hdfl.lowpass_mt(r, '3_1')
ds1dx_7 = hdfg.gradient_nr(s1, dx, 0.1, 7)
s2 = hdfl.lowpass_mt(s1, '3_1')
ds2dx_7 = hdfg.gradient_nr(s2, dx, 0.1, 7)

plt.title("Hamming nonrecursive differentiation filter fc=0.1 N=7")
plt.xlabel("x")
plt.ylabel("sin(x)+random()/50, cos(x),  d(sin(x))/dx")
plt.plot(x, r, x, z, x[7:-7], drdx_7[7:-7])
plt.show()

plt.title("Hamming Low-Pass monotone filter p=3 q=1")
plt.xlabel("x")
plt.ylabel("sin(x)+random()/50")
plt.plot(x, r, x, s1)
plt.show()

plt.title("Hamming nonrecursive differentiation filter fc=0.1 N=7")
plt.xlabel("x")
plt.ylabel("sin(x)+random()/50, cos(x),  d(sin(x))/dx")
plt.plot(x, s1, x, z, x[7:-7], ds1dx_7[7:-7])
plt.show()

plt.title("Hamming nonrecursive differentiation filter fc=0.1 N=7")
plt.xlabel("x")
plt.ylabel("sin(x)+random()/50, cos(x),  d(sin(x))/dx")
plt.plot(x, s2, x, z, x[7:-7], ds2dx_7[7:-7])
plt.show()
