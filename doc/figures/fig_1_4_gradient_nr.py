# All rights reserved (c) 2022, Vladimir Botka <vbotka@gmail.com>
# Simplified BSD License, https://opensource.org/licenses/BSD-2-Clause

from hamming_digital_filters import gradient as hdfg

import numpy as np
import matplotlib.pyplot as plt
import math
import random

dx = 0.01
x = np.array([i * dx for i in range(2 * 314)])
y = np.array([math.sin(xi) for xi in x])
z = np.array([math.cos(xi) for xi in x])

dydx_5 = hdfg.gradient_nr(y, dx, 0.1, 5)
dydx_7 = hdfg.gradient_nr(y, dx, 0.1, 7)
dydx_9 = hdfg.gradient_nr(y, dx, 0.1, 9)

r = np.array([yi + (random.random() - 0.5) / 50 for yi in y])
drdx_7 = hdfg.gradient_nr(r, dx, 0.1, 7)
drdx_n = np.gradient(r, dx)

plt.title("Hamming nonrecursive differentiation filter fc=0.1 N=5,7,9")
plt.xlabel("x")
plt.ylabel("sin(x), cos(x),  d(sin(x))/dx")
plt.plot(x, y, x[5:-5], dydx_5[5:-5], x[7:-7], dydx_7[7:-7], x[9:-9], dydx_9[9:-9])
plt.show()

plt.title("Hamming nonrecursive differentiation filter fc=0.1 N=7")
plt.xlabel("x")
plt.ylabel("sin(x), cos(x),  d(sin(x))/dx")
plt.plot(x, y, x, z, x[7:-7], dydx_7[7:-7])
plt.show()

plt.title("Hamming nonrecursive differentiation filter fc=0.1 N=7")
plt.xlabel("x")
plt.ylabel("sin(x)+random()/50, cos(x),  d(sin(x))/dx")
plt.plot(x, r, x, z, x[7:-7], drdx_7[7:-7])
plt.show()

plt.title("NumPy gradient filter")
plt.xlabel("x")
plt.ylabel("sin(x)+random()/50, cos(x),  d(sin(x))/dx")
plt.plot(x, r, x, z, x, drdx_n)
plt.show()
