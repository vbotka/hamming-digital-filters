# All rights reserved (c) 2022, Vladimir Botka <vbotka@gmail.com>
# Simplified BSD License, https://opensource.org/licenses/BSD-2-Clause

from hamming_digital_filters import gradient as hdfg
from hamming_digital_filters import lowpass as hdfl

import numpy as np
import matplotlib.pyplot as plt
import math
import random
import io


def print_to_string(x):
    output = io.StringIO()
    print("%.3f" % x, file=output, end='')
    contents = output.getvalue()
    output.close()
    return contents


def distance(x, y):
    s = sum([(xi - yi) ** 2 for xi, yi in zip(x, y)])
    return round(s / len(x) * 1000, 3)


fc = 0.1
nfp = 7
nnfp = (2 * nfp)
noise = 1 / 100

dx = 0.01
x = np.array([i * dx for i in range(2 * 314)])
y = np.array([math.sin(xi) for xi in x])
z = np.array([math.cos(xi) for xi in x])

r = np.array([yi + (random.random() - 0.5) * noise for yi in y])
dist_r = distance(y, r)

drdx = hdfg.gradient_nr(r, dx, fc, nfp)
dist_d = distance(z[nfp:-nfp], drdx[nfp:-nfp])

s1 = hdfl.lowpass_mt(r, '3_1')
dist_s1 = distance(y, s1)

ds1dx = hdfg.gradient_nr(s1, dx, fc, nfp)
dist_d1 = distance(z[nnfp:-nnfp], ds1dx[nnfp:-nnfp])

s2 = hdfl.lowpass_mt(s1, '3_1')
dist_s2 = distance(y, s2)

ds2dx = hdfg.gradient_nr(s2, dx, fc, nfp)
dist_d2 = distance(z[nnfp:-nnfp], ds2dx[nnfp:-nnfp])

plt.title("gradient_nr fc=0.1, N=7, dist f="
          + print_to_string(dist_r) + ", dist df="
          + print_to_string(dist_d))
plt.xlabel("x")
plt.ylabel("sin(x)+random()/100, cos(x),  d(sin(x))/dx")
plt.plot(x, r, x, z, x[nfp:-nfp], drdx[nfp:-nfp])
plt.show()

plt.title("lowpass_mt p=3, q=1, dist f="
          + print_to_string(dist_s1))
plt.xlabel("x")
plt.ylabel("sin(x)+random()/100")
plt.plot(x, r, x, s1)
plt.show()

plt.title("gradient_nr fc=0.1, N=7, dist f="
          + print_to_string(dist_s1) + ", dist df="
          + print_to_string(dist_d1))
plt.xlabel("x")
plt.ylabel("sin(x)+random()/50, cos(x),  d(sin(x))/dx")
plt.plot(x[nfp:-nfp], s1[nfp:-nfp], x, z, x[nnfp:-nnfp], ds1dx[nnfp:-nnfp])
plt.show()

plt.title("gradient_nr fc=0.1, N=7, dist f="
          + print_to_string(dist_s2) + ", dist df="
          + print_to_string(dist_d2))
plt.xlabel("x")
plt.ylabel("sin(x)+random()/100, cos(x),  d(sin(x))/dx")
plt.plot(x[nfp:-nfp], s2[nfp:-nfp], x, z, x[nnfp:-nnfp], ds2dx[nnfp:-nnfp])
plt.show()
