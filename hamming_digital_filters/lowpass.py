# All rights reserved (c) 2022, Vladimir Botka <vbotka@gmail.com>
# Simplified BSD License, https://opensource.org/licenses/BSD-2-Clause

# Digital Filters by R.W.Hamming,third edition
# Dover publications, INC., Mineola, New York
# ISBN 0-486-65088-X

import math


def lowpass_nr(f, fc=0.2, nfp=5):
    ''' Chapter 6.2 A Low-Pass Filter Design. Nonrecursive. '''

    coef = [0.0] * (2 * nfp + 1)
    coef[nfp] = 2.0 * fc
    sf = [0.0] * len(f)

    for k in range(1, nfp):
        xd = math.pi * k
        yd = xd / nfp
        c1 = math.sin(yd) / yd
        c2 = math.sin(2.0 * fc * xd) / xd
        coef[nfp + k] = c1 * c2
        coef[nfp - k] = c1 * c2

    for i in range(0, len(f) - 2 * nfp):
        s = 0.0
        for j in range(0, 2 * nfp + 1):
            s = s + coef[j] * f[i + j]
        sf[nfp + i] = s

    return sf


def lowpass_mt(f, pq='3_1'):
    ''' Chapter 7.5 The Design Of A Smooth Filter. Monotone. '''

    frac = (1 / 16) ** 2
    coef_3_1 = [frac * c for c in [-1, -5, -5, 20, 70, 98, 70, 20, -5, -5, -1]]
    coef = {'3_1': coef_3_1}
    nfp = (len(coef[pq]) - 1) // 2
    sf = [0.0] * len(f)

    for i in range(0, len(f) - 2 * nfp):
        s = 0.0
        for j in range(0, 2 * nfp + 1):
            s = s + coef[pq][j] * f[i + j]
        sf[nfp + i] = s

    return sf
