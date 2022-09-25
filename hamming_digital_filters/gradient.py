# All rights reserved (c) 2022, Vladimir Botka <vbotka@gmail.com>
# Simplified BSD License, https://opensource.org/licenses/BSD-2-Clause

# Digital Filters by R.W.Hamming,third edition
# Dover publications, INC., Mineola, New York
# ISBN 0-486-65088-X

import math


def gradient_nr(f, dx, fc=0.2, nfp=5):
    ''' Chapter 6.4 Differentiation filter nonrecursive '''

    PI = 3.141592654
    omega = 2.0 * PI * fc
    f3 = nfp + 1
    coef = [0.0] * (2 * nfp + 1)
    dfx = [0.0] * len(f)
    for i in range(1, nfp):
        f2 = i * omega
        smarg = PI * i / f3
        coef[nfp + i] = (math.sin(f2) / (i * i) - omega * math.cos(f2) / i) / \
            PI * math.sin(smarg) / smarg
        coef[nfp - i] = - coef[nfp + i]
    coef[nfp] = 0.0

    for i in range(0, len(f) - 2 * nfp):
        sdfx = 0.0
        for j in range(0, 2 * nfp):
            sdfx = sdfx + coef[j] * f[i + j]
        dfx[nfp + i] = sdfx / dx
    return dfx
