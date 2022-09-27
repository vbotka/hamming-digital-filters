# All rights reserved (c) 2022, Vladimir Botka <vbotka@gmail.com>
# Simplified BSD License, https://opensource.org/licenses/BSD-2-Clause

# Digital Filters by R.W.Hamming,third edition
# Dover publications, INC., Mineola, New York
# ISBN 0-486-65088-X

import math


def gradient_nr(f, dx, fc=0.2, nfp=5):
    ''' Chapter 6.4 Differentiation filter nonrecursive '''

    omega = 2.0 * math.pi * fc
    f3 = nfp + 1
    coef = [0.0] * (2 * nfp + 1)
    coef[nfp] = 0.0
    dfx = [0.0] * len(f)
    for k in range(1, nfp):
        f2 = k * omega
        smarg = math.pi * k / f3
        coef[nfp + k] = (math.sin(f2) / (k * k) - omega * math.cos(f2) / k) / \
            math.pi * math.sin(smarg) / smarg
        coef[nfp - k] = - coef[nfp + k]

    for i in range(0, len(f) - 2 * nfp):
        sdfx = 0.0
        for j in range(0, 2 * nfp):
            sdfx = sdfx + coef[j] * f[i + j]
        dfx[nfp + i] = sdfx / dx
    return dfx
