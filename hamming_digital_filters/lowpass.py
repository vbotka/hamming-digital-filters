# All rights reserved (c) 2022, Vladimir Botka <vbotka@gmail.com>
# Simplified BSD License, https://opensource.org/licenses/BSD-2-Clause

# Digital Filters by R.W.Hamming,third edition
# Dover publications, INC., Mineola, New York
# ISBN 0-486-65088-X

import math


def lowpass_nr(f, fc=0.2, nfp=5):
    ''' Chapter 6.2 A Low-Pass Filter Design '''

    coef = [0.0] * (2 * nfp + 1)
    coef[nfp] = 2.0 * fc
    sfx = [0.0] * len(f)
    for k in range(1, nfp):
        xd = math.pi * k
        yd = xd / nfp
        c1 = math.sin(yd) / yd
        c2 = math.sin(2.0 * fc * xd) / xd
        coef[nfp + k] = c1 * c2
        coef[nfp - k] = c1 * c2

    for i in range(0, len(f) - 2 * nfp):
        ssfx = 0.0
        for j in range(0, 2 * nfp):
            ssfx = ssfx + coef[j] * f[i + j]
        sfx[nfp + i] = ssfx
    return sfx
