# -*- coding: utf-8 -*-
# test integrator in practice

import numpy as np
import matplotlib.pyplot as plt 
from numpy import pi, sin, cos, exp, log 
from scipy.integrate import nquad
from scipy.special import gamma as spgamma # factorial

from gaussquad import GaussQuad as GQ


def get_cos_sinc(delta, am):
    am = np.array(am)
    def cos_sinc(t):
        return cos(delta*t)*np.prod(np.sinc(t*am/pi))
    return cos_sinc


def spPrde(delta, cbar, Q, k, h):
    am = [cbar*Q**(k+1+m) for m in range(h)]
    i, _ = nquad(get_cos_sinc(delta,am), [[0,np.inf]], opts={'limit':200})
    return i/pi


def glPrde(delta, cbar, Q, k, h):
    am = [cbar*Q**(k+1+m) for m in range(h)]
    i, _ = GQ.QNGL(get_cos_sinc(delta,am), [0,np.inf], division=200)
    return i/pi


def gkPrde(delta, cbar, Q, k, h):
    am = [cbar*Q**(k+1+m) for m in range(h)]
    i, _ = GQ.QAGK(get_cos_sinc(delta,am), [0,np.inf], limit=200)
    return i/pi


def bwPrde(delta, cbar, Q, k, h):
    d = delta/(cbar*Q**(k+1))
    a = [Q**m for m in range(h)]
    
    if abs(d) > sum(a):
        return 0

    asum = 0
    s = [-1 for _ in range(h)]
    for _ in range(1<<h):
        b = np.dot(s, a)
        c = np.prod(s)

        asum += c*(b+d)**(h-1)*np.sign(b+d)

        for i in range(h):
            s[i] = -s[i]
            if s[i] == 1:
                break

    asum /= (2<<h) * spgamma(h) * np.prod(a) * cbar* Q**(k+1)

    return asum



def main():
    Q, k, h = (.511, 5, 10)
    
    cbars = [.01, .1, 1., 10., 100., 1000.]
    delta = np.linspace(0, .05*Q**(k+1))
    cbar = .01
    dlta = Q**(k+1)
    
    fig, ax = plt.subplots()
    ax.plot(delta, [spPrde(d, cbar, Q, k, h) for d in delta],
            'x', label='Scipy Quad')
    ax.plot(delta, [gkPrde(d, cbar, Q, k, h) for d in delta],
            'o', label='Gauss Kronord')
    ax.plot(delta, [glPrde(d, cbar, Q, k, h) for d in delta],
            '^', label='Gauss Legendre')
    ax.plot(delta, [bwPrde(d, cbar, Q, k, h) for d in delta],
            'd', label='Borwein Integral')
    ax.legend()
    plt.show()

    return


if __name__ == '__main__':
    main()
