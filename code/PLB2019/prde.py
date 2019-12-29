# -*- coding: utf-8 -*-

import numpy as np
from numpy import pi, sqrt, exp, log
from scipy.integrate import nquad, quad
from scipy.special import gamma, gammaincc

import priors

cbar_le, cbar_ge = (.001, 1000)

def prde_AB_cbar(delta, cbar, Q, h, k):
    i, _ = nquad(priors._A_delta_if_cbar_f, [[0, np.inf]],
                 args=[delta, cbar, Q, h, k], opts={'limit':200})
    return i/pi


def prde_A(delta, ccck, Q, h, k):
    def weight(cbar):
        res = priors._pr_cbar_A(cbar)
        for cn in ccck:
            res *= priors._pr_cn_if_cbar_A(cn, cbar)
        return res

    lower = np.amax(np.abs(ccck))
    numerator, _ = quad(lambda cbar: prde_AB_cbar(delta,cbar,Q,h,k)*(
        weight(cbar)), lower, cbar_ge)
    denominator, _ = quad(weight, lower, cbar_ge)
    return numerator/denominator


def prde_B(delta, ccck, Q, h, k):
    def weight(cbar):
        res = priors._pr_cbar_B(cbar)
        for cn in ccck:
            res *= priors._pr_cn_if_cbar_B(cn, cbar)
        return res

    lower = np.amax(np.abs(ccck))
    numerator, _ = quad(lambda cbar: prde_AB_cbar(delta,cbar,Q,h,k)*(
        weight(cbar)), lower, np.inf)
    denominator, _ = quad(weight, lower, np.inf)
    return numerator/denominator


def prde_C(delta, ccck, Q, h, k):
    n_c = len(ccck)
    q = Q**(k+1) * sqrt((1-Q**(2*h))/(1-Q**2))
    ckb = np.sum(np.square(ccck))
    alpha = .5 * ckb
    beta = .5 * (ckb + delta*delta/(q*q))
    numerator = .5 * beta**(-n_c*.5-.5) * gamma(.5*n_c+.5) * (
        gammaincc(.5*n_c+.5, beta/(cbar_ge*cbar_ge)) -
        gammaincc(.5*n_c+.5, beta/(cbar_le*cbar_le)))
    denominator = .5 * alpha**(-n_c*.5) * gamma(.5*n_c) * (
        gammaincc(.5*n_c,alpha/(cbar_ge*cbar_ge)) -
        gammaincc(.5*n_c,alpha/(cbar_le*cbar_le)))
    return numerator/(denominator*q*sqrt(2*pi))
    

prde = {'A': prde_A,
        'B': prde_B,
        'C': prde_C}


def test_prde():

    Q, h, k = (.5, 10, 2)
    ccck = [1., .5, .1]

    import matplotlib.pyplot as plt

    deltas = np.linspace(0, 2*Q**(k+1))
    A = [prde_A(delta, ccck[1:], Q, h, k) for delta in deltas]
    B = [prde_B(delta, ccck[1:], Q, h, k) for delta in deltas]
    C = [prde_C(delta, ccck[1:], Q, h, k) for delta in deltas]


    plt.plot(deltas, A, 'o', label='Set A')
    plt.plot(deltas, B, 'o', label='Set B')
    plt.plot(deltas, C, 'o', label='Set C')
    
    plt.legend()
    plt.show()
    return



if __name__ == '__main__':
    import cProfile
    cProfile.run('test_prde()', filename='prde.out')
