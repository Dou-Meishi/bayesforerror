# -*- coding: utf-8 -*-
# Calculate PRC2015 TALBE. VIII

import numpy as np
from numpy import sqrt
from scipy.optimize import root_scalar


# first column is T_lab, not p_rel
# others are cross sections (See PRC2015 TABLE. IV)
Lambda = 600                    # R=0.9 fm
dataIV = [[50, 183.6, 166.5, 167.0, 166.8, 167.5],
          [96, 84.8, 75.1, 78.3, 77.5, 78.0],
          [143, 52.5, 49.1, 54.2, 53.7, 53.9],
          [200, 34.9, 35.9, 42.6, 43.2, 42.7]]


def _get_ccck(T, Lambda, sigmas):
    '''Extract ccck from experimental data (See Eq. (37)).'''
    Q = sqrt(.5*939*T)/Lambda   # 939 MeV is the mass of nucleon
    sigmas = [sigma/sigmas[0] for sigma in sigmas]
    ccck = [(sigmas[n]-sigmas[n-1])/Q**(n+1) for n in
            range(1, len(sigmas))]
    ccck.insert(0, 1.0)         # choose sigma_{LO} as refs
    ccck.insert(1, 0.0)         # c_1 is absent
    return ccck


def dkp_A_eps(p, *, Q, ccck, k, n_c, err=.001):
    '''Return the half width of p% DOB interval.'''
    cbar_ge = 1000
    cbar_k = np.amax(np.absolute(ccck))
    p_t = n_c/(n_c+1) * cbar_k
    p_t *= ((1/cbar_k)**(n_c+1)-(1/cbar_ge)**(n_c+1))
    p_t /= ((1/cbar_k)**(n_c)-(1/cbar_ge)**(n_c))
    if p <= p_t:
        return cbar_k/p_t * p * Q**(k+1)
    else:
        def f(dkp):
            y = 1/dkp**n_c - 1/(cbar_ge*Q**(k+1))**n_c
            y *= Q**((k+1)*(n_c+1))/n_c
            y -= (cbar_ge*Q**(k+1) - dkp)/cbar_ge**(n_c+1)
            y *= cbar_ge**n_c*cbar_k**n_c/(cbar_ge**n_c - cbar_k**n_c)
            y *= n_c/(n_c+1)/Q**(k+1)
            return 1-y-p
        return root_scalar(f, bracket=[cbar_k*Q**(k+1), cbar_ge*Q**(k+1)],
                           method='brentq').root


def main():
    LO_p = {'k': 0, 'n_c': 1}
    LO = {'k': 1, 'n_c': 1}     # c_1 is absent
    NLO = {'k': 2, 'n_c': 2}
    N2LO = {'k': 3, 'n_c': 3}
    N3LO = {'k': 4, 'n_c': 4}
    N4LO = {'k': 5, 'n_c': 5}
    X_p = [LO_p, LO, NLO, N2LO, N3LO, N4LO]

    args = {'err': .001}
    fmt = '{:<5.3g} {:<6.3g} {:<6.3g} {:<9.3g} {:<9.3g} {:<10.3g}'
    res = np.zeros(len(X_p))
    for p in [.68, .95]:
        print('\n{:.0%} DOB:'.format(p))
        for data in dataIV:
            sigma_ref = data[1]
            T = data[0]
            Q = sqrt(.5*939*T)/Lambda
            print('\nT: {:d}, Q: {:.3f},'.format(T, Q))
            ccck = _get_ccck(data[0], Lambda, data[1:])
            for n in range(len(X_p)):
                res[n] = dkp_A_eps(p, Q=Q, ccck=ccck[:n+1], **X_p[n])
            print(fmt.format(*res))
            print(fmt.format(*(res*sigma_ref)))
    return



if __name__ == '__main__':
    import cProfile
    cProfile.run('main()', filename='tab8.out')
