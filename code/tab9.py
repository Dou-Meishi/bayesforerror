# -*- coding: utf-8 -*-
# Calculate PRC2015 TALBE. IX

import numpy as np
from numpy import sqrt
import prde, plpr


# first column is T_lab, not p_rel
# others are cross sections (See PRC2015 TABLE. IV)
Lambda = 600                    # R=0.9 fm
dataIV = [[50, 183.6, 166.5, 167.0, 166.8, 167.5],
          [96, 84.8, 75.1, 78.3, 77.5, 78.0],
          [143, 52.5, 49.1, 54.2, 53.7, 53.9],
          [200, 34.9, 35.9, 42.6, 43.2, 42.7]]
X_p = [dict(zip(['k', 'n_c'], [0, 1])), # LO_prime
       dict(zip(['k', 'n_c'], [1, 1])), # LO
       dict(zip(['k', 'n_c'], [2, 2])), # NLO
       dict(zip(['k', 'n_c'], [3, 3])), # N2LO
       dict(zip(['k', 'n_c'], [4, 4])), # N3LO
       dict(zip(['k', 'n_c'], [5, 5]))] # N4LO

fmt = '{SET}_Q-{Q:.3g}-k-{k:d}'
res_fmt = '{:<5.3g} {:<6.3g} {:<6.3g} {:<9.3g} {:<9.3g} {:<10.3g}'


def _get_ccck(T, Lambda, sigmas):
    '''Extract ccck from experimental data (See Eq. (37)).'''
    Q = sqrt(.5*939*T)/Lambda   # 939 MeV is the mass of nucleon
    sigmas = [sigma/sigmas[0] for sigma in sigmas]
    ccck = [(sigmas[n]-sigmas[n-1])/Q**(n+1) for n in
            range(1, len(sigmas))]
    ccck.insert(0, 1.0)         # choose sigma_{LO} as refs
    ccck.insert(1, 0.0)         # c_1 is absent
    return ccck


def tab9C(p, data, *, err=.001):
    T = data[0]
    Q = sqrt(.5*939*T)/Lambda
    ccck = _get_ccck(T, Lambda, data[1:])
    res = np.zeros(len(X_p))
    depth = 10
    for n in range(len(X_p)):
        prde.setting_args(Q=Q, ccck=ccck[:n+1], **X_p[n])
        plpr.setting_args(SET='C', Q=Q, **X_p[n], fmt=fmt)
        res[n] = plpr.get_dkp(p, depth=10, err=err,
                              pr_delta=prde.pr_delta_C)
    return res


def main():
    # p = .68
    # data = dataIV[0]
    # T = data[0]
    # Q = sqrt(.5*939*T)/Lambda
    # ccck = _get_ccck(T, Lambda, data[1:])
    # prde.setting_args(Q=Q, ccck=ccck, k=5, n_c=5)
    # plpr.setting_args(SET='C', Q=Q, k=5, fmt=fmt)
    # plpr.calc(prde.pr_delta_C, 'refine')
    # print(plpr.find_width(p))
    # plpr.plot_pr(show=True)
    for p in [.68, .95]:
        print('\n{:.0%} DOB:'.format(p))
        for data in dataIV:
            print('\nT: {}'.format(data[0]))
            print(res_fmt.format(*tab9C(p, data)))
    return


if __name__ == '__main__':
    import cProfile
    cProfile.run('main()', filename='tab9.out')
