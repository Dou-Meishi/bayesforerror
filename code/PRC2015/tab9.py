# -*- coding: utf-8 -*-
# Calculate PRC2015 TALBE. IX

import numpy as np
from numpy import sqrt
import prde, plpr
import logging.config, yaml     # logging module


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


def tab9(p, data, SET):
    T = data[0]
    Q = sqrt(.5*939*T)/Lambda
    ccck = _get_ccck(T, Lambda, data[1:])
    res = np.zeros(len(X_p))
    for n in range(len(X_p)):
        prde.setting_args(Q=Q, ccck=ccck[:n+1], **X_p[n])
        plpr.setting_args(SET=SET, Q=Q, **X_p[n], fmt=fmt)
        d_kw={'A': {'pr_delta' : prde.pr_delta_A,
                    'depth'    : 10,
                    'max_length': 300},
              'B': {'pr_delta' : prde.pr_delta_B,
                    'depth'    : 10,
                    'max_length': 300},
              'C': {'pr_delta' : prde.pr_delta_C,
                    'depth'    : 15,
                    'max_length': 1500}}
        res[n] = plpr.interpolate_dkp(p, **d_kw[SET])
    return res


def ini_log():
    with open('./logconf.yaml', 'r') as f:
        logconf = yaml.safe_load(f.read())
    logging.config.dictConfig(logconf)
    return


def main():
    ini_log()
    
    for p in [.68]:
        print('\n{:.0%} DOB:'.format(p))
        for data in dataIV:
            print('\nT: {}'.format(data[0]))
            print(res_fmt.format(*tab9(p, data, 'A')))
            print(res_fmt.format(*tab9(p, data, 'C')))
            print(res_fmt.format(*tab9(p, data, 'B')))

    logger = logging.getLogger('dmslog').getChild(__name__)
    logger.info("Main programe finished.")
    return


if __name__ == '__main__':
    import cProfile
    cProfile.run('main()', filename='tab9.out')
