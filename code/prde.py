# -*- coding: utf-8 -*-
# Filename: prde.py
# Calculate p.d.f. of Delta in TABLE. III
# for Set A, C without first-omitted-term approximation

import logging
import numpy as np
from numpy import sqrt, pi
from scipy.integrate import quad
from scipy.special import gamma, gammaincc
import c_priors as priors       # write with Cython

cbar_le = .001
cbar_ge = 1000
limit = 200                      # arg passed to scipy.integrate.quad()
_args = {}


def setting_args(*, Q, ccck, k=2, n_c=3, h=4, **kw):
    _args['Q'] = Q
    _args['ccck'] = ccck
    _args['h'] = h
    _args['k'] = k
    _args['n_c'] = n_c
    return


def _A_delta_if_cbar(delta, cbar):
    '''pr(Delta|cbar). See PRC2017 Eq. (A4).'''
    def f(t):
        return priors._A_delta_if_cbar_f(t,delta, cbar,
            _args['Q'], _args['h'], _args['k'])
    res, _ = quad(f, 0, np.inf, limit=limit)
    return res/np.pi


def pr_delta_A(delta):
    '''pr(Delta|ccck) for Set A. See PRC2017 Eq. (6).'''
    def _delta_if_cbar(cbar):
        return _A_delta_if_cbar(delta, cbar)
    def _ccck_if_cbar(cbar):
        y = 1
        for cn in _args['ccck']:
            y *= priors._pr_cn_if_cbar_A(cn, cbar)
        return y
    numerator, _ = quad(lambda cbar:_delta_if_cbar(
        cbar)*_ccck_if_cbar(cbar)*priors._pr_cbar_A(cbar),
                        cbar_le, cbar_ge, limit=limit)
    denominator, _ = quad(lambda cbar:_ccck_if_cbar(cbar
    )*priors._pr_cbar_A(cbar), cbar_le, cbar_ge, limit=limit)
    res = numerator/denominator
    logger = logging.getLogger('dmslog').getChild(__name__)
    logger.debug('delta: {}\npr: {}'.format(delta,res))
    return res


def pr_delta_C(delta):
    '''pr(Delta|ccck) for Set C. See PRC2017 Eq. (A9).'''
    Q, k, h, n_c = (_args['Q'], _args['k'], _args['h'], _args['n_c'])
    q = Q**(k+1) * sqrt((1-Q**(2*h))/(1-Q**2))
    ckb = np.sum(np.square(_args['ccck']))
    alpha = .5 * ckb
    beta = .5 * (ckb + delta*delta/(q*q))
    numerator = .5 * beta**(-n_c*.5-.5) * gamma(.5*n_c+.5) * (
        gammaincc(.5*n_c+.5, beta/(cbar_ge*cbar_ge)) -
        gammaincc(.5*n_c+.5, beta/(cbar_le*cbar_le)))
    denominator = .5 * alpha**(-n_c*.5) * gamma(.5*n_c) * (
        gammaincc(.5*n_c,alpha/(cbar_ge*cbar_ge)) -
        gammaincc(.5*n_c,alpha/(cbar_le*cbar_le)))
    return numerator/(denominator*q*sqrt(2*pi))


def main():
    args = {'SET': 'Cc', 'Q': .5}
    import plpr
    setting_args(**args)
    plpr.setting_args(**args)
    # plpr.calc(pr_delta_C, 'refine')
    print(plpr.find_width(.95, err=.001))
    plpr.plot_pr(show=True)
    return


if __name__ == '__main__':
    import cProfile
    cProfile.run('main()', filename='prde.out')
