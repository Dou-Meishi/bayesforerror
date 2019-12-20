# -*- coding: utf-8 -*-
# Filename: prde.py
# Calculate p.d.f. of Delta in TABLE. III
# for Set A, C without first-omitted-term approximation

import logging
import numpy as np
from numpy import sqrt, pi, exp, log
from scipy.integrate import quad, nquad
from scipy.special import gamma, gammaincc, erf, erfc
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


def pr_delta_A(delta):
    '''pr(Delta|ccck) for Set A. See PRC2017 Eq. (6).'''
    def _delta_if_cbar(cbar):
        res, _ = nquad(priors._A_delta_if_cbar_f,
                      [[0, np.inf]],
                      args=[delta, cbar,
                            _args['Q'], _args['h'], _args['k']],
                      opts={'limit': limit})
        return res/np.pi
    ccck = np.array(_args['ccck'])
    cbar_k = np.amax(np.abs(ccck[ccck.nonzero()]))

    numerator, _ = quad(lambda cbar:_delta_if_cbar(
        cbar)/cbar**(_args['n_c']+1), cbar_k, cbar_ge, limit=limit)
    denominator = 1/cbar_k**_args['n_c'] - 1/cbar_ge**_args['n_c']
    res = numerator/denominator*_args['n_c']
    logger = logging.getLogger('dmslog').getChild(__name__)
    logger.debug('delta: {}\npr: {}'.format(delta,res))
    return res


def pr_delta_B(delta):
    '''pr(Delta|ccck) for Set B. See PRC2017 Eq. (A9).'''
    def _delta_if_cbar(cbar):
        res, _ = nquad(priors._A_delta_if_cbar_f, # same as Set A
                       [[0, np.inf]],
                       args=[delta, cbar,
                             _args['Q'], _args['h'], _args['k']],
                       opts={'limit': limit})
        return res/np.pi
    ccck = np.array(_args['ccck'])
    cbar_k = np.amax(np.abs(ccck[ccck.nonzero()]))
    numerator, _ = quad(lambda cbar:_delta_if_cbar(
        cbar)*priors._weight(cbar,_args['n_c']),
                        cbar_k, cbar_ge, limit=limit)
    denominator = (erfc((log(cbar_k)+_args['n_c'])/sqrt(2)) -
                   erfc((log(cbar_ge)+_args['n_c'])/sqrt(2)))
    res = numerator/denominator*sqrt(2/pi)*exp(
        -.5*_args['n_c']*_args['n_c'])
    logger = logging.getLogger('dmslog').getChild(__name__)
    logger.debug('delta: {}\npr: {}'.format(delta, res))
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


def test_denominator(ccck, delta, n_c):
    ccck = np.array(ccck)
    def _ccck_if_cbar(cbar):
        y = 1
        for cn in ccck[ccck.nonzero()]:
            y *= priors._pr_cn_if_cbar_B(cn, cbar)
        return y
    denominator, _ = quad(lambda cbar:_ccck_if_cbar(cbar
    )*priors._pr_cbar_B(cbar), cbar_le, cbar_ge, limit=limit)
    print(denominator)

    cbar_k = np.amax(np.abs(ccck[ccck.nonzero()]))
    denominator, _ = quad(lambda t: sqrt(2)*exp(-t**2+.5*n_c*n_c),
                          (n_c+log(cbar_k))/sqrt(2),
                          (n_c+log(cbar_ge))/sqrt(2),
                          limit=limit)
    denominator /= 2**n_c * sqrt(2*pi)
    print(denominator)

    denominator = exp(.5*n_c*n_c)/2**(n_c+1)*(
        erfc((log(cbar_k)+n_c)/sqrt(2)) -
        erfc((log(cbar_ge)+n_c)/sqrt(2)))
    print(denominator)
    return


def main():
    Q = .511
    k = 5
    n_c = 5
    ccck = [1., .0, .11, 1.44, .25, -.41]
    delta = Q**(k+1)
    test_denominator(ccck, delta, n_c)
    return


if __name__ == '__main__':
    import cProfile
    cProfile.run('main()', filename='prde.out')
