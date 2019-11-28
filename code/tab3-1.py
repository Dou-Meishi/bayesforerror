# -*- coding: utf-8 -*-
# Calculate table III

import numpy as np
from scipy.integrate import quad, nquad
from scipy.optimize import root_scalar

k = 2
cbar_le = .001
cbar_ge = 1000
limit = 50                      # arg passed to scipy.integrate.quad()

def find_width(p, pr, *, bound=1):
    '''
    Find a interval (-d,+d) on which the integral of pr(x) equals p.
    Assuming the p.d.f. pr(x) is symmetric about x=0.

    Return d.
    '''
    integral = lambda d: 2*quad(pr, 0, d)[0]
    D = integral(bound)
    t, max_t = (0, 1000)
    while p >= integral(bound) and t < max_t:
        t += 1
        bound *= 10
    if p < integral(bound):
        f = lambda d: p-integral(d)
        return root_scalar(f, bracket=[0, bound], method='brentq').root
    else:
        raise Exception('Upper bound of interval width not found!\
        \nThe last tried interval is (-{d},+{d}) on which the integral\
        is {prb} while given probability is {p}'.format(
            d=bound, prb=integral(bound), p=p))

def _dkp_eps(p, pr_cn_if_cbar, pr_cbar, *, Q, cns, limit=limit):
    def pr_delta(delta):
        def numerator_f(cbar):
            y = pr_cn_if_cbar(delta/Q**(k+1), cbar)
            for cn in cns:
                y *= pr_cn_if_cbar(cn, cbar)
            return y*pr_cbar(cbar)
        def denominator_f(cbar):
            y = Q**(k+1)
            for cn in cns:
                y *= pr_cn_if_cbar(cn, cbar)
            return y*pr_cbar(cbar)
        numerator, _ = quad(numerator_f, 0, np.inf, limit=limit)
        denominator, _ = quad(denominator_f, 0, np.inf, limit=limit)
        return numerator/denominator
    return find_width(p, pr_delta, bound=Q**(k+1))


def pr_cn_if_cbar(cn, cbar):
    return np.exp(-1/2 * cn**2/cbar**2)/np.sqrt(2*np.pi)/cbar

def pr_cbar(cbar):
    theta = lambda x:np.heaviside(x,0)
    return 1/cbar/np.log(cbar_ge/cbar_le) * theta(
            cbar-cbar_le)* theta(cbar_ge-cbar)

def main():
    for cns in [[1., 1., 1.], [1., .5, .1], [1., .1, .1]]:
        print('cns is [{:.2f},{:.2f},{:.2f}]:'.format(*cns))
        for p in [.68, .95]:
            print('\tDOB is {:.2%}:'.format(p))
            for Q in [.20, .33, .50]:
                print('\t\tQ is {:.2f}:'.format(Q))
                i3 = _dkp_eps(p, pr_cn_if_cbar, pr_cbar, Q=Q, cns=cns)
                print('{:.3g}'.format(i3))


if __name__ == '__main__':
    main()
