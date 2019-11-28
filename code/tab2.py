# -*- coding: utf-8 -*-
# Calculate table II

# todo:
# + round off rule as uncertainty

import numpy as np
from scipy.optimize import root_scalar

def dkp_A_eps(p, *, Q, k, cbar_ge, cbar_k, err=.001):
    '''Return the half width of p% DOB interval.'''
    n_c = int(k)+1
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
    args = [[.20, .33, .50],       # Q
            [1000, 4., 2.],        # cbar_ge
            [0, 1, 2, 3, 4],       # k
            [1]]                   # cbar_k
    # get the Cartesian product of args
    args = np.stack(np.meshgrid(*args),axis=-1).reshape(-1,len(args))
    keys = ['Q','cbar_ge','k','cbar_k']

    for p in [.68,.95]:
        print('{:.0%} DOB width:'.format(p))
        for arg in args:
            kw = dict(zip(keys,arg))
            d = dkp_A_eps(.68,**kw)
            print('\t{:.3g}'.format(d))


if __name__ == '__main__':
    main()
