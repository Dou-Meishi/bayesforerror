# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt 
from numpy import pi, exp, sqrt
from scipy.integrate import nquad
from scipy.optimize import root_scalar

class PDF(object):
    '''p.d.f., symmetric about x=0.'''

    path = './data'
    
    def __init__(self, func, basename=None):
        self.func = func
        self.basename = basename


    def __call__(self, *args):
        return self.func(*args)


    def get_dkp_bisect(self, p, *, args, upper):
        def f(d):
            i, _ = nquad(self.func, [[0,d]], args=args)
            return i*2-p
        
        for _ in range(20):
            f_b = f(upper)
            upper *= 2
            if f_b > 0:
                return root_scalar(f, method='brentq',
                                   bracket=[0, upper/2]).root
        else:
            raise RuntimeError("Upper bound not found!")


    def get_dkp_newton(self, p, *, args, x0):
        def f(d):
            i, _ = nquad(self.func, [[0,d]], args=args)
            return i*2-p
        def fp(d):
            return 2*self.func(d, *args)

        return root_scalar(f, method='newton', x0=x0, fprime=fp).root


def test():
    prde = PDF(lambda t, sigma: exp(-.5*t*t/(sigma*sigma))/(
        sigma*sqrt(2*pi)))
    print(prde.get_dkp_bisect(.6826, args=[1.], upper=1.))
    print(prde.get_dkp_newton(.6826, args=[1.], x0=1.))

    print(prde.get_dkp_bisect(.9544, args=[1.], upper=1.))
    print(prde.get_dkp_newton(.9544, args=[1.], x0=1.))

    print(prde.get_dkp_bisect(.9974, args=[1.], upper=1.))
    print(prde.get_dkp_newton(.9974, args=[1.], x0=1.))    
    
    return
    


if __name__ == '__main__':
    import cProfile
    cProfile.run('test()', filename='dmsPDF.out')
