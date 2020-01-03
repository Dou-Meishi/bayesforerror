# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from numpy import pi, exp, sqrt
from scipy.integrate import nquad, quad
from scipy.optimize import root_scalar

import prde                     # define pr_delta

class PDF(object):
    '''single variable p.d.f., symmetric about x=0.'''

    dirpath = './data'
    
    def __init__(self, func, basename=None):
        self.func = func
        self.basename = basename # basename of data file


    def __call__(self, t, *args, **kw):
        return self.func(t, *args, **kw)


    def get_dkp_bisect(self, p, *, upper, args=[], kw={}):
        def f(d):
            return 2*quad(lambda t: self.func(t,*args,**kw), 0, d)[0] - p

        for _ in range(20):
            upper *= 2
            if f(upper*.5) > 0:
                return root_scalar(f, method='brentq',
                                   bracket=[0, upper*.5]).root
        else:
            raise RuntimeError("Upper bound not found!")


    def get_dkp_newton(self, p, *, x0, args=[], kw={}):
        def f(d):
            return 2*quad(lambda t: self.func(t,*args,**kw), 0, d)[0] - p
        def fp(d): return 2*self.func(d, *args, **kw)

        return root_scalar(f, method='newton', x0=x0, fprime=fp).root


def test_PDF():
    prde = PDF(lambda t, sigma: exp(-.5*t*t/(sigma*sigma))/(
        sigma*sqrt(2*pi)))
    print(prde.get_dkp_bisect(.6826, args=[1.], upper=1.))
    print(prde.get_dkp_newton(.6826, args=[1.], x0=1.))

    print(prde.get_dkp_bisect(.9544, args=[1.], upper=1.))
    print(prde.get_dkp_newton(.9544, args=[1.], x0=1.))

    print(prde.get_dkp_bisect(.9974, args=[1.], upper=1.))
    print(prde.get_dkp_newton(.9974, args=[1.], x0=1.))
    
    return


class PrDelta(PDF):
    fmt = "{SET}_Q-{Q}_k-{k}_h-{h}"

    def __init__(self, SET, Q, c0ck, abcent_id, h=10):
        self.SET = SET
        self.Q = Q
        self.c0ck = c0ck.copy()
        self.abcent_id = abcent_id.copy()
        self.h = h

        self.k = len(c0ck) - 1
        self.ccck = [c0ck[i] for i in range(self.k+1)
                     if i not in abcent_id]
        self.n_c = len(self.ccck)

        super().__init__(lambda t: prde.prde[SET](t, self.ccck,
                                                  self.Q, self.h, self.k),
                         self.fmt.format(SET=SET, Q=self.Q,
                                         k=self.k, h=self.h))


    def update_ck(self, ck1, abcent=False):
        c0ck = self.c0ck.copy()
        abcent_id = self.abcent_id.copy()

        c0ck.append(ck1)
        if abcent:
            abcent_id.append(self.k+1)
        return PrDelta(self.SET, self.Q, c0ck, abcent_id, self.h)


    def get_dkp(self, p, *, method):
        if method == 'bisect':
            return self.get_dkp_bisect(p, upper=self.Q**(self.k + 1))
        elif method == 'newton':
            return self.get_dkp_newton(p, x0=self.Q**(self.k +1))
        else:
            raise KeyError(
                "This method '{}' is not available.".format(method))
        return
        

def test_PrDelta():
    A = PrDelta(SET='A', Q=.5, c0ck=[1.,.5,.1], abcent_id=[0])
    B = PrDelta(SET='B', Q=.5, c0ck=[1.,.5,.1], abcent_id=[0])
    C = PrDelta(SET='C', Q=.5, c0ck=[1.,.5,.1], abcent_id=[0])

    deltasA = np.linspace(0, 2*A.Q**(A.k+1))
    deltasB = np.linspace(0, 2*B.Q**(B.k+1))
    deltasC = np.linspace(0, 2*C.Q**(C.k+1))

    plt.plot(deltasA, [A(delta) for delta in deltasA], 'o', label=A.SET)
    plt.plot(deltasB, [B(delta) for delta in deltasB], 'o', label=B.SET)
    plt.plot(deltasC, [C(delta) for delta in deltasC], 'o', label=C.SET)

    plt.legend()
    plt.show()

    # A = PrDelta(SET='A', Q=.511, c0ck=[1., .0, .11, 1.44, .25, -.41],
    #             abcent_id=[1], h=4)
    # B = PrDelta(SET='B', Q=.511, c0ck=[1., .0, .11, 1.44, .25, -.41],
    #             abcent_id=[1], h=4)
    # C = PrDelta(SET='C', Q=.511, c0ck=[1., .0, .11, 1.44, .25, -.41],
    #             abcent_id=[1], h=4)

    # print(A.basename, A.get_dkp(.95, method='newton'))
    # print(B.basename, B.get_dkp(.95, method='newton'))
    # print(C.basename, C.get_dkp(.95, method='newton'))
    return



if __name__ == '__main__':
    import cProfile
    cProfile.run('test_PrDelta()', filename='dmsPDF.out')
