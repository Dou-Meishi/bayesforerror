# -*- coding: utf-8 -*-
'''Test performance.'''


import numpy as np
import time, functools
from math import factorial
from numpy import pi, sin, cos, exp, log



class Probe(object):
    '''a collection of decorators.'''

    def __init__(self, repeat=3, start=0):
        self.repeat= repeat    # timer
        self.count = start     # counter
        self.tlist = np.zeros(repeat)


    def counter(self, func, recount=True):
        '''counter wrapper.'''
        if recount:
            self.count = 0
        @functools.wraps(func)
        def new_func(*args, **kws):
            self.count += 1
            return func(*args, **kws)
        return new_func


    def timer(self, func):
        '''timer wrapper.'''
        @functools.wraps(func)
        def new_func(*args, **kws):
            for i in range(self.repeat):
                ts = time.time()
                res = func(*args, **kws)
                self.tlist[i] = time.time() - ts
            return res
        return new_func



class InteDriver(object):

    fmtstr = '\n'.join(["   integrand: {func_name}",
               "     interval: [{interval[0]}, {interval[1]}]",
               "    quad rule: {quad_name}",
               " n evaluation: {neval:9}     (total of {repeat} times)",
               "time per call: {tmin:9.8f} sec (best of {repeat} times)",
               "integral, err: {i:9.8f}, {e:9.8f}",
               "          ref: {ii:9.8f}, {ee:9.8f}",
               "     add info: {add info}"])


    def __init__(self, integrator, repeat=3):
        self.integrator = integrator
        self.prb = Probe(repeat=repeat)


    def __call__(self, *args, **kws):
        return self.integrator(*args, **kws)


    def _test(self, func, interval, *args, **kws):
        new_func = self.prb.counter(func)
        new_inte = self.prb.timer(self.integrator)
        i, e = new_inte(new_func, interval, *args, **kws)
        return {'func_name': func.__name__,
                 'interval': interval,
                'quad_name': self.integrator.__name__,
                    'neval': self.prb.count,
                     'tmin': np.amin(self.prb.tlist),
                   'repeat': self.prb.repeat,
                'i': i, 'e': e, 'ii': np.nan, 'ee': np.nan,
                 'add info': None}


    def test_x29(self):
        def xto29(x): return x**29
        odict = self._test(xto29, [0,1])
        odict['ii'] = .1/3
        odict['ee'] = abs(odict['i'] - odict['ii'])
        return self.fmtstr.format(**odict)


    def test_sin(self):
        odict = self._test(sin, [0, pi/2])
        odict['ii'] = 1
        odict['ee'] = abs(odict['i'] - odict['ii'])
        return self.fmtstr.format(**odict)


    def test_exp(self):
        odict = self._test(exp, [-np.inf, 0])
        odict['ii'] = 1
        odict['ee'] = abs(odict['i'] - odict['ii'])
        return self.fmtstr.format(**odict)


    def runtest(self):
        for attr in dir(self):
            if not attr.startswith('test_'):
                continue

            f = getattr(self, attr)
            if callable(f):
                print('\n', f())

        return



class InteSincDriver(InteDriver):

    @staticmethod
    def get_cos_sinc(delta, am):
        am = np.array(am)
        def cos_sinc(t):
            return cos(delta*t) * np.prod(np.sinc(t*am/pi))
        return cos_sinc

    @staticmethod
    def borw_cos_sinc(delta, am):
        d = delta
        a = np.array(am)
        h = len(a)

        if abs(d) > sum(a):
            return 0
        else:
            pass

        asum = 0
        s = [-1 for _ in range(h)]
        for _ in range(1<<h):
            c = np.prod(s)
            b = np.dot(s, a)
            asum += c*(b+d)**(h-1)*np.sign(b+d)

            for i in range(h):
                s[i] = -s[i]    # update s
                if s[i] == 1:
                    break

        asum /= (2<<h) * factorial(h-1) * np.prod(a)
        return asum*pi
            

    def test_0(self):
        delta, am = (0, [1., 1.])
        odict = self._test(self.get_cos_sinc(delta, am), [0,np.inf])
        odict['ii'] = self.borw_cos_sinc(delta,am)
        odict['ee'] = abs(odict['i'] - odict['ii'])
        odict['add info'] = "delta: {:.3f} am: {}".format(delta,str(am))
        return self.fmtstr.format(**odict)

        
    def test_1(self):
        delta, am = (0, [1., 1/2, 1/4, 1/8, 1/16, 1/32])
        odict = self._test(self.get_cos_sinc(delta, am), [0,np.inf])
        odict['ii'] = self.borw_cos_sinc(delta,am)
        odict['ee'] = abs(odict['i'] - odict['ii'])
        odict['add info'] = "delta: {:.3f} am: {}".format(delta,str(am))
        return self.fmtstr.format(**odict)


    def test_2(self):
        delta, am = (np.random.rand(), np.random.rand(10)+.001)
        odict = self._test(self.get_cos_sinc(delta, am), [0,np.inf])
        odict['ii'] = self.borw_cos_sinc(delta,am)
        odict['ee'] = abs(odict['i'] - odict['ii'])
        odict['add info'] = "delta: {:.3f} am: {}".format(delta, str(am))
        return self.fmtstr.format(**odict)


    def test_3(self, cbar_scale=10):
        delta = np.random.rand()*2*pi
        cbar  = np.random.rand()*cbar_scale + .01
        Q     = np.random.rand()*.5 + .5
        k     = int(np.random.rand()*4+2)
        h     = 10
        am = [cbar*Q**(k+1+m) for m in range(h)]
        odict = self._test(self.get_cos_sinc(delta, am), [0,np.inf])

        odict['ii'] = self.borw_cos_sinc(delta/(cbar*Q**(k+1)), [Q**m for m in range(h)])/(cbar*Q**(k+1))
        odict['ee'] = abs(odict['i'] - odict['ii'])
        odict['add info'] = "delta: {:.3f} cbar: {:.3f} Q: {:.3f} k: {} am: {}".format(delta, cbar, Q, k, str(am))
        return self.fmtstr.format(**odict)




def main():
    from scipy.integrate import quad, nquad
    def mquad(func, interval, *args, **kws):
        return nquad(func, [interval], *args, **kws, opts={'limit':50})
    
    mq_driver = InteSincDriver(mquad)
    mq_driver.runtest()
    return



if __name__ == '__main__':
    main()
