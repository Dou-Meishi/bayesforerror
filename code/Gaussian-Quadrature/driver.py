# -*- coding: utf-8 -*-
'''Test performance.'''


import numpy as np
import time, functools
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

    fmtstr = '\n'.join(["     integrand: {func_name}",
               "     interval: [{interval[0]}, {interval[1]}]",
               "    quad rule: {quad_name}",
               " n evaluation: {neval:9}     (total of {repeat} times)",
               "time per call: {tmin:9.6f} sec (best of {repeat} times)",
               "integral, err: {i:9.6f}, {e:9.6f}",
               "          ref: {ii:9.6f}, {ee:9.6f}",
               "     add info: {add info}"])


    def __init__(self, integrator):
        self.integrator = integrator
        self.prb = Probe(repeat=3)


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
                 'add info': np.nan}


    def test_x3(self):
        def xto3(x): return x**3
        odict = self._test(xto3, [0,1])
        odict['ii'] = .25
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



class InteSincDriver(InteDriver):

    @staticmethod
    def get_cos_sinc(delta, am):
        am = np.array(am)
        def cos_sinc(t):
            return cos(delta*t) * np.prod(np.sinc(t*am/pi))
        return cos_sinc


    def test_0(self):
        delta, am = (0, [1., 1.])
        odict = self._test(self.get_cos_sinc(delta, am), [0,np.inf])
        odict['ii'] = pi/2
        odict['ee'] = abs(odict['i'] - odict['ii'])
        odict['add info'] = "delta: {:.3f} am: {}".format(delta,str(am))
        return self.fmtstr.format(**odict)

        
    def test_1(self):
        delta, am = (0, [1., 1/3, 1/5, 1/7])
        odict = self._test(self.get_cos_sinc(delta, am), [0,np.inf])
        odict['ii'] = pi/2
        odict['ee'] = abs(odict['i'] - odict['ii'])
        odict['add info'] = "delta: {:.3f} am: {}".format(delta,str(am))
        return self.fmtstr.format(**odict)



def main():
    from scipy.integrate import quad, nquad
    def mquad(func, interval, *args, **kws):
        return nquad(func, [interval], *args, **kws)
    

    mq_driver = InteSincDriver(mquad)

    for f in dir(mq_driver):
        if not f.startswith('test_'):
            continue

        func = getattr(mq_driver, f)
        if callable(func):
            print('\n', func())
        
    return



if __name__ == '__main__':
    main()
