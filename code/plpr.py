# -*- coding: utf-8 -*-
# Filename: plpr.py
# plot p.d.f. of Delta (without first-term-approximation)
# for SET A or C with given function pr_delta()
# and find p*100% DOB interval width using that p.d.f.

import os
import numpy as np
from scipy.integrate import trapz
from matplotlib import pyplot as plt

_args = {}
fmt = '{}_Q-{:3g}_{}.npy'

def setting_args(*,  SET, Q, datapath='data', k=2, **kw):
    _args['datapath'] = datapath
    _args['fn_x'] = os.path.join(datapath, fmt.format(SET, Q, 'x'))
    _args['fn_y'] = os.path.join(datapath, fmt.format(SET, Q, 'y'))
    _args['SET'] = SET
    _args['Q'] = Q
    _args['k'] = k
    return 

def calc(pr_delta, mode='relax'):
    try:
        print('Getting previous calculated data...')
        x = np.load(_args['fn_x'])
        y = np.load(_args['fn_y'])
    except FileNotFoundError:
        print('Previous data not found! Initialize...')
        x = np.linspace(0,_args['Q']**(_args['k']+1),10)
        y = [pr_delta(xx) for xx in x]

    if mode == 'refine':
        step = (x[1] - x[0])/2
        if abs(x[2]-step*2-x[1]) > .001*step:
            raise ValueError('Please check stepsize of x!')
        else:
            print('Calculating new {} points...'.format(len(x)))
        new_x = x + step
        new_y = [pr_delta(xx) for xx in new_x]
        x = np.stack((x, new_x),axis=-1).reshape(-1)[:-1]
        y = np.stack((y, new_y),axis=-1).reshape(-1)[:-1]
    elif mode == 'extend':
        step = x[-1] - x[0] + x[1] - x[0]
        if abs(x[1] + step - x[-1]) < .999 * (x[1]-x[0]):
            raise ValueError('Please check stepsize of x!')
        else:
            print('Calculating new {} points...'.format(len(x)))
        new_x = x + step
        new_y = [pr_delta(xx) for xx in new_x]
        x = np.hstack((x, new_x))
        y = np.hstack((y, new_y))
    else:
        pass

    if os.path.exists(_args['datapath']):
        print('Saving calculated data...')
    else:
        print('Data dir not found! Create it and saving data...')
        os.mkdir(_args['datapath'])

    np.save(_args['fn_x'], x)
    np.save(_args['fn_y'], y)
    return x, y

def find_width(p, *, err=.005):
    x = np.load(_args['fn_x'])
    y = np.load(_args['fn_y'])

    if 2*trapz(y,x) < p:
        raise ValueError('Need more (extend) data points!\n\
The integration on ({},{}) is {} while Given probability\
is {}.'.format(x[0],x[-1],2*trapz(y,x),p))
    else:
        D, d = (0, 0)
    for i in range(len(x)-1):
        d = 2*trapz(y[i:i+2], x[i:i+2])
        D += d
        if D > p:
            if d < err:
                return x[i], x[i+1]
            else:
                pass
                raise ValueError('Need more (refine) data points!\n\
Upper error is {} while required accuracy is {}.'.format(d, err))
        else:
            pass
    return -1                     # won't happen usually

def plot_pr():
    x = np.load(_args['fn_x'])
    y = np.load(_args['fn_y'])
    plt.scatter(x, y)
    plt.show()
    return 'Integration:', 2*trapz(y, x)

def _main():
    setting_args(Q=.33, SET='A')
    # calc(pr_delta, 'refine')
    print(find_width(.95))
    print(plot_pr())

if __name__ == '__main__':
    import cProfile
    cProfile.run('_main()', filename='plpr.out', sort='cumulative')
