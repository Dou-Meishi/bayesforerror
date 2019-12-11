# -*- coding: utf-8 -*-
# Filename: plpr.py
# plot p.d.f. of Delta at specific points
# for SET A or C with given function pr_delta()
# and read the p*100% DOB interval width from figures

import os, logging
import numpy as np
from scipy.integrate import trapz
from matplotlib import pyplot as plt

_args = {}


def setting_args(*,  SET, Q, k=2, fmt, datapath='data', **kw):
    _args['datapath'] = datapath
    _args['fn'] = os.path.join(datapath, fmt.format(
        SET=SET, Q=Q, k=k, **kw))
    _args['fn_x'] = _args['fn'] + '_x.npy'
    _args['fn_y'] = _args['fn'] + '_y.npy'
    _args['SET'] = SET
    _args['Q'] = Q
    _args['k'] = k
    return


def calc(pr_delta, mode='relax'):
    logger = logging.getLogger('dmslog').getChild(__name__)
    try:
        logger.info('Getting previous calculated data...')
        x = np.load(_args['fn_x'])
        y = np.load(_args['fn_y'])
    except FileNotFoundError:
        logger.info('Previous data not found! Initialize...')
        x = np.linspace(0,_args['Q']**(_args['k']+1),10)
        y = [pr_delta(xx) for xx in x]

    if mode == 'refine':
        step = (x[1] - x[0])/2
        logger.info('Calculating new {} points...'.format(len(x)))
        new_x = x + step
        new_y = [pr_delta(xx) for xx in new_x]
        x = np.stack((x, new_x),axis=-1).reshape(-1)[:-1]
        y = np.stack((y, new_y),axis=-1).reshape(-1)[:-1]
    elif mode == 'extend':
        step = x[-1] - x[0] + x[1] - x[0]
        logger.info('Calculating new {} points...'.format(len(x)))
        new_x = x + step
        new_y = [pr_delta(xx) for xx in new_x]
        x = np.hstack((x, new_x))
        y = np.hstack((y, new_y))
    else:
        pass

    if os.path.exists(_args['datapath']):
        logger.info('Saving calculated data...')
    else:
        logger.warning('Data dir not found! Create it and saving data...')
        os.mkdir(_args['datapath'])

    np.save(_args['fn_x'], x)
    np.save(_args['fn_y'], y)
    return x, y


class ThinIntervalError(Exception):
    pass



class LowAccuracyError(Exception):
    pass



def find_width(p, *, err=.005):
    x = np.load(_args['fn_x'])
    y = np.load(_args['fn_y'])

    if 2*trapz(y,x) < p:
        raise ThinIntervalError('Need more (extend) data points!\n\
 The integration on ({},{}) is {} while Given probability\
 is {}.'.format(x[0],x[-1],2*trapz(y,x),p))
    else:
        D, d = (0, 0)
    for i in range(len(x)-1):
        d = 2*trapz(y[i:i+2], x[i:i+2])
        D += d
        if D > p:
            if (D-p) < err:
                return x[i], x[i+1]
            else:
                raise LowAccuracyError('Need more (refine) data points!\n\
Upper error is {} while required accuracy is {}.'.format(d, err))
        else:
            pass
    return -1                     # won't happen usually


def get_dkp(p, *, depth, err, pr_delta):
    calc(pr_delta, 'relax')
    for _ in range(depth):
        try:
            find_width(p, err=err)
        except ThinIntervalError:
            calc(pr_delta, 'extend')
        except LowAccuracyError:
            calc(pr_delta, 'refine')
        else:
            break

    try:
        return find_width(p,err=err)[0]
    except LowAccuracyError:
        logger = logging.getLogger('dmslog')
        logger.warning('Failed to find DOB width within\
 {} search!'.format(depth))
        intv = find_width(p,err=np.inf)
        logger.warning("It's assumed to be in ({:.3g},{:.3g}) in the last\
 search.".format(intv[0], intv[1]))
        return intv[0]



def plot_pr(*, show=False):
    x = np.load(_args['fn_x'])
    y = np.load(_args['fn_y'])
    plt.scatter(x, y, label='SET: {}\nQ: {}'.format(
                _args['SET'], _args['Q']))
    plt.legend()
    if show:
        plt.show()
    else:
        plt.savefig(_args['fn']+'.png')
    return 'Totoal Integration:', 2*trapz(y, x)


def _main():
    setting_args(Q=.33, SET='Ca')
    # calc(pr_delta, 'refine')
    print(find_width(.95))
    print(plot_pr())
    setting_args(Q=.33, SET='A')
    # calc(pr_delta, 'refine')
    print(find_width(.95))
    print(plot_pr())
    plt.show()



if __name__ == '__main__':
    import cProfile
    cProfile.run('_main()', filename='plpr.out', sort='cumulative')
