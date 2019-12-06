# -*- coding: utf-8 -*-
# Filename: tab3-2.py
# print result of TABLE. III

import numpy as np
import plpr, prde
import logging.config, yaml # logging module

def tab3C(p, **args):
    prde.setting_args(**args)
    plpr.setting_args(**args)
    plpr.calc(prde.pr_delta_C, 'relax')
    depth = 10

    for _ in range(depth):
        try:
            plpr.find_width(p, err=.001)
        except plpr.ThinIntervalError:
            plpr.calc(prde.pr_delta_C, 'extend')
        except plpr.LowAccuracyError:
            plpr.calc(prde.pr_delta_C, 'refine')
        else:
            break

    try:
        return plpr.find_width(p,err=.001)[0]
    except plpr.LowAccuracyError:
        logger = logging.getLogger('dmslog')
        logger.warning('Failed to find DOB width within\
 {} search!'.format(depth))
        intv = plpr.find_width(p,err=np.inf)
        logger.warning("It's assumed to be in ({:.3g},{:.3g}) in the last\
 search.".format(intv[0], intv[1]))
        return intv[0]


def tab3A(p, **args):
    plpr.setting_args(**args)
    prde.setting_args(**args)
    plpr.calc(prde.pr_delta_A, 'relax')
    depth = 0

    for _ in range(depth):
        logger = logging.getLogger('dmslog').getChild(__name__)
        try:
            plpr.find_width(p, err=.001)[0]
        except plpr.ThinIntervalError:
            logger.info('Extending interval...')
            plpr.calc(prde.pr_delta_A, 'extend')
        except plpr.LowAccuracyError:
            logger.info('Refining interval...')
            plpr.calc(prde.pr_delta_A, 'refine')
        else:
            break

    try:
        return plpr.find_width(p,err=.001)[0]
    except plpr.LowAccuracyError:
        logger = logging.getLogger('dmslog')
        logger.warning('Failed to find DOB width within\
 {} search!'.format(depth))
        intv = plpr.find_width(p,err=np.inf)
        logger.warning("It's assumed to be in ({:.3g},{:.3g}) in the last\
 search.".format(intv[0], intv[1]))
        return intv[0]


def ini_log():
    with open('./logconf.yaml', 'r') as f:
        logconf = yaml.safe_load(f.read())
    logging.config.dictConfig(logconf)
    return

def main():
    ini_log()
    
#     for SET in ['A']:
#         p = .68
#         for Q in [.20, .33, .50]:
#             print('{:.0%} DOB width of Set {} with Q {} is \n\t\
# {:.3g}'.format(p, SET, Q, tab3A(p, SET=SET, Q=Q)))
#         p = .95
#         for Q in [.20, .33, .50]:
#             print('{:.0%} DOB width of Set {} with Q {} is \n\t\
# {:.3g}'.format(p, SET, Q, tab3A(p, SET=SET, Q=Q)))

    for p in [.68, .95]:
        for Q in [.20, .33, .50]:
            for SET in ['Ca', 'Cb', 'Cc']:
                print('{:.3g}'.format(tab3C(p, SET=SET, Q=Q)))

    logger = logging.getLogger('dmslog').getChild(__name__)
    logger.info('Main program finished.')

if __name__ == '__main__':
    import cProfile
    cProfile.run('main()', filename='tab3-2.out')
