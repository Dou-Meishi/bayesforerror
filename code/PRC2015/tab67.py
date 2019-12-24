# -*- coding: utf-8 -*-
# Calculate PRC2015 TABLE. VI, VII
# using data of PRC2015 TABLE. IV, V respectively

import numpy as np
from numpy import sqrt

# first column is p_rel, not T_lab
dataIV = [[153, 183.6, 166.5, 167.0, 166.8, 167.5],
          [212, 84.8, 75.1, 78.3, 77.5, 78.0],
          [259, 52.5, 49.1, 54.2, 53.7, 53.9],
          [307, 34.9, 35.9, 42.6, 43.2, 42.7]]
dataV =  [[153, 159.4, 164.8, 165.6, 167.2, 167.9],
          [212, 60.2, 68.9, 71.3, 78.1, 78.5],
          [259, 30.8, 38.6, 41.4, 52.6, 52.7],
          [307, 17.2, 22.5, 25.0, 38.6, 38.3]]

def _get_ccck(p, Lambda, sigmas):
    '''Extract ccck from experimental data (See Eq. (37)).'''
    Q = p/Lambda
    sigmas = [sigma/sigmas[0] for sigma in sigmas]
    ccck = [(sigmas[n]-sigmas[n-1])/Q**(n+1) for n in
            range(1, len(sigmas))]
    ccck.insert(0, 1.0)         # choose sigma_{LO} as refs
    ccck.insert(1, 0.0)         # c_1 is absent
    return ccck


def main():
    fmt = 'p: {:d}\nccck: {:.3g}, {:.3g}, {:.3g}, {:.3g}, \
{:.3g}, {:.3g}'

    print('\nTABLE. VI\n')
    for data in dataIV:
        p, Lambda = (data[0], 600) # R=0.9 fm
        print(fmt.format(p, *_get_ccck(p, Lambda, data[1:])))

    print('\nTABLE. V\n')
    for data in dataV:
        p, Lambda = (data[0], 400) # R=1.2 fm
        print(fmt.format(p, *_get_ccck(p, Lambda, data[1:])))
    return



if __name__ == '__main__':
    main()
