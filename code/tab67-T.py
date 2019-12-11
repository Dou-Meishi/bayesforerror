# -*- coding: utf-8 -*-
# Calculate PRC2015 TABLE. VI, VII
# using data of PRC2015 TABLE. IV, V respectively

import numpy as np
from numpy import sqrt

# first column is T_lab, not p_rel
dataIV = [[50, 183.6, 166.5, 167.0, 166.8, 167.5],
          [96, 84.8, 75.1, 78.3, 77.5, 78.0],
          [143, 52.5, 49.1, 54.2, 53.7, 53.9],
          [200, 34.9, 35.9, 42.6, 43.2, 42.7]]
dataV =  [[50, 159.4, 164.8, 165.6, 167.2, 167.9],
          [96, 60.2, 68.9, 71.3, 78.1, 78.5],
          [143, 30.8, 38.6, 41.4, 52.6, 52.7],
          [200, 17.2, 22.5, 25.0, 38.6, 38.3]]

def _get_ccck(T, Lambda, sigmas):
    '''Extract ccck from experimental data (See Eq. (37)).'''
    Q = sqrt(.5*939*T)/Lambda   # 939 MeV is the mass of nucleon
    sigmas = [sigma/sigmas[0] for sigma in sigmas]
    ccck = [(sigmas[n]-sigmas[n-1])/Q**(n+1) for n in
            range(1, len(sigmas))]
    ccck.insert(0, 1.0)         # choose sigma_{LO} as refs
    ccck.insert(1, 0.0)         # c_1 is absent
    return ccck


def main():
    fmt = 'T: {:d}\nccck: {:.3g}, {:.3g}, {:.3g}, {:.3g}, \
{:.3g}, {:.3g}'

    print('\nTABLE. VI\n')
    for data in dataIV:
        T, Lambda = (data[0], 600) # R=0.9 fm
        print(fmt.format(T, *_get_ccck(T, Lambda, data[1:])))

    print('\nTABLE. V\n')
    for data in dataV:
        T, Lambda = (data[0], 400) # R=1.2 fm
        print(fmt.format(T, *_get_ccck(T, Lambda, data[1:])))
    return



if __name__ == '__main__':
    main()
