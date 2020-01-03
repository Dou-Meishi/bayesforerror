# -*- coding: utf-8 -*-
# use bayesian approach to estimate simple series

import numpy as np
import matplotlib.pyplot as plt
from numpy import cos, sin, exp, pi, log

import dmsPDF


def test1():
    '''1/(1-Q) = 1 + Q + Q^2 + Q^3 + ...'''
    Q = .5
    abcent_id = [0]

    for k in range(2,6):
        c0ck = [1 for i in range(k+1)]

        pdfB = dmsPDF.PrDelta('B', Q, c0ck, abcent_id)

        stddelta = 1/(1-Q)-sum([Q**m for m in range(k+1)])
        print("Standard delta: {}".format(stddelta))
        print(pdfB.get_dkp(.95, method='newton'))
        print('')

    return


def testcos():
    '''cos Q = 1 + (-1)^m/(2m)! Q^2m, m=1, 2, ...'''
    def _cos_cm():
        m, res = (0, 1)
        yield res
        while True:
            yield 0             # odd term is 0
            m += 1
            res *= -1./(4*m*m - 2*m)
            yield res

    cos_cm = _cos_cm()
    Q = .5
    c0ck = [next(cos_cm) for i in range(3)]
    abcent_id = [0, 1]

    pdfC = dmsPDF.PrDelta('C', Q, c0ck, abcent_id)
    stddelta = cos(Q) - sum([c0ck[m]*Q**m for m in range(len(c0ck))])

    for k in range(len(c0ck), 7):
        ck = next(cos_cm)
        if ck == 0:
            pdfC = pdfC.update_ck(ck, True)
        else:
            pdfC = pdfC.update_ck(ck, False)

        stddelta -= ck*Q**k
        print("  Standard delta: {:.3g}".format(stddelta))
        print("width of .68 DOB: {:.3g}".format(pdfC.get_dkp(.68, method='bisect')))
        print('')


    return


def main():
    # test1()
    testcos()

    return


if __name__ == '__main__':
    main()
    
