# -*- coding: utf-8 -*-
# write in Cython
# prior p.d.f. for Set A, B, C

from libc.stdlib cimport malloc, free

cdef extern from 'math.h':
    double exp(double x)
    double log(double x)
    double cos(double x)
    double sin(double x)
    # double abs(double x)
    double sqrt(double x)

# setting hyper parameters
# cdef int h = 4
# cdef int k = 2
# cdef int n_c = 3
cdef double cbar_le = .001
cdef double cbar_ge = 1000
cdef double sigma = 1.

# cdef double Q = .2
# cdef double * ccck = [1., 1., 1.]

# pre-calulated constant
cdef double pi = 3.14159265358979323846
cdef double _one_over_sqrt_two_pi = .39894228        # 1/sqrt(2*pi)
cdef double _one_over_log_ge_over_le = .07238241365  # 1/log(1000*1000)

cdef double _theta(double x):
    if x >= 0:
        return 1
    else:
        return 0

def _pr_cn_if_cbar_A(double cn, double cbar):
    return _theta(cbar-abs(cn))/cbar*.5

def _pr_cbar_A(double cbar):
    return 1/cbar*_one_over_log_ge_over_le * _theta(
            cbar-cbar_le)* _theta(cbar_ge-cbar)

def _pr_cn_if_cbar_B(double cn, double cbar):
    return _theta(cbar-abs(cn))/cbar*.5

def _pr_cbar_B(double cbar):
    return exp(-log(cbar)*log(cbar)/(2*sigma*sigma)
    ) /(cbar*sigma) * _one_over_sqrt_two_pi

def _pr_cn_if_cbar_C(double cn, double cbar):
    return exp(-.5 * cn*cn/(cbar*cbar))/cbar * _one_over_sqrt_two_pi

def _pr_cbar_C(double cbar):
    return 1/cbar*_one_over_log_ge_over_le * _theta(
            cbar-cbar_le)* _theta(cbar_ge-cbar)

def _A_delta_if_cbar_f(double t, double delta,
                       double cbar, double Q, int h, int k):
    cdef double y = cos(t*delta)
    cdef double x = Q**(k+1)*t*cbar
    cdef int n = 0
    if x == 0:
        return 1
    else:
        while n < h:
            y *= sin(x)/x
            x *= Q
            n += 1
        return y

def _weight(double cbar, int n_c):
    return exp(-.5*log(cbar)*log(cbar))/cbar**(n_c+1)
