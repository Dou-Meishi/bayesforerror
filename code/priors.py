# -*- coding: utf-8 -*-
# prior p.d.f. for Set A, B, C

from numpy import exp, log, sqrt

# setting hyper parameters
cbar_le, cbar_ge = (.001, 1000)
sigma = 1.

# pre-calulated constant
_one_over_sqrt_two_pi = .39894228        # 1/sqrt(2*pi)
_one_over_log_ge_over_le = .07238241365  # 1/log(1000*1000)

def _theta(x):
    if x >= 0:
        return 1
    else:
        return 0

def _pr_cn_if_cbar_A(cn, cbar):
    return _theta(cbar-abs(cn))/cbar*.5

def _pr_cbar_A(cbar):
    return 1/cbar*_one_over_log_ge_over_le * _theta(
            cbar-cbar_le)* _theta(cbar_ge-cbar)

def _pr_cn_if_cbar_B(cn, cbar):
    return _theta(cbar-abs(cn))/cbar*.5

def _pr_cbar_B(cbar):
    return exp(-log(cbar)*log(cbar)/(2*sigma*sigma)
    ) /(cbar*sigma) * _one_over_sqrt_two_pi

def _pr_cn_if_cbar_C(cn, cbar):
    return exp(-1/2 * cn*cn/(cbar*cbar))/cbar * _one_over_sqrt_two_pi

def _pr_cbar_C(cbar):
    return 1/cbar*_one_over_log_ge_over_le * _theta(
            cbar-cbar_le)* _theta(cbar_ge-cbar)
