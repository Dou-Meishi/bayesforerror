#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pstats
import sys

p = len(sys.argv)-1

if p < 1:
    print('Please give file name!')
    sys.exit()
else:
    fn = sys.argv[1]
    sort_method = ['tottime']

if p >= 2:
    sort_method = sys.argv[2:]
else:
    pass

p = pstats.Stats(fn)
p.sort_stats(*sort_method).print_stats(10)
