#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 10:32:14 2020

@author: bioc1642
"""

import numpy as np
import matplotlib.pyplot as plt

aa = 10
bb = 0.01
cc = np.arange(2.0,10.0,1)
a = aa
b = bb
c = cc
d = 10

for c in cc:
    x = np.linspace(1,400,10000)
    y = a * c**(-b*x) + d
    plt.plot(x,y,label='c=%1.2f' % c)

plt.legend()
plt.show()
