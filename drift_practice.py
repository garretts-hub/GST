# -*- coding: utf-8 -*-
"""
Created on Mon Aug 20 11:41:29 2018

@author: GA28573
"""

import Drift
from Drift import Drift
import numpy as np
import pylab as plt


def make_time(start, stop, num, sigma):
    t = np.linspace(start, stop, num)
    for i in range(len(t)):
        t[i] += np.random.normal(0,sigma)
        while t[i] < start or t[i] > stop:
            t[i] += np.random.normal(0,sigma)
    return t

freq = 4
start_t = 0
stop_t = 5
num_points = 1001
sigma = (stop_t/num_points)*0.0
x = make_time(start_t, stop_t, num_points, sigma)
f = np.sin(2*np.pi*freq*x)
x0 = np.linspace(start_t, stop_t, 500)
f0 = np.sin(2*np.pi*freq*x0)

#drifted = Drift(x0, f0)
drifted = Drift(x, f)

drifted.plot_input()
drifted._manual_ndft(print_details=True)
drifted.plot_power_spectrum()