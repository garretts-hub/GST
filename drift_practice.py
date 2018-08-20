# -*- coding: utf-8 -*-
"""
Created on Mon Aug 20 11:41:29 2018

@author: GA28573
"""

import Drift
from Drift import Drift
import numpy as np
import pylab as plt

times = np.linspace(0, 3, 400)
freq = 6
vals = np.sin(2*np.pi*freq*times)

def make_time(start, stop, num, sigma):
    t = np.linspace(start, stop, num)
    for i in range(len(t)):
        t[i] += np.random.normal(0,sigma)
        while t[i] < start or t[i] > stop:
            t[i] += np.random.normal(0,sigma)
    return t

freq = 3
start_t = 0
stop_t = 3
num_points = 40
sigma = (stop_t/num_points)*0.00
x = make_time(start_t, stop_t, num_points, sigma)
f = np.sin(2*np.pi*freq*x)
x0 = np.linspace(start_t, stop_t, 500)
f0 = np.sin(2*np.pi*freq*x0)

#drifted = Drift(x0, f0)
drifted = Drift(x, f)

drifted.plot_input()
drifted._manual_ndft()
print("average timestep {:.4f} s".format(drifted.avg_timestep))
print("average sample rate is {:.3f} Hz".format(1/drifted.avg_timestep))
print("Actual max frequency is {} Hz".format(drifted.frequencies[-1]))
print("Effective max frequency is {} Hz".format(drifted.frequencies[drifted.samples//2]))
drifted.plot_power_spectrum()