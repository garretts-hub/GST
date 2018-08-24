# -*- coding: utf-8 -*-
"""
Created on Mon Aug 20 11:41:29 2018

@author: GA28573
"""

import Drift
import NoiseSignal2 as _ns
from data_list_creator import gate_string_to_list, gate_list_to_string, create_data
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

freq = 1
start_t = 0
stop_t = 30
num_points = 300
sigma = (stop_t/num_points)*0.2
x = make_time(start_t, stop_t, num_points, sigma)
f = np.sin(2*np.pi*freq*x)
x0 = np.linspace(start_t, stop_t, 500)
f0 = np.sin(2*np.pi*freq*x0)

gate_list = gate_string_to_list("GxGxGxGxGxGxGxGxGx")
nSamples = 1000  #total samples to take for each measurement
nCounts = 1      #total shots to take at one; =nSamples: same noise, probabilities for all repeats; =1, new experiment & noise for each count
time_per_count = 0.016 #seconds
time_units = 1e-3 #seconds
amp = 0.05
noise_type='Sine' #Sine, Random Walk, Telegraph
plot_noise=True
walking_amp = 0.001
telegraph_amp = 0.02
res = 1
freq_list=(4,)
amp_list=(0.1,)
phase_list=(0,)
add_noise=None
start_f = 0.1
stop_f = 2
fluctuators= 40


one_counts, zero_counts, timestamps, probs = \
    create_data(time_per_count, nSamples, nCounts, gate_list, time_units, noise_type, walking_amp, telegraph_amp, \
            res, freq_list, amp_list, phase_list, start_f, stop_f, fluctuators,plot_noise,add_noise)
    


#drifted = Drift(x0, f0)
#drifted = Drift(x, f)
drifted = Drift(timestamps, one_counts, nCounts = nCounts)

drifted.plot_input()
drifted._manual_ndft(print_details=True)
#drifted._dct(print_details = True)
drifted.plot_power_spectrum()