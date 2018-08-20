# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 20:24:35 2018

@author: Garrett
"""

import numpy as np
import pylab as plt
import nfft

def make_time(start, stop, num, sigma):
    t = np.linspace(start, stop, num)
    for i in range(len(t)):
        t[i] += np.random.normal(0,sigma)
        while t[i] < start or t[i] > stop:
            t[i] += np.random.normal(0,sigma)
    return t

N = 50
freq = 4.5
start_t = 0
stop_t = 3
num_points = 30
sigma = 3
x = make_time(start_t, stop_t, N, sigma)
x0 = np.linspace(start_t, stop_t,800)
f = np.sin(2*np.pi*freq*x)
f0 = np.sin(2*np.pi*freq*x0)

plt.plot(x,f, label="Irregular",marker='.',markersize=8,ls="None")
plt.plot(x0,f0, label="Regular")
plt.grid()
plt.legend(loc='lower right')
plt.xlabel("Time, seconds")
plt.title("Data Points")
plt.show()

k = np.arange(-N//2, N//2)
powers = np.absolute(nfft.nfft_adjoint(x,f,N))
plt.plot(k,powers.real,label="Real",marker='.')
#plt.plot(k, powers.imag, label="Imaginary")
plt.title("Non-Equispaced Fast Fourier Transform")
plt.legend()
plt.grid()
plt.xlim(0, k[-1])
plt.show()