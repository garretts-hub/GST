# -*- coding: utf-8 -*-
"""
Created on Wed Aug 15 15:25:23 2018

@author: GA28573
"""

import numpy as np
from scipy.fftpack import dct, idct
import pylab as plt


class Drift(object):
    def __init__(self, timestamps_array, data_array, nCounts=None):
        '''Input an array of of the num of 1 counts per timestamp'''
        self.data = data_array
        self.times = timestamps_array
        self.samples = len(data_array)
        self.counts = nCounts
        
        summed_timestep = 0
        for i in range(1, self.samples):
            summed_timestep += self.times[i] - self.times[i-1]
        self.avg_timestep = summed_timestep/self.samples
        
        #check if input is valid
        for point in self.data:
            if nCounts != None:
                if point > self.counts:
                    raise ValueError("Cannot have a sample with greater than nCounts!")
                elif point < 0:
                    raise ValueError("Cannot have a sample with negative counts!")
                    
        if len(self.data) != len(self.times):
            raise ValueError("Must have the same number of data points and timestamps!")
        
        #initialize variables
        self.timesteps = [0]*(len(self.times)-1) #will list the time spacing between all 
        for i in range(0, len(self.times)-1 ):
            self.timesteps[i] = self.times[i+1] - self.times[i]
        self.avg_timestep = np.mean(self.timesteps)
        self.frequencies = np.zeros(self.samples)
        self.modes = np.zeros(self.samples)
        self.powers = np.zeros(self.samples)
        
    def __repr__(self):
        return "Drift(Len={}, start_time={}s, end_time={}s, counts_per_sample={}, average_timestep={:.3}s".format(self.samples, self.times[0], self.times[-1], self.counts, self.avg_timestep)
    
    def _dct(self, null_hypothesis=None):
        """    
        x : array; Data string, on which the normalization and discrete cosine transformation is performed. If
            counts is not specified, this must be a bit string.
        null_hypothesis : array, optional
            If not None, an array to use in the normalization before the DCT. If None, it is
            taken to be an array in which every element is the mean of x.
        """
        x = self.data
        nCounts = self.counts
        x_mean = np.mean(x)
        N = len(x) 
        # If the null hypothesis is not specified, we take our null hypothesis to be a constant bias
        # coin, with the bias given by the mean of the data / number of counts.
        if null_hypothesis is None:    
            null_hypothesis = x_mean/nCounts
            if null_hypothesis <= 0 or null_hypothesis >= 1:
                return np.zeros(N)
        input_array= (x - nCounts*null_hypothesis)/np.sqrt(nCounts*null_hypothesis * (1 - null_hypothesis))
        return dct(input_array, norm='ortho')
    
    def _manual_ndft(self):
        times = self.times
        vals = self.data
        N = len(times)
        T = times[-1]
        print("Will return {}-frequencies, spaced at {} Hz".format(N, (1/T)))
        self.frequencies = np.arange(N)*(1/T)
        for m in range(len(self.frequencies)):
            summed = 0
            for n in range(N):
                summed += vals[n]*np.exp(-2j*np.pi*m/T*times[n])
            self.modes[m] = summed.real
        self.powers = self.modes**2
        return self.frequencies, self.powers
    
    def plot_input(self):
        plt.plot(self.times, self.data,marker='.')
        plt.grid()
        plt.xlabel("Time, s")
        plt.ylabel("Measured Value, a.u.")
        plt.title("Data Set to be Transformed")
        plt.show()
    
        
    def plot_power_spectrum(self):
        plt.plot(self.frequencies, self.powers, marker='.')
        plt.grid()
        plt.xlabel("Frequency, Hz")
        plt.ylabel("Power, a.u.")
        plt.title("Power Spectrum")
        plt.show()
    
    
    