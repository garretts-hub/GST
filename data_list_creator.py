# -*- coding: utf-8 -*-
"""
Created on Fri Aug  3 14:07:03 2018

@author: GA28573
"""

import sys
#sys.path.append("../../transmon-sim/physical_sim")
sys.path.append("C:/Users/GA28573/AppData/Local/Continuum/anaconda32/Lib/site-packages/qutip-4.3.1-py3.6-win-amd64.egg/qutip")


import NoiseSignal2 as _ns
import qutip as _qt
import numpy as np
import pylab as plt



def gate_string_to_list(gate_string):
    #returns gate sequence as a list
    #ex, 'GxGxGx' --> ['Gx', 'Gx', 'Gx', 'Gx']
    #"Gx(Gy)^2Gz" --> ['Gx', 'Gy', 'Gy', 'Gz']
    #do not raise anything to the 1 exponent, only 2 or more; always use parenthesis
    #currently can only handle exponents up to 99
    #do not put more than one gate in parentheses for right now
    gate_list = []
    for i in range(len(gate_string)):
        char = gate_string[i]
        if char == "G":
            gate_list.append('G'+ gate_string[i+1])
        elif char.isnumeric():
            if i != (len(gate_string)-1): 
                next_char = gate_string[i+1]
            else:
                next_char = ''
            prev_char = gate_string[i-1]
            if prev_char.isnumeric():
                pass
            if next_char.isnumeric():
                char = int(gate_string[i] + gate_string[i+1])
            for count in range(int(char)-1):
                gate_list.append('G' + gate_string[i-3])
    return gate_list

def gate_list_to_string(gate_list):
    gate_string = ''
    for gate in gate_list:
        gate_string = gate_string + gate
    return gate_string

'''a Garrett-made function. I made this just to export and make adding 1/f noise to my time-dependent data simulator a little easier'''
def telegraph_noise(total_time, amplitude, total_fluctuators, start_f, stop_f, time_unit):
    '''if initialization_time != 0:
        iterations = total_time//initialization_time'''
    signal = _ns.NoiseSignalTelegraph(initial_seed=1, time_unit=time_unit)
    signal.configure_noise(1.0, amplitude, total_fluctuators, start_f, stop_f, total_time)
    
    sample_val = []
    for ii in range(1):
        signal.init()
        signal.interpolation_settings(True)
    
        tvals = [t for t in range(total_time)]
        for t in tvals:
            sample_val += [signal[t]]
        tt = tvals

    plt.plot(tt, sample_val)
    plt.grid()
    plt.title("Telegraph Noise from {} to {} Hz".format(start_f, stop_f))
    plt.xlabel("Time in time_unit integers".format(time_unit))
    plt.show()
    return tt, sample_val

'''This simulates data based on an overrotation based on 1/f noise'''
def create_1f_data(time_per_count, nSamples, nCounts, gate_list, amp, total_fluctuators, start_f, stop_f, time_unit, time_between_samples=0, recalibration_time=None):
    '''To-do: make a way to recalibrate error after nCounts (i.e. each sample) even if the timebetweensamples is zero.'''
    time_per_sample = time_per_count*nCounts
    timestep = time_per_sample + time_between_samples
    timestamps = np.arange(timestep, nSamples*timestep, timestep)
    rho0 = _qt.operator_to_vector(_qt.ket2dm(_qt.basis(2,0)))
    rho1 = _qt.operator_to_vector(_qt.ket2dm(_qt.basis(2,1)))
    zero_counts = []
    one_counts = []
    probs = []
    probs_ideal = []
    angles = []
    
    biggest_stamp = int(timestamps[len(timestamps)-1]/time_unit)
    max_time = biggest_stamp + 1
    
    noise_time, noise_vals = telegraph_noise(max_time, amp, total_fluctuators, start_f, stop_f, time_unit)
    print("Greatest timestamp is {}\nMax time index in the noise list is {} (both in time_units)".format(biggest_stamp, max_time))
    
    for time in timestamps:
        time_in_timeunits = int(time/time_unit) #finds the time as an integer number of timeunits to use as an index for the noise list
        noise_at_time = noise_vals[time_in_timeunits]
        #print("Starting timestamp {:.3f} s".format(time))
        rho = rho0
        rho_ideal = rho0
        for gate in gate_list:
            if gate == 'Gx':
                angle = np.pi/2 + noise_at_time
                angles.append(angle)
                rho = _qt.to_super(_qt.rx(angle)) * rho
                rho_ideal = _qt.to_super(_qt.rx(np.pi/2)) * rho_ideal
            elif gate == 'Gi':
                angle = noise_at_time
                angles.append(angle)
                rho = _qt.to_super(_qt.rx(angle)) * rho
                rho_ideal = _qt.to_super(_qt.rx(0)) * rho_ideal
            elif gate == 'Gy':
                np.pi/2 + noise_at_time
                angles.append(angle)
                rho = _qt.to_super(_qt.ry(angle)) * rho
                rho_ideal = _qt.to_super(_qt.ry(np.pi/2)) * rho_ideal
            elif gate == 'Gz':
                np.pi/2 + noise_at_time
                angles.append(angle)
                rho = _qt.to_super(_qt.rz(angle)) * rho
                rho_ideal = _qt.to_super(_qt.rz(np.pi/2)) * rho_ideal
            
        #print(rho)
        #calculate probabilities of being in 1 after the experiment has been applied
        p1 = (rho.dag()*rho1).norm()
        p1_ideal = (rho_ideal.dag()*rho1).norm()
        if p1 >= 1:
            p1 = 1
        elif p1<0:
            p1 = 0
        #print("*****Time {:.3f}, prob {:.3f}".format(time, p1))
        probs.append(p1)
        probs_ideal.append(p1_ideal)
        one_count = np.random.binomial(nCounts, p1) #simulates a summation of the number of 1-counts you get in one bitstring sample
        zero_count = nCounts - one_count #simulates summation of the number of 0-counts in one bitstring sample
        one_counts.append(one_count)
        zero_counts.append(zero_count)
   
    plt.plot(timestamps, probs_ideal, label="Ideal")
    plt.plot(timestamps, probs, ls='none', marker='.', color="black", label="Noisy")
    plt.ylim(0,1)
    plt.xlabel("Time, seconds")
    plt.ylabel("Probability of Measuring State {1}")
    plt.title("Simulated {} with 1/f Overrotation Error".format(gate_list_to_string(gate_list)))
    plt.grid()
    plt.legend()
    plt.show()
        
    return (np.asarray(one_counts), np.asarray(zero_counts), np.asarray(timestamps), probs)



def create_data(time_per_count, num_samples, num_counts, gate_list, time_unit, noise_type=None, walking_amp=None, telegraph_amp=None, \
                res=None, freq_list=None, amp_list=None, phase_list=None, start_f=None, stop_f=None, fluctuators=None, plot_noise=False, add_noise=False):
    #time_per_shot: time in seconds for a single (prep-gate-measure+delay)
    #num_samples: how many timestamps and strings of counts you want to have
    #num_counts: how many data points to create (how many ones and zeros) per sample (i.e. per timestamp) --> affects both time of a sample and precision
    #num_shots: how many times (shots) you apply (prep-gate_list-meas) to get one count (a single 0 or 1) --> determines the time of one count, but won't affect precision
    #gate_list: gates you want to do for your operation, entered as a list of strings
    #xerr,yerr,zerr: 2D tuples with overrotation amplitude in radians and frequency in Hz
    rho0 = _qt.operator_to_vector(_qt.ket2dm(_qt.basis(2,0)))
    rho1 = _qt.operator_to_vector(_qt.ket2dm(_qt.basis(2,1)))
    zero_counts = []
    one_counts = []
    timestep = num_counts*time_per_count #the time to get a full bitstring of zeros and ones for one sample, i.e. one timestamp
    timestamps = np.arange(timestep, num_samples*timestep, timestep) #array of 1*timestep, 2*timestep,....(num_samples)*timestep
    probs = []
    angles = []
    total_time = (time_per_count*num_counts*num_samples)/time_unit
    
    sig = 0
    if noise_type == "Sine":
        sig = _ns.NoiseSignalSine(time_unit=time_unit)
        sig.configure_noise(resolution_factor=res, freq_list=freq_list, amp_list=amp_list, phase_list=phase_list, total_time=total_time)
        sig.init()
        if add_noise != None:  
            sig.add_random_noise(add_noise) #add normal noise with specified std deviation if requested
    elif noise_type == "Random Walk":
        sig = _ns.NoiseSignalRandomWalk(initial_seed=1234, time_unit=time_unit)
        sig.configure_noise(walking_amp, res, total_time)
        sig.init() 
    elif noise_type == "Telegraph":
        sig = _ns.NoiseSignalTelegraph(initial_seed=1234, time_unit=time_unit)
        sig.configure_noise(exponent=1, amplitude=telegraph_amp, total_fluctuators=fluctuators, start_freq=start_f, stop_freq=stop_f, total_time=total_time)
        sig.init()
        sig.interpolation_settings(do_interpolation=True, resolution_factor=res)
    
    if plot_noise==True:
        sig.plot_noise_signal()
    
    for time in timestamps:
        noise_at_time = 0
        if noise_type != None:
            #print(time/time_unit)
            noise_at_time = sig[time/time_unit]
        rho = rho0
        
        for gate in gate_list:
            '''
            Next step: calculate change in rotation error within each shot. Currently takes the time at the start of the experiment shot
            and applies that to all gates in one shot. Depending on the timescale of the error and time per shot, this simplification may need
            to be addressed so that each gate, say each Gx in (Gx)^11, has an error associated with its specific time, not the same error for
            all 11 Gx gates.
            '''
            if gate == 'Gx':
                angle = np.pi/2 + noise_at_time
                angles.append(angle)
                rho = _qt.to_super(_qt.rx(angle)) * rho
            elif gate == 'Gy':
                angle = np.pi/2 + noise_at_time
                angles.append(angle)
                rho = _qt.to_super(_qt.ry(angle)) * rho
            elif gate == 'Gz':
                angle = np.pi/2 + noise_at_time
                angles.append(angle)
                rho = _qt.to_super(_qt.rz(angle)) * rho
        #calculate probabilities of being in 1 after the experiment has been applied
        p1 = (rho.dag()*rho1).norm()
        probs.append(p1)
        one_count = np.random.binomial(num_counts, p1) #simulates a summation of the number of 1-counts you get in one bitstring sample
        zero_count = num_counts - one_count #simulates summation of the number of 0-counts in one bitstring sample
        one_counts.append(one_count)
        zero_counts.append(zero_count)
    
    if plot_noise == True:
        plt.plot(timestamps, probs)
        plt.ylim(0,1)
        plt.xlabel("Time, seconds")
        plt.ylabel("Probability of Measuring State {1}")
        plt.title("Simulated {} with {} Noise".format(gate_list_to_string(gate_list), noise_type))
        plt.grid()
        plt.show()
        
    return (np.asarray(one_counts), np.asarray(zero_counts), np.asarray(timestamps), probs)



if __name__=='__main__':

    gate_list = gate_string_to_list("GxGxGxGxGx")
    #print(gate_string_to_list("GxGxGxGxGxGxGxGxGxGxGx"))
    #print(gate_string_to_list("(Gx)^5"))
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
    freq_list=(0.8,)
    amp_list=(0.04,)
    phase_list=(0,)
    add_noise=0.01
    start_f = 0.1
    stop_f = 2
    fluctuators= 40
    
    
    
    #create_1f_data(time_per_count, nSamples, nCounts, gate_list, amp, fluctuators, start_f, stop_f, time_units)
    #create_data(time_per_count, nSamples, nCounts, gate_list, time_units, noise_type, walking_amp, telegraph_amp, \
    #            res, freq_list, amp_list, phase_list, start_f, stop_f, fluctuators,plot_noise,add_noise)
    
    
    
    