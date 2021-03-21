# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 16:48:00 2021

@author: r8891
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import pi
from scipy.fftpack import fft


sample_num = 2000 # 取樣點數
total_time = 2 # 總時間
sampling_rate = sample_num / total_time # 取樣頻率

fs = [(2, 20), (8, 5), (16, 2)] # sin 波的頻率與振幅組合。 (Hz, Amp)

noise_mag = 0
time = np.linspace(0, total_time, sample_num, endpoint=False)
vib_data = [amp * np.sin(2*pi*hz*time) for hz, amp in fs]

max_time = int(sample_num / 4)

plt.figure(figsize=(12, 8))
# Show seperated signal
for idx, vib in enumerate(vib_data):
    plt.subplot(2, 2, idx+1)
    plt.plot(time[0:max_time], vib[0:max_time])
    plt.xlabel('time')
    plt.ylabel('vib_' + str(idx))
    plt.ylim((-24, 24))

vib = sum(vib_data) 
plt.subplot(2, 2, 4)
plt.plot(time[0:max_time], vib[0:max_time])
plt.xlabel('time')
plt.ylabel('vib(combine)')
plt.ylim((-24, 24))


#反轉過去，確定scipy.fftpack這個套件的轉換可不可行
fd = np.linspace(0.0, sampling_rate, int(sample_num), endpoint=False)
vib_fft = fft(vib)
mag = 2/sample_num * np.abs(vib_fft) # Magnitude

plt.figure()
plt.plot(fd[0:int(sample_num/2)], mag[0:int(sample_num/2)])
plt.xlabel('Hz')
plt.ylabel('Amp')