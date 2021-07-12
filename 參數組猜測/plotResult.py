# -*- coding: utf-8 -*-
"""
Created on Sun Jul 11 23:54:58 2021

@author: r8891
"""

import numpy as np
import matplotlib.pyplot as plt

G0=27.6
a0=62.8
h00=0.3
b=22000008690

G=[26.772,36.71,16.339]
a=[55.264,59.032,64.684]
h0=[0.3,0.891,0.0004]

fd=np.linspace(0,250,2000)
wd=fd*2*np.pi

for i in range(3):
    T=G0*wd**2/(a0**4-2*a0**2*wd**2+wd**4+4*(h00+G0**2/(a0*b))**2*a0**2*wd**2)**0.5
    Ts=G[i]*wd**2/(a[i]**4-2*a[i]**2*wd**2+wd**4+4*(h0[i]+G[i]**2/(a[i]*b))**2*a[i]**2*wd**2)**0.5
    
    plt.figure(i)
    plt.title("FreqResponse_Function")
    plt.loglog(basex=10,basey=10)
    plt.xlim(1,1000)
    plt.ylim(1,100)
    plt.xlabel("Freq(Hz)")
    plt.ylabel("Amp_T(Volt/m/s)")
    
    plt.plot(fd,T)
    plt.plot(fd,Ts)
    
    plt.grid(which='major', axis='both', linewidth=0.75)
    plt.grid(which='minor', axis='both', linewidth=0.75)
                                                                     
                                                            