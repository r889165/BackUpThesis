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




G=[G0*0.96,G0*1.19,G0*1.05]
a=[a0*0.89,a0*0.93,a0*0.86]
h0=[h00*1.09,h00*1.23,h00*1.239]

G2=[G0*0.95,G0*1.169,G0*1.301]
a2=[a0*1.006,a0*0.874,a0*1.155]
h02=[h00*1.099,h00*1.166,h00*1.195]

fd=np.linspace(0,250,2000)
wd=fd*2*np.pi

for i in range(3):
    T=G0*wd**2/(a0**4-2*a0**2*wd**2+wd**4+4*(h00+G0**2/(a0*b))**2*a0**2*wd**2)**0.5
    #Ts=G[i]*wd**2/(a[i]**4-2*a[i]**2*wd**2+wd**4+4*(h0[i]+G[i]**2/(a[i]*b))**2*a[i]**2*wd**2)**0.5
    Ts2=G2[i]*wd**2/(a2[i]**4-2*a2[i]**2*wd**2+wd**4+4*(h02[i]+G2[i]**2/(a2[i]*b))**2*a2[i]**2*wd**2)**0.5
    
    
    
    plt.title("FreqResponse_Function")
    plt.loglog(basex=10,basey=10)
    plt.xlim(1,1000)
    plt.ylim(1,100)
    plt.xlabel("Freq(Hz)")
    plt.ylabel("Amp_T(Volt/m/s)")
    
    plt.plot(fd,T,color="black")
    #plt.plot(fd,Ts,color="black")
    #plt.plot(fd,Ts2)
    
    plt.grid(which='major', axis='both', linewidth=0.75)
    plt.grid(which='minor', axis='both', linewidth=0.75)
                                                                     
                                                            