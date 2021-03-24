# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 22:08:57 2017
Author: Shih-Chao Wei
------------------------------------------------------------------------------
Poles-Zeros Calculation for Passive Velocity Transducer (Geophone) Responses
Ref: https://ds.iris.edu/NRL/sensors/sercel/passive_responses.html
------------------------------------------------------------------------------
Definition of Parameters: 
    fo = natural frequency (Hz)
    fn = normalization frequency (Hz)     
    m = mass (Kg)
    Go = instrinsic sensitivity = sqrt(4*pi*fo*m*Rt*bc)
    Ge = effective sensitivity = Go*(Rload/Rt)

    Rc = feedback coil resistance
    Rs = shunt resistor (if there is no Rs, enter -9999)
    Zamp = input impedance of the preamp or datalogger     
    Rload = parallel sum of Rs and Zamp = (Rs*Zamp) / (Rs+Zamp)
                                        = Zamp, if there is no Rs
    Rt = total resistance = Rc + Rload
    bo = open circuit damping due to the mechanical properties of the sensor
    bc = current damping due to the electrical properties of the sensor
    bt = total damping = bo + bc 
    
    Res = Resolution of datalogger (bits)
    MSR = Max. Sampling Rate of datalogger (Sampling/sec.)
    
    A0 = normalization factor = gain
    k1 = Sensor sensitivity: Output ground motion to voltage (Volt/M/S) (or Volt/M/S/S)
    k2 = Datalogger conversion factor: Convert input voltage into digital count (Count/Volt)
"""
import numpy as np
import matplotlib.pyplot as plt

## GS-32CT data (395 Ohm)
#fo = 10 # (Hz)
#fn = 200 # (Hz)
#m = 0.0112 # (kg) = 11.2 (g) = 0.395 (oz)
#Go = 27.5 # (V/m/s) = 0.275 (V/cm/s) = 0.698 (V/in/s)
#Rc = 395 # (Ohm)
#Rs = -9999
#bo = 0.316 # = 31.6 (%)

## GS-32CT data (395 Ohm) In
#fo = 10 # (Hz)
#fn = 200 # (Hz)
#m = 0.395 # (oz)
#Go = 0.698 # (V/in/s)
#Rc = 395 # (Ohm)
#Rs = -9999
#bo = 0.316 # = 31.6 (%)

## GS-20DX data (395 Ohm)
#fo = 10 # (Hz)
#fn = 200 # (Hz)
#m = 0.011 # (kg) = 11.0 (g) 
#Go = 28 # (V/m/s) = 0.28 (V/cm/s) = 0.70 (V/in/s)
#Rc = 395 # (Ohm)
#Rs = -9999
#bo = 0.3 # = 30 (%)

# GS-20DX data (395 Ohm) In
fo = 10 # (Hz)
fn = 200 # (Hz)
m = 0.388 # (oz) 
Go = 0.7 # 0.70 (V/in/s)
Rc = 395 # (Ohm)
Rs = -9999
bo = 0.3 # = 30 (%)

## USB-4716
#Res = 16 # (bits)
#MSR = 200*1000 # (S/s) = 200 (kS/s) max.
#Zamp = 100*10**6 + 1/(1j*MSR*100*10**(-12)) # 100 MΩ/100 pF

# PCI-1713U
Res = 12 # (bits)
MSR = 100*1000 # (S/s) = 100 (kS/s) max.
Zamp = 10**9 # 1 GΩ

if Rs == -9999:
    Rload = Zamp
else:
    Rload = (Rs*Zamp) / (Rs+Zamp) # parallel sum of Rs and Zamp
Rt = Rc + Rload # total resistance
bc = (Go**2)/(4*np.pi*fo*m*Rt) 
bt = bo + bc # total damping

p1 = (-2*np.pi*fo*bt) + 1j*(2*np.pi*fo*np.sqrt(1-bt**2))
p2 = (-2*np.pi*fo*bt) - 1j*(2*np.pi*fo*np.sqrt(1-bt**2))
z1 = 0.0 + 0j
z2 = 0.0 - 0j

fl_p, fu_p, nof = (-3, 3, 10000) # range: 10^-3~10^3, sampling: 1000
df_p = (fu_p-fl_p) / nof # interval
fre = np.linspace(0, nof, nof+1)
for i in range(0, nof+1, 1):
    fre[i] = 10**(fl_p + df_p * i)
wn = 2*np.pi*fn
Ao = 1/np.abs( (1j*wn-z1)*(1j*wn-z2)/(1j*wn-p1)/(1j*wn-p2) )
k1 = Go
k2 = (2**(Res)-1)/5 # full scale = +/- 2.5V for 16 bit datalogger
const = Ao*k1
w = 2*np.pi*fre
Amp = np.abs( const * (1j*w-z1)*(1j*w-z2)/(1j*w-p1)/(1j*w-p2) )
Pha = np.angle( const * (1j*w-z1)*(1j*w-z2)/(1j*w-p1)/(1j*w-p2), deg=True)




Fig, ax = plt.subplots(2, 1, figsize=(6,8))
ax[0].plot(fre, Amp, 'r-')
ax[0].grid(True)
ax[0].set_ylim(0.05, 2.0)
#ax[0].set_ylim(5, 100)
ax[0].set_xlim(5, 500)
#ax[0].set_ylabel("$Amplitude (V/m/s)$", fontsize=14)
ax[0].set_ylabel("$Amplitude (V/in/s)$", fontsize=14)
ax[1].semilogx(fre, Pha, 'k-')
ax[1].grid(True)
ax[1].set_ylim(0, 180)
ax[1].set_xlim(5, 500)
ax[1].set_ylabel("$Phase (degree)$", fontsize=14)
ax[1].set_xlabel("$Frequency (Hz)$", fontsize=14)
plt.show()
#plt.savefig('Instrument Response Figure_GS-32CT_In.pdf')
#plt.savefig('Instrument Response Figure_GS-20DX_In.pdf')