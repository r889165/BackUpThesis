import numpy as np
import matplotlib.pyplot as plt

total_time= 12
dt=0.001
sample_rate=1/dt
sample_num= int(sample_rate*total_time)
fd = np.linspace(0.0, int(sample_rate), int(sample_num), endpoint=False)


# GS-20DX
fo = 10 # (Hz)
fn = 200 # (Hz)
m = 11 # (g) 
Go = 27.6 # (V/m/s)
Rc = 395 # (Ohm)
Rs = -9999
bo = 0.3 # = 30 (%)


# PCI-1713U
Res = 12 # (bits)
MSR = 100*1000 # (S/s) = 100 (kS/s) max.
Zamp = 10**9 # 1 GΩ
if Rs == -9999:
    Rload = Zamp
else:
    Rload = (Rs*Zamp) / (Rs+Zamp) # parallel sum of Rs and Zamp
    
# 未知類比轉換器、假裝沒有 
Rload=0    

Rt = Rc + Rload # total resistance
bc = (Go**2)/(4*np.pi*fo*m*Rt) 
bt = bo + bc # total damping

p1 = (-2*np.pi*fo*bt) + 1j*(2*np.pi*fo*np.sqrt(1-bt**2))
p2 = (-2*np.pi*fo*bt) - 1j*(2*np.pi*fo*np.sqrt(1-bt**2))
z1 = 0.0 + 0j
z2 = 0.0 - 0j

Ao = 1
k1 = -Go
k2 = (2**(Res)-1)/5 # full scale = +/- 2.5V for 16 bit datalogger
const = Ao*k1

w = 2*np.pi*fd
Amp = np.abs( const* (1j*w-z1)*(1j*w-z2)/(1j*w-p1)/(1j*w-p2) )
Pha = np.angle( const* ((1j*w-z1)*(1j*w-z2))/((1j*w-p1)*(1j*w-p2)), deg=True)*(-1)


Fig, ax = plt.subplots(2, 1, figsize=(6,8))
ax[0].loglog(fd, Amp, 'r-')
ax[0].grid(which='major', axis='both', linewidth=0.75)
ax[0].grid(which='minor', axis='both', linewidth=0.75)
ax[0].set_ylim(1, 100)
ax[0].set_xlim(1, 1000)
ax[0].set_ylabel("$Amplitude (V/m/s)$", fontsize=14)

ax[1].semilogx(fd, Pha, 'k-')
ax[1].grid(which='major', axis='both', linewidth=0.75)
ax[1].grid(which='minor', axis='both', linewidth=0.75)
ax[1].set_ylim(0, 180)
ax[1].set_xlim(1, 1000)
ax[1].set_ylabel("$PhaseLag (degree)$", fontsize=14)
ax[1].set_xlabel("$Frequency (Hz)$", fontsize=14)
plt.show()
#plt.savefig('Instrument Response Figure_GS-32CT_In.pdf')
#plt.savefig('Instrument Response Figure_GS-20DX_In.pdf')