# -*- coding: utf-8 -*-
"""
Created on Mon Jun 28 21:24:42 2021

@author: r8891
"""
import numpy as np
import matplotlib.pyplot as plt
import math


G0=27.6
a0=62.8
h00=0.3
b=22000008690

G= np.zeros((11,1))
a= np.zeros((11,1))
h0= np.zeros((11,1))

ansbox= np.zeros((11,1))


for i in range(11):
    G[i]=G0*(0.5+0.1*i)
    a[i]=a0*(0.5+0.1*i)
    h0[i]=h00*(0.5+0.1*i)
    



#全部點
fd=np.linspace(0,250,2000)
wd=fd*2*np.pi

#分隔點
p1=2
p2=5
p3=40
p4=100
fds=np.array([p1,p2,p3,p4])
wds=fds*2*np.pi

#第一段
fd1=np.linspace(math.log(p1,10),math.log(p2,10),5)
fd1=10**fd1
#第二段
fd2=np.arange(p2+1,p3,0.5)
#第三段
fd3=np.linspace(math.log(p3,10),math.log(p4,10),5)
fd3=10**fd3
#合併
fd_select= np.concatenate([fd1,fd2,fd3])
wd_select= fd_select*2*np.pi


num= len(fd_select)




for i in range(0,11):
    Z=0
    j=5
    k=5
    
    #線
    T=G[i]*wd**2/(a[j]**4-2*a[j]**2*wd**2+wd**4+4*(h0[k]+G[i]**2/(a[j]*b))**2*a[j]**2*wd**2)**0.5
    
    #分隔點
    Ts=G[i]*wds**2/(a[j]**4-2*a[j]**2*wds**2+wds**4+4*(h0[k]+G[i]**2/(a[j]*b))**2*a[j]**2*wds**2)**0.5
    
    #三段的所有點
    T_select=G[i]*wd_select**2/(a[j]**4-2*a[j]**2*wd_select**2+wd_select**4+4*(h0[k]+G[i]**2/(a[j]*b))**2*a[j]**2*wd_select**2)**0.5
    T_correct=G0*wd_select**2/(a0**4-2*a0**2*wd_select**2+wd_select**4+4*(h00+G0**2/(a0*b))**2*a0**2*wd_select**2)**0.5
    
    Z=((T_select-T_correct)**2).sum()/num
    ansbox[i]=Z
    
    
    plt.title("FreqResponse_Function")
    plt.loglog(basex=10,basey=10)
    plt.xlim(1,1000)
    plt.ylim(1,100)
    plt.xlabel("Freq(Hz)")
    plt.ylabel("Amp_T(Volt/m/s)")
    
    plt.plot(fd,T)
    

    plt.scatter(fds,Ts,color="red", s=14)
    
    plt.grid(which='major', axis='both', linewidth=0.75)
    plt.grid(which='minor', axis='both', linewidth=0.75)
    
ansbox=ansbox.T
    
    
    



