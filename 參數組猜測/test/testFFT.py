# -*- coding: utf-8 -*-
"""
Created on Sun Jul 18 11:35:13 2021

@author: r8891
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.fft import fft

#設計一個紀錄一秒收集1400個值(頻率:1400Hz)
#已知是三個sin波疊加，週期分別是180、390、600
x=np.linspace(0,1,1400)

y=7*np.sin(2*np.pi*180*x)+2.8*np.sin(2*np.pi*390*x)+5.1*np.sin(2*np.pi*600*x)
yy=fft(y) 

fd = np.linspace(0, 1400, 1400, endpoint=False)
amp_y= np.abs(2/1400 *yy[0:int(1400/2)])
fd= fd[0:int(1400/2)]
             
plt.plot(fd,amp_y)