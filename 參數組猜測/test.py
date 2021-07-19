# -*- coding: utf-8 -*-

import xlwings as xw
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
from scipy.fft import fft

workbook1= xw.Book("knockData_vertical.xlsx")
sheet1= workbook1.sheets[0]
start_time=0.002
end_time=134.2
dt=0.002

start_index= round(24+(start_time/dt)-1)
sample_rate=1/dt
total_time= end_time-start_time
sample_num= round(((end_time-start_time)/dt))+1
time= np.linspace(start_time,end_time,sample_num,endpoint=True)

X=[1,4]
Y=[2,5]
Z=[3,6]
direction= Z

Data_Geo1= sheet1.range((start_index,direction[0]),(start_index+sample_num-1,direction[0])).options(np.array).value
Data_Geo2= sheet1.range((start_index,direction[1]),(start_index+sample_num-1,direction[1])).options(np.array).value


plt.figure(1)
plt.suptitle("E(Time_Domain)")
down= min(min(Data_Geo1),min(Data_Geo2))
up= max(max(Data_Geo1),max(Data_Geo2))

plt.subplot(2,1,1)
plt.xlabel("Time(Sec)")
plt.ylabel("E(Volt)")
plt.ylim(down,up)
plt.plot(time,Data_Geo1,color='black')
plt.subplot(2,1,2)
plt.xlabel("Time(Sec)")
plt.ylabel("E(Volt)")
plt.ylim(down,up)
plt.plot(time,Data_Geo2,color='black')


Data_Geo1_fft= fft(Data_Geo1)
Data_Geo2_fft= fft(Data_Geo2)
fd = np.linspace(0, int(sample_rate), sample_num, endpoint=False)

amp_E_Geo1= np.abs(2/sample_num *Data_Geo1_fft[0:int(sample_num/2)])
amp_E_Geo2= np.abs(2/sample_num *Data_Geo2_fft[0:int(sample_num/2)])
fd=fd[0:int(sample_num/2)]
wd=fd*2*np.pi

plt.figure(2)
plt.suptitle("E(Freq_Domain)")
down= min(min(amp_E_Geo1),min(amp_E_Geo2))
up= max(max(amp_E_Geo1),max(amp_E_Geo2))

plt.subplot(2,1,1)
plt.xlabel("Freq(Hz)")
plt.ylabel("Amp_E(Volt)")
plt.xlim(0,250)
plt.ylim(down,up)
plt.plot(fd, amp_E_Geo1,color='black')
plt.subplot(2,1,2)
plt.xlabel("Freq(Hz)")
plt.ylabel("Amp_E(Volt)")
plt.xlim(0,250)
plt.ylim(down,up)
plt.plot(fd, amp_E_Geo2,color='black')