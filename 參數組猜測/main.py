# -*- coding: utf-8 -*-
"""
Created on Fri Jun 18 12:58:24 2021

@author: r8891
"""
import xlwings as xw
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.fft import fft

##需要手動輸入的資料
workbook= xw.Book("knockData_vertical.xlsx")
sheet1= workbook.sheets[0]
total_time= 134.2
dt= 0.002
sample_rate= 1/dt
sample_num= int(sample_rate*total_time)
time= np.linspace(0,total_time, sample_num, endpoint=False)
X=[1,4]
Y=[2,5]
Z=[3,6]
direction= X
freq_lowerlimit= 1
freq_upperlimit= 249

##電壓時域資料
Data_Geo1= sheet1.range((24,direction[0]),(24+sample_num-1,direction[0])).options(np.array).value
Data_Geo2= sheet1.range((24,direction[1]),(24+sample_num-1,direction[1])).options(np.array).value

plt.figure(1)
plt.suptitle("E(Time_Domain)")
down= min(min(Data_Geo1),min(Data_Geo2))
up= max(max(Data_Geo1),max(Data_Geo2))
plt.subplot(2,1,1)
plt.xlabel("Time(Sec)")
plt.ylabel("E(Volt)")
plt.ylim(down,up)
plt.plot(time,Data_Geo1)
plt.subplot(2,1,2)
plt.xlabel("Time(Sec)")
plt.ylabel("E(Volt)")
plt.ylim(down,up)
plt.plot(time,Data_Geo2)

##傅立葉轉換兩筆電壓時域資料to頻域資料
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
plt.ylim(down,up)
plt.plot(fd, amp_E_Geo1)
plt.subplot(2,1,2)
plt.xlabel("Freq(Hz)")
plt.ylabel("Amp_E(Volt)")
plt.ylim(down,up)
plt.plot(fd, amp_E_Geo2)

##需要手動輸入標準地聲儀器參數
##頻率響應函數
#地聲檢知器:GS20-DX
f0= 10 #自然頻率(Hz or 1/s)
a0= f0*2*np.pi #自然角頻率(rad/s)
m= 11 #(g)
G0= 27.6 #(V/m/s)
Rc= 395 #(Ohm)
Rs= -9999
b0= 0.3 #開路阻尼(%)

#類比訊號轉換器:PCI-1713U
Res = 12 # (bits)
MSR = 100*1000 # (S/s) = 100 (kS/s) max.
Zamp = 10**9 # 1 GΩ

#運算其他參數
if Rs == -9999:
    Rload = Zamp
else:
    Rload = (Rs*Zamp) / (Rs+Zamp) # parallel sum of Rs and Zamp

Rt = Rc + Rload #總電阻(Ohm)
bc = (G0**2)/(2*a0*m*Rt) 
bt = b0 + bc # total damping

#地聲儀器響應標準規範寫法
p1= (-a0*bt) + 1j*(a0*np.sqrt(1-bt**2))
p2= (-a0*bt) - 1j*(a0*np.sqrt(1-bt**2))
z1= 0.0 + 0j
z2= 0.0 - 0j
A0= 1
k1= -G0
k2= 1
const= A0*k1*k2
Amp= np.abs( const* (1j*wd-z1)*(1j*wd-z2)/(1j*wd-p1)/(1j*wd-p2) )
Pha= np.angle( const* ((1j*wd-z1)*(1j*wd-z2))/((1j*wd-p1)*(1j*wd-p2)), deg=True)
amp_T_Geo1= Amp

plt.figure(3)
plt.title("FreqResponse_Function")
plt.loglog(basex=10,basey=10)
plt.xlim(1,1000)
plt.ylim(1,100)
plt.xlabel("Freq(Hz)")
plt.ylabel("Amp_T(Volt/m/s)")
plt.plot(fd,Amp)
plt.grid(which='major', axis='both', linewidth=0.75)
plt.grid(which='minor', axis='both', linewidth=0.75)

##地表震動速度頻域資料
#df 是上面所有地聲資料(總整理)
df= pd.DataFrame(fd,columns=['freq(1/sec)'])
df['anglefreq(rad/sec)']= wd
df['amp_E_Geo1(volt)']= amp_E_Geo1
df['amp_T_Geo1(volt/m/s)']= amp_T_Geo1
df['amp_V_Geo1(m/s)']= df['amp_E_Geo1(volt)']/df['amp_T_Geo1(volt/m/s)']
df['amp_E_Geo2(volt)']= amp_E_Geo2
amp_V_Geo1= df['amp_V_Geo1(m/s)'].values

#第一項沒有值，直接另成0
amp_V_Geo1[0]=0

plt.figure(4)
plt.title("Velocity(Freq_Domain)")
down= min(amp_V_Geo1)
up= max(amp_V_Geo1)
plt.plot(fd,amp_V_Geo1)
plt.xlabel("Freq(Hz)")
plt.ylabel("Amp_V(m/s)")
plt.ylim(down,0.00015)

##擷取要分析的頻率段資料，傳到excel
#df_select是篩選過後的所有地聲資料
df_select= df[(df['freq(1/sec)'] >= freq_lowerlimit) & (df['freq(1/sec)'] <= freq_upperlimit)]
df_select.set_index('freq(1/sec)', inplace = True)
df_select.to_excel("Result.xlsx",sheet_name='data')
wd= df_select["anglefreq(rad/sec)"].values
amp_E_Geo1= df_select["amp_E_Geo1(volt)"].values
amp_T_Geo1= df_select["amp_T_Geo1(volt/m/s)"].values
amp_V_Geo1= df_select["amp_V_Geo1(m/s)"].values
amp_E_Geo2= df_select["amp_E_Geo2(volt)"].values

##最陡梯度法(SteepestDescentMethod)猜待測地聲儀器參數組
G= G0
a= a0
b= 2*m*Rt
h0= b0

learning_rate= -0.9
guess_time= 10000
end_condition= 0
recordingBox= np.zeros((guess_time,9))
title= ["G","a","b","h0","z","下一步(G)","下一步(a)","下一步(b)","下一步(h0)"]

z=0
nextstep_G=0
nextstep_a=0
nextstep_b=0
nextstep_h0=0
#用矩陣運算
for n in range(guess_time):
    temp1= G**(-2)*a**4
    temp2= (-2)*G**(-2)*a**2+4*G**(-2)*a**2*h0**2+8*a*b**(-1)*h0+4*G**2*b**(-2)
    temp3= (-2)*G**(-3)*a**4
    temp4= 4*G**(-3)*a**2-8*G**(-3)*a**2*h0**2+8*G*b**(-2)
    temp5= 4*G**(-2)*a**3
    temp6= 4*G**(-2)*a+8*G**(-2)*a*h0**2+8*h0*b**(-1)
    temp7= (-8)*a*h0*b**(-2)-8*G**2*b**(-3)
    temp8= 8*G**(-2)*a**2*h0+8*a*b**(-1)

    z= ((amp_V_Geo1)**2-(2*amp_V_Geo1*amp_E_Geo2/(wd**2))*(temp1+temp2*wd**2+G**(-2)*wd**4)**(0.5)+(amp_E_Geo2**2/wd**4)*(temp1+temp2*wd**2+G**(-2)*wd**4)).sum()
    nextstep_G= learning_rate*((amp_E_Geo2**2/wd**4)*(temp3+temp4*wd**2-2*G**(-3)*wd**4)-(amp_V_Geo1*amp_E_Geo2/(wd**2))*(temp1+temp2*wd**2+G**(-2)*wd**4)**(-0.5)*(temp3+temp4*wd**2-2*G**(-3)*wd**4)).sum()
    nextstep_a= learning_rate*((amp_E_Geo2**2/wd**4)*(temp5+temp6*wd**2)-(amp_V_Geo1*amp_E_Geo2/(wd**2))*(temp1+temp2*wd**2+G**(-2)*wd**4)**(-0.5)*(temp5+temp6*wd**2)).sum()
    nextstep_b= learning_rate*((amp_E_Geo2**2/wd**4)*temp7-(amp_V_Geo1*amp_E_Geo2)*(temp1+temp2*wd**2+G**(-2)*wd**4)**(-0.5)*temp7).sum()  
    nextstep_h0= learning_rate*((amp_E_Geo2**2/wd**4)*temp8-(amp_V_Geo1*amp_E_Geo2)*(temp1+temp2*wd**2+G**(-2)*wd**4)**(-0.5)*temp8).sum()

    recordingBox[n,0]= G
    recordingBox[n,1]= a
    recordingBox[n,2]= b
    recordingBox[n,3]= h0
    recordingBox[n,4]= z
    recordingBox[n,5]= nextstep_G
    recordingBox[n,6]= nextstep_a
    recordingBox[n,7]= nextstep_b
    recordingBox[n,8]= nextstep_h0
    
    #走向下一步
    G+= nextstep_G
    a+= nextstep_a
    b+= nextstep_b
    h0+= nextstep_h0
    
print('Ok')

wb= xw.Book("Result.xlsx")
sheet= wb.sheets.add()
sheet.range("A1").value= title
sheet.range("A2").value= recordingBox