# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 16:09:31 2021

@author: r8891
"""

import xlwings as xw
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
from scipy.fft import fft

#輸入的資訊
workbook1= xw.Book("mudflow.xlsx")
sheet1= workbook1.sheets[0]
start_time=80
end_time=250
dt=0.002

X=[1,4]
Y=[2,5]
Z=[3,6]
direction= Z

#其他資料
start_index= round(24+(start_time/dt)-1)
sample_rate=1/dt
total_time= end_time-start_time
sample_num= round(((end_time-start_time)/dt))+1
time= np.linspace(start_time,end_time,sample_num,endpoint=True)


Data_Geo1= sheet1.range((start_index,direction[0]),(start_index+sample_num-1,direction[0])).options(np.array).value
Data_Geo1_fft= fft(Data_Geo1)

fd = np.linspace(0, int(sample_rate), sample_num, endpoint=False)
amp_E_Geo1= np.abs(2/sample_num *Data_Geo1_fft[0:int(sample_num/2)])
fd=fd[0:int(sample_num/2)]
wd=fd*2*np.pi

#標準地聲檢知器:GS20-DX
f0= 10 #自然頻率(Hz or 1/s)
a0= f0*2*np.pi #自然角頻率(rad/s)
m= 0.011 #(Kg)
G0= 27.6 #(V/m/s)
Rc= 395 #(Ohm)
h00= 0.3 #開路阻尼(%)

#待測地聲檢知器:GS20-DX
f1= 20
a1= f1*2*np.pi
m=0.011
G1= 40
Rc= 395
h01= 0.41


#類比訊號轉換器:PCI-1713U
Rs= -9999
Zamp= 10**9 # 1 GΩ
k2= 1 #加倍係數(count/Volt)

if Rs == -9999:
    Rload = Zamp
else:
    Rload = (Rs*Zamp) / (Rs+Zamp) 

Rt = Rc + Rload #總電阻(Ohm)
bc = (G0**2)/(2*a0*m*Rt) 
bt = h00+ bc # total damping

T_Geo1=G0*wd**2/(a0**2+2*1j*(h00+G0**2/(2*m*Rt*a0))*a0*wd-wd**2)
T_Geo2=G1*wd**2/(a1**2+2*1j*(h01+G1**2/(2*m*Rt*a1))*a1*wd-wd**2)
amp_T_Geo1= abs(T_Geo1)
amp_T_Geo2= abs(T_Geo2)

#df 是上面所有地聲資料(總整理)
df= pd.DataFrame(fd,columns=['1頻率(1/sec)'])
df['2角頻率(rad/sec)']= wd
df['3標準儀器:電壓(volt)']= amp_E_Geo1
df['4標準儀器:振幅響應(volt/m/s)']= amp_T_Geo1
df['5標準儀器:震動速度(m/s)']= df['3標準儀器:電壓(volt)']/df['4標準儀器:振幅響應(volt/m/s)']
df['5標準儀器:震動速度(m/s)'].iloc[0]=0
df['6待測儀器:振幅響應(volt/m/s)']= amp_T_Geo2

amp_V_Geo1= df['5標準儀器:震動速度(m/s)'].values
df['7待測儀器:電壓(volt)']= amp_T_Geo2*amp_V_Geo1

##輸入值
p1= 2
p2= 5
p3= 40
p4= 100

n1=20
n2=5
n3=20


##df_select1、2、3分成三段
df_select1= df[(df['1頻率(1/sec)'] >= p1) & (df['1頻率(1/sec)'] <= p2)]
df_select2= df[(df['1頻率(1/sec)'] > p2) & (df['1頻率(1/sec)'] < p3)]
df_select3= df[(df['1頻率(1/sec)'] >= p3) & (df['1頻率(1/sec)'] <= p4)]

df_select1.reset_index(inplace=True, drop=True)
df_select2.reset_index(inplace=True, drop=True)
df_select3.reset_index(inplace=True, drop=True)

#第一段
fd1=np.linspace(math.log(p1,10),math.log((p2),10),n1,endpoint=True)
fd1= (10**fd1)
choose1=list(range(n1))

for i in range(len(fd1)):
    choose1[i]= int((len(df_select1)-1)*(fd1[i]-p1)/(p2-p1)) 

df_1=pd.DataFrame()
for i in range(len(choose1)):
    df_1= df_1.append(df_select1.iloc[choose1[i]],ignore_index=True)

#第二段
choose2= list(range(n2*int(p3-p2)))
interval= len(df_select2)/((p3-p2)*n2)

for i in range(len(choose2)):
    choose2[i]=math.floor(i*interval)

df_2=pd.DataFrame()
for i in range(len(choose2)):
    df_2= df_2.append(df_select2.iloc[choose2[i]],ignore_index=True)

#第三段
fd3=np.linspace(math.log(p3,10),math.log((p4),10),n3,endpoint=True)
fd3= (10**fd3)
choose3=list(range(n3))

for i in range(len(fd3)):
    choose3[i]= int((len(df_select3)-1)*(fd3[i]-p3)/(p4-p3)) 

df_3=pd.DataFrame()
for i in range(len(choose3)):
    df_3= df_3.append(df_select3.iloc[choose3[i]],ignore_index=True)

df_select=df_1.append(df_2)
df_select=df_select.append(df_3)

df_select.to_excel("Result.xlsx",sheet_name='data')

fd= df_select["1頻率(1/sec)"].values
wd= df_select["2角頻率(rad/sec)"].values
amp_E_Geo1= df_select["3標準儀器:電壓(volt)"].values
amp_T_Geo1= df_select["4標準儀器:振幅響應(volt/m/s)"].values
amp_V_Geo1= df_select["5標準儀器:震動速度(m/s)"].values
amp_T_Geo2= df_select["6待測儀器:振幅響應(volt/m/s)"].values
amp_E_Geo2= df_select["7待測儀器:電壓(volt)"].values
num=len(df_select)

#正規化
W0=100
V0=10**(-7)
E0=10**(-5)
b0=2*m*Rt

wd= wd/W0
amp_V_Geo1= amp_V_Geo1/V0
amp_E_Geo2= amp_E_Geo2/E0

const1= round(E0**2*a0**4/(W0**4*V0**2*G0**2),3)
const2= round((-2)*E0**2*a0**2/(W0**2*V0**2*G0**2),3)
const3= round(4*E0**2*a0**2*h00**2/(W0**2*V0**2*G0**2),3)
const4= 8*E0**2*a0*h00/(W0**2*V0**2*b0)
const5= 4*E0**2*G0**2/(W0**2*V0**2*b0**2)
const6= round(E0**2/(V0**2*G0**2),3)


G= 1.5
a= 1.5
h0= 1.5
learning_rate= -0.0001
learning_rate_h0= 1*learning_rate

guess_time= 30000
end_condition= 0
recordingBox= np.zeros((guess_time,10))
title= ["G","a","h0","Z","dZdG","dZda","dZdh0","下一步G","下一步a","下一步h0"]


Z=0
nextstep_G=0
nextstep_a=0
nextstep_h0=0

#用矩陣運算
for n in range(guess_time):
    
    dZdG=0
    dZda=0
    dZdh0=0
    
    
    temp1= const1*a**4*G**(-2)
    temp2= const2*a**2*G**(-2)+const3*a**2*h0**2*G**(-2)
    temp3= const6*G**(-2)
    
    temp4= const1*(-2)*a**4*G**(-3)
    temp5= const2*(-2)*a**2*G**(-3)+const3*(-2)*a**2*h0**2*G**(-3)
    
    temp6= const6*(-2)*G**(-3)
    temp7= const1*4*a**3*G**(-2)
    
    temp8= const2*2*a*G**(-2)+const3*2*a*h0**2*G**(-2)
    temp9= const3*2*a**2*h0*G**(-2)
   

    Z= (((amp_V_Geo1-amp_E_Geo2*(temp1+temp2*wd**2+temp3*wd**4)**0.5/(wd**2))**2).sum())/num
    
    
    dZdG= ((amp_E_Geo2**2/(wd**4)*(temp4+temp5*wd**2+temp6*wd**4)-(amp_V_Geo1*amp_E_Geo2/(wd**2))*(temp1+temp2*wd**2+temp3*wd**4)**(-0.5)*(temp4+temp5*wd**2+temp6*wd**4)).sum())/num
    dZda= ((amp_E_Geo2**2/(wd**4)*(temp7+temp8*wd**2)-(amp_V_Geo1*amp_E_Geo2/(wd**2))*(temp1+temp2*wd**2+temp3*wd**4)**(-0.5)*(temp7+temp8*wd**2)).sum())/num
    dZdh0= ((amp_E_Geo2**2/(wd**2)*temp9-amp_V_Geo1*amp_E_Geo2*(temp1+temp2*wd**2+temp3*wd**4)**(-0.5)*temp9).sum())/num
    
    
    nextstep_G= learning_rate*dZdG
    nextstep_a= learning_rate*dZda
    nextstep_h0= learning_rate_h0*dZdh0

    recordingBox[n,0]= G
    recordingBox[n,1]= a
    recordingBox[n,2]= h0
    recordingBox[n,3]= Z
    
    recordingBox[n,4]= dZdG
    recordingBox[n,5]= dZda
    recordingBox[n,6]= dZdh0
    
    recordingBox[n,7]= nextstep_G
    recordingBox[n,8]= nextstep_a
    recordingBox[n,9]= nextstep_h0

    #走向下一步
    G+= nextstep_G
    a+= nextstep_a
    h0+= nextstep_h0
    
print('Ok')

wb= xw.Book("Result.xlsx")
sheet= wb.sheets.add()
sheet.range("A1").value= title
sheet.range("A2").value= recordingBox


G=G0*G
a=a0*a
h0=h00*h0

fd_x=np.linspace(0,250,2000)
wd_x=fd_x*2*np.pi

T1=G0*wd_x**2/(a0**4-2*a0**2*wd_x**2+wd_x**4+4*(h00+G0**2/(a0*b0))**2*a0**2*wd_x**2)**0.5
T2=G*wd_x**2/(a**4-2*a**2*wd_x**2+wd_x**4+4*(h0+G**2/(a*b0))**2*a**2*wd_x**2)**0.5

plt.figure(5)
plt.title("FreqResponse_Function")
plt.loglog(basex=10,basey=10)
plt.xlim(1,1000)
plt.ylim(1,100)
plt.xlabel("Freq(Hz)")
plt.ylabel("Sensitivity(Volt/m/s)")

#標準
plt.plot(fd_x,T1,color='black')

#迴歸出來的線
plt.plot(fd_x,T2,color='red')

#取樣的點
plt.scatter(fd,amp_T_Geo2)

plt.grid(which='major', axis='both', linewidth=0.75)
plt.grid(which='minor', axis='both', linewidth=0.75)





