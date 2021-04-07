#-*- coding: utf-8 -*-
"""
Created on Thu Mar 18 13:20:39 2021

@author: r8891
"""
import xlrd
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft

##1.先處理一開始得到的資料
#輸入必要資訊
File= xlrd.open_workbook('knock_vertical.xlsx')
Sheet= File.sheet_by_index(0)

total_time= 134
dt=0.002
sample_rate=1/dt
sample_num= int(sample_rate*total_time)
time=np.linspace(0, total_time, sample_num, endpoint=False)

#輸入要分析哪個方向的震動
X_direction=[0,3,6]
Y_direction=[1,4,7]
Z_direction=[2,5,8]
directionYouWantToAnaylze= Z_direction

Geo1Data= Sheet.col_values(directionYouWantToAnaylze[0],23,23+sample_num)
Geo1Data= np.array(Geo1Data)
Geo2Data= Sheet.col_values(directionYouWantToAnaylze[1],23,23+sample_num)
Geo2Data= np.array(Geo2Data)
Geo3Data= Sheet.col_values(directionYouWantToAnaylze[2],23,23+sample_num)
Geo3Data= np.array(Geo3Data)

#輸入要分析的頻率區段
freq_lowerlimit= 20
freq_upperlimit= 200


#電壓時域圖
plt.figure(1)
plt.suptitle("E(time)")
down= min(min(Geo1Data),min(Geo2Data),min(Geo3Data))
up= max(max(Geo1Data),max(Geo2Data),max(Geo3Data))

plt.subplot(3,1,1)
plt.xlabel("time(s)")
plt.ylabel("voltage(V)")
plt.ylim(down,up)
plt.plot(time,Geo1Data)

plt.subplot(3,1,2)
plt.xlabel("time(s)")
plt.ylabel("voltage(V)")
plt.ylim(down,up)
plt.plot(time,Geo2Data)

plt.subplot(3,1,3)
plt.xlabel("time(s)")
plt.ylabel("voltage(V)")
plt.ylim(down,up)
plt.plot(time,Geo3Data)


##2.電壓時域圖轉頻域圖、使用scipy.fft套件
#要注意Nyquist效應(對稱答案有一半會不見)。基於能量不變，所以振幅全部都乘兩倍
Geo1Data_fft= fft(Geo1Data)
Geo2Data_fft= fft(Geo2Data)
Geo3Data_fft= fft(Geo3Data)
fd = np.linspace(0.0, int(sample_rate), int(sample_num), endpoint=False)

E_Geo1Data= 2/sample_num *Geo1Data_fft[0:int(sample_num/2)]
E_Geo2Data= 2/sample_num *Geo2Data_fft[0:int(sample_num/2)]
E_Geo3Data= 2/sample_num *Geo3Data_fft[0:int(sample_num/2)]
fd=fd[0:int(sample_num/2)]

amp_E_Geo1Data= np.abs(E_Geo1Data)
amp_E_Geo2Data= np.abs(E_Geo2Data)
amp_E_Geo3Data= np.abs(E_Geo3Data)


plt.figure(2)
plt.suptitle("E(freq)")
down= min(min(amp_E_Geo1Data),min(amp_E_Geo2Data),min(amp_E_Geo3Data))
up= max(max(amp_E_Geo1Data),max(amp_E_Geo2Data),max(amp_E_Geo3Data))

plt.subplot(3,1,1)
plt.xlabel("freq(Hz)")
plt.ylabel("amp(V)")
plt.ylim(down,up)
plt.plot(fd, amp_E_Geo1Data)

plt.subplot(3,1,2)
plt.xlabel("freq(Hz)")
plt.ylabel("amp(V)")
plt.ylim(down,up)
plt.plot(fd, amp_E_Geo2Data)

plt.subplot(3,1,3)
plt.xlabel("freq(Hz)")
plt.ylabel("amp(V)")
plt.ylim(down,up)
plt.plot(fd, amp_E_Geo3Data)

##3.選擇準確的儀器編號:求響應函數
#範例:GS-20DX
fo = 10 # (Hz)
fn = 200 # (Hz)
m = 11 # (g) 
Go = 27.6 # (V/m/s)
Rc = 395 # (Ohm)
Rs = -9999
bo = 0.3 # = 30 (%)

#類比訊號轉換器:PCI-1713U
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

Ao= 1
k1 = -Go
k2 = 1
const = Ao*k1*k2
wd= 2*np.pi*fd
Amp = np.abs( const* (1j*wd-z1)*(1j*wd-z2)/(1j*wd-p1)/(1j*wd-p2) )
Pha = np.angle( const* ((1j*wd-z1)*(1j*wd-z2))/((1j*wd-p1)*(1j*wd-p2)), deg=True)

plt.figure(3)
Fig, ax = plt.subplots(2, 1, figsize=(6,8))
ax[0].loglog(fd, Amp, 'r-')
ax[0].grid(which='major', axis='both', linewidth=0.75)
ax[0].grid(which='minor', axis='both', linewidth=0.75)
ax[0].set_ylim(1, 100)
ax[0].set_xlim(1, 1000)
ax[0].set_ylabel("$Amplitude (V/m/s)$", fontsize=14)

ax[1].semilogx(fd, Pha*(-1), 'k-')
ax[1].grid(which='major', axis='both', linewidth=0.75)
ax[1].grid(which='minor', axis='both', linewidth=0.75)
ax[1].set_ylim(0, 180)
ax[1].set_xlim(1, 1000)
ax[1].set_ylabel("$PhaseLag (degree)$", fontsize=14)
ax[1].set_xlabel("$Frequency (Hz)$", fontsize=14)
plt.show()


##4.過濾想要的頻率區間，去轉成速度頻域圖
V_Geo1Data=[0]*len(fd)
for i in range(len(fd)):
    if fd[i]<freq_upperlimit and fd[i]>freq_lowerlimit:
        V_Geo1Data[i]=E_Geo1Data[i]/(const*(1j*wd[i]-z1)*(1j*wd[i]-z2)/(1j*wd[i]-p1)/(1j*wd[i]-p2))
    else:
        V_Geo1Data[i]=0
V_Geo2Data=V_Geo1Data
V_Geo3Data=V_Geo1Data

amp_V_Geo1Data= np.abs(V_Geo1Data)
amp_V_Geo2Data= np.abs(V_Geo2Data)
amp_V_Geo3Data= np.abs(V_Geo3Data)


plt.figure(4)
plt.suptitle("velocity(freq)")
down= min(min(amp_V_Geo1Data),min(amp_V_Geo2Data),min(amp_V_Geo3Data))
up= max(max(amp_V_Geo1Data),max(amp_V_Geo2Data),max(amp_V_Geo3Data))

plt.subplot(3,1,1)
plt.xlabel("freq(Hz)")
plt.ylabel("amp(in/s)")
plt.xlim(0,freq_upperlimit)
plt.ylim(down,up)
plt.plot(fd, amp_V_Geo1Data)

plt.subplot(3,1,2)
plt.xlabel("freq(Hz)")
plt.ylabel("amp(in/s)")
plt.xlim(0,freq_upperlimit)
plt.ylim(down,up)
plt.plot(fd, amp_V_Geo2Data)

plt.subplot(3,1,3)
plt.xlabel("freq(Hz)")
plt.ylabel("amp(in/s)")
plt.xlim(0,freq_upperlimit)
plt.ylim(down,up)
plt.plot(fd, amp_V_Geo3Data)


##5.得到待測儀器的響應函數
#畫正確跟待測儀器兩者的響應函數
T_Geo2Data=[0]*len(fd)
for i in range(len(fd)):
    if V_Geo2Data[i]==0:
        T_Geo2Data[i]=0
    else:
        T_Geo2Data[i]= E_Geo2Data[i]/V_Geo2Data[i]

T_Geo1Data=[0]*len(fd)
for i in range(len(fd)):
    if V_Geo2Data[i]==0:
        T_Geo1Data[i]=0
    else:
        T_Geo1Data[i]=(const*(1j*wd[i]-z1)*(1j*wd[i]-z2)/(1j*wd[i]-p1)/(1j*wd[i]-p2))
        
T_Geo3Data=[0]*len(fd)
for i in range(len(fd)):
    if V_Geo3Data[i]==0:
        T_Geo3Data[i]=0
    else:
        T_Geo3Data[i]= E_Geo3Data[i]/V_Geo3Data[i]

amp_T_Geo1Data= np.abs(T_Geo1Data)
amp_T_Geo2Data= np.abs(T_Geo2Data)
amp_T_Geo3Data= np.abs(T_Geo3Data)

gap=amp_T_Geo2Data-amp_T_Geo1Data

plt.figure(5)
plt.suptitle("response(freq)")
down= min(min(amp_T_Geo1Data),min(amp_T_Geo2Data),min(amp_T_Geo3Data))
up= (min(max(amp_T_Geo1Data),max(amp_T_Geo2Data),max(amp_T_Geo3Data))+max(max(amp_T_Geo1Data),max(amp_T_Geo2Data),max(amp_T_Geo3Data)))/2

plt.subplot(3,1,1)
plt.xlabel("freq(Hz)")
plt.ylabel("amp")
plt.xlim(0,freq_upperlimit)
plt.ylim(down,up)
plt.plot(fd, amp_T_Geo1Data)

plt.subplot(3,1,2)
plt.xlabel("freq(Hz)")
plt.ylabel("amp")
plt.xlim(0,freq_upperlimit)
plt.ylim(down,up)
plt.plot(fd, amp_T_Geo2Data)

plt.subplot(3,1,3)
plt.xlabel("freq(Hz)")
plt.ylabel("amp")
plt.xlim(0,freq_upperlimit)
plt.ylim(down,up)
plt.plot(fd, amp_T_Geo3Data)



plt.figure(6)
plt.title("$G1-G2(ResponseFreq)$")
plt.xlim(0,freq_upperlimit)
plt.xlabel("freq(Hz)")
plt.ylabel("G2-G1")
plt.plot(fd,gap)


##6.迴歸求儀器參數(擷取T_testData有值的位置、跟對應的頻率(W))
#a+bi的倒數=>(a/(a^2+b^2))+(-b/(a^2+b^2))i
start=0
for i in range(len(fd)):
    if T_Geo2Data[i]!=0:
        start=i
        break
end=1
for i in range(start,len(fd)):
    if T_Geo2Data[i]==0:
        end=i
        break

T_Geo2Data_select=T_Geo2Data[start:end]
a=np.real(T_Geo2Data_select)
b=np.imag(T_Geo2Data_select)
W=wd[start:end]
Fre=fd[start:end]

#(1)簡單回歸:Y1=b0+b1*X1
Y1=a/(a**2+b**2)
X1=1/W**2

#(2)強制過原點簡單迴歸:Y2=0+b2*X2
Y2=-b/(a**2+b**2)
X2=1/W


#迴歸資料寫進新的xls檔案裡面
import xlwt
linearRegressionData= xlwt.Workbook()
sheet1= linearRegressionData.add_sheet('data')

sheet1.write(0,0,"Freq")
sheet1.write(0,1,"AngularFreq")
sheet1.write(0,2,"X1")
sheet1.write(0,3,"Y1")
sheet1.write(0,4,"X2")
sheet1.write(0,5,"Y2")

k=1
for i in Fre:
    sheet1.write(k,0,i)
    k=k+1
k=1
for i in W:
    sheet1.write(k,1,i)
    k=k+1
k=1
for i in X1:
    sheet1.write(k,2,i)
    k=k+1
k=1
for i in Y1:
    sheet1.write(k,3,i)
    k=k+1
k=1
for i in X2:
    sheet1.write(k,4,i)
    k=k+1
k=1
for i in Y2:
    sheet1.write(k,5,i)
    k=k+1
    
linearRegressionData.save('linearRegressionData.xls')


    


plt.figure(100)
a=amp_E_Geo1Data/Amp
plt.plot(fd,a)
plt.xlim(0,200)
plt.ylim(0,10**(-5))

plt.figure(101)
b=amp_E_Geo2Data/amp_V_Geo2Data
plt.plot(fd,b)
plt.xlim(0,200)
plt.ylim(0,2500)

    

    
    


    



    





