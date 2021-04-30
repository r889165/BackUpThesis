#-*- coding: utf-8 -*-
"""
Created on Thu Mar 18 13:20:39 2021

@author: r8891
"""
import xlrd
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft

##1.先處理一開始得到的資料(電壓時域圖)
#取樣頻率:(一秒取幾個電壓值)
sample_num= 25000
total_time= 5
sample_rate= sample_num/total_time
time=np.linspace(0, total_time, sample_num, endpoint=False)

correctFile= xlrd.open_workbook('correctData.xlsx')
correctSheet= correctFile.sheet_by_index(0)
correctData= correctSheet.col_values(0,0,sample_num)
correctData= np.array(correctData)

testFile= xlrd.open_workbook('testData.xls')
testSheet=  testFile.sheet_by_index(0)
testData= testSheet.col_values(0,0,sample_num)
testData= np.array(testData)

print(max(correctData))
print(max(testData))

plt.figure(1)
plt.suptitle("E(TimeDomain)")
plt.subplot(2,1,1)
plt.xlabel("time(s)")
plt.ylabel("voltage(V)")
plt.plot(time,correctData)
plt.subplot(2,1,2)
plt.xlabel("time(s)")
plt.ylabel("voltage(V)")
plt.plot(time,testData)



##2.電壓時域圖轉頻域圖、使用scipy.fft套件
#要注意Nyquist效應(對稱答案有一半會不見)。基於能量不變，所以振幅全部都乘兩倍
correctData_fft= fft(correctData)
testData_fft= fft(testData)
fd = np.linspace(0.0, int(sample_rate), int(sample_num), endpoint=False)

E_correctData= 2/sample_num *correctData_fft[0:int(sample_num/2)]
E_testData= 2/sample_num *testData_fft[0:int(sample_num/2)]
fd=fd[0:int(sample_num/2)]

amp_E_correctData= np.abs(E_correctData)
amp_E_testData= np.abs(E_testData)
print(max(amp_E_correctData))
print(max(amp_E_testData))
print("確認電壓最大振福時域頻域答案不會差太多")
print("---------------------------")

plt.figure(2)
plt.suptitle("E(W)")
plt.subplot(2,1,1)
plt.xlabel("freq(Hz)")
plt.ylabel("amp(V)")
plt.plot(fd, amp_E_correctData)
plt.subplot(2,1,2)
plt.xlabel("freq(Hz)")
plt.ylabel("amp(V)")
plt.plot(fd, amp_E_testData)


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


##4.過濾想要的頻率區間，去轉成速度頻域圖
#範例:correct跟test會一樣
freq_lowerlimit= 1
freq_upperlimit= 249


V_correctData=[0]*len(fd)
for i in range(len(fd)):
    if fd[i]<freq_upperlimit and fd[i]>freq_lowerlimit:
        V_correctData[i]=E_correctData[i]/(const*(1j*wd[i]-z1)*(1j*wd[i]-z2)/(1j*wd[i]-p1)/(1j*wd[i]-p2))
    else:
        V_correctData[i]=0

V_testData=V_correctData
      
amp_V_correctData= np.abs(V_correctData)
amp_V_testData= np.abs(V_testData)
 
plt.figure(4)
plt.suptitle("velocity(W)")
plt.subplot(2,1,1)
plt.xlabel("freq(Hz)")
plt.ylabel("amp(in/s)")
plt.xlim(0,freq_upperlimit)
plt.ylim(0,0.0005)
plt.plot(fd, amp_V_correctData)

plt.subplot(2,1,2)
plt.xlabel("freq(Hz)")
plt.ylabel("amp(in/s)")
plt.xlim(0,freq_upperlimit)
plt.ylim(0,0.0005)
plt.plot(fd, amp_V_testData)


            
##5.得到待測儀器的響應函數
#畫正確跟待測儀器兩者的響應函數
T_testData=[0]*len(fd)
for i in range(len(fd)):
    if V_testData[i]==0:
        T_testData[i]=0
    else:
        T_testData[i]= E_testData[i]/V_testData[i]

T_correctData=[0]*len(fd)
for i in range(len(fd)):
    if V_testData[i]==0:
        T_correctData[i]=0
    else:
        T_correctData[i]=(const*(1j*wd[i]-z1)*(1j*wd[i]-z2)/(1j*wd[i]-p1)/(1j*wd[i]-p2))

amp_T_correctData= np.abs(T_correctData)
amp_T_testData= np.abs(T_testData)
gap=amp_T_correctData-amp_T_testData

plt.figure(5)
plt.suptitle("response(W)")
plt.subplot(2,1,1)
plt.xlabel("freq(Hz)")
plt.ylabel("amp")
plt.xlim(0,freq_upperlimit)
plt.plot(fd, amp_T_correctData)
plt.subplot(2,1,2)
plt.xlabel("freq(Hz)")
plt.ylabel("amp")
plt.xlim(0,freq_upperlimit)
plt.plot(fd, amp_T_testData)

plt.figure(6)
plt.xlim(0,freq_upperlimit)
plt.xlabel("freq(Hz)")
plt.ylabel("相減得到")
plt.plot(fd,gap)


##6.迴歸求儀器參數(擷取T_testData有值的位置、跟對應的頻率(W))
#a+bi的倒數=>(a/(a^2+b^2))+(-b/(a^2+b^2))i
for i in range(len(fd)):
    if T_testData[i]!=0:
        start=i
        break
end=1
for i in range(start,len(fd)):
    if T_testData[i]==0:
        end=i
        break

T_testData_select=T_testData[start:end]
a=np.real(T_testData_select)
b=np.imag(T_testData_select)
W= wd[start:end]
Fre= fd[start:end]

#(1)簡單回歸:Y1=b0+b1*X1
Y1=a/(a**2+b**2)
X1=1/W**2

#(2)強制過原點簡單迴歸:Y2=0+b2*X2
Y2=-b/(a**2+b**2)
X2=1/W

#(3)做振幅的回歸: Y3=b3+b4*X3+b5*X3^2
Y3=1/(a**2+b**2)
X3=X1


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
sheet1.write(0,6,"X3")
sheet1.write(0,7,"Y3")

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

k=1
for i in X3:
    sheet1.write(k,6,i)
    k=k+1   
    
k=1
for i in Y3:
    sheet1.write(k,7,i)
    k=k+1   

    
linearRegressionData.save('linearRegressionData.xls')






    

    
    


    



    





