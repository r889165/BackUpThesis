# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft
import xlrd

##1.先寫好SDM的方法
##這個在這邊沒有用到
def gradient_descent(gradient, start, learn_rate, n_iter=50, tolerance= 1e-06):
    vector= start
    for _ in range(n_iter):
        diff= -learn_rate * gradient(vector)
        if np.all(np.abs(diff) <= tolerance):
            break
        vector+= diff
    return vector

#gradient: python內建函數，將函數取梯度
#start: 起始位置(給array)
#learn_rate: learning rate(給自訂值)
#n_iter: 重複走的次數
#tolerance: 太小了，已接近最小值



##2. 輸入一些基本資料
#geo1: correct
#geo2: test
#讀取xlsx檔名、分析哪個方向資料
File= xlrd.open_workbook('testData.xlsx')
Sheet= File.sheet_by_index(0)

total_time=50
dt=0.002
sample_rate=1/dt
sample_num= int(sample_rate*total_time)
time=np.linspace(0, total_time, sample_num, endpoint=False)

X_direction=[0,3,6]
Y_direction=[1,4,7]
Z_direction=[2,5,8]
directionYouWantToAnaylze= X_direction

#之後分析的頻率段
freq_lowerlimit= 1
freq_upperlimit= 249

Geo1Data= Sheet.col_values(directionYouWantToAnaylze[0],23,23+sample_num)
Geo1Data= np.array(Geo1Data)
Geo2Data= Sheet.col_values(directionYouWantToAnaylze[1],23,23+sample_num)
Geo2Data= np.array(Geo2Data)

#plot 第一手資料:電壓時域圖E(t)
plt.figure(1)
plt.suptitle("E(time)")
down= min(min(Geo1Data),min(Geo2Data))
up= max(max(Geo1Data),max(Geo2Data))

plt.subplot(2,1,1)
plt.xlabel("time(s)")
plt.ylabel("voltage(V)")
plt.ylim(down,up)
plt.plot(time,Geo1Data)
plt.subplot(2,1,2)
plt.xlabel("time(s)")
plt.ylabel("voltage(V)")
plt.ylim(down,up)
plt.plot(time,Geo2Data)



##3. 傅立葉轉換(使用scipy.fft套件)
Geo1Data_fft= fft(Geo1Data)
Geo2Data_fft= fft(Geo2Data)
fd = np.linspace(0.0, int(sample_rate), int(sample_num), endpoint=False)

#要注意Nyquist效應(對稱答案有一半會不見)。基於能量不變，所以振幅全部都乘兩倍
E_Geo1Data= 2/sample_num *Geo1Data_fft[0:int(sample_num/2)]
E_Geo2Data= 2/sample_num *Geo2Data_fft[0:int(sample_num/2)]
fd=fd[0:int(sample_num/2)]

#取振幅
amp_E_Geo1Data= np.abs(E_Geo1Data)
amp_E_Geo2Data= np.abs(E_Geo2Data)

#plot 轉換後資料:電壓頻域圖E1(f)、E2(f)
plt.figure(2)
plt.suptitle("E(freq)")
down= min(min(amp_E_Geo1Data),min(amp_E_Geo2Data))
up= max(max(amp_E_Geo1Data),max(amp_E_Geo2Data))

plt.subplot(2,1,1)
plt.xlabel("freq(Hz)")
plt.ylabel("amp(V)")
plt.ylim(down,up)
plt.plot(fd, amp_E_Geo1Data)
plt.subplot(2,1,2)
plt.xlabel("freq(Hz)")
plt.ylabel("amp(V)")
plt.ylim(down,up)
plt.plot(fd, amp_E_Geo2Data)

#plot correct E1(f) v.s test E2(f)
plt.figure(3)
plt.title("E1(freq)-E2(freq)")
plt.xlabel("freq(Hz)")
plt.ylabel("amp")
plt.plot(fd, amp_E_Geo1Data-amp_E_Geo2Data)



##4. 地聲儀、類比轉換器型號
#地聲儀:GS-20DX
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

#plot 儀器響應曲線
plt.figure(4)
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



##5. 算出地表震動速度頻域圖V1(f) 
V_Geo1Data=[0]*len(fd)
for i in range(len(fd)):
    if fd[i]<freq_upperlimit and fd[i]>freq_lowerlimit:
        V_Geo1Data[i]=E_Geo1Data[i]/(const*(1j*wd[i]-z1)*(1j*wd[i]-z2)/(1j*wd[i]-p1)/(1j*wd[i]-p2))
    else:
        V_Geo1Data[i]=0

#取振幅
amp_V_Geo1Data= np.abs(V_Geo1Data)
        
#plot地表震動速度頻域圖V1(f)
plt.figure(5)
plt.title("velocity(freq)")
plt.xlabel("freq(Hz)")
plt.ylabel("amp(in/s)")
plt.xlim(0,freq_upperlimit)
plt.plot(fd, amp_V_Geo1Data)

for i in range(len(fd)):
    if amp_V_Geo1Data[i]!=0:
        start=i
        break
    
for i in range(start+1,len(fd)):
    if amp_V_Geo1Data[i]==0:
        end=i
        break

wd= wd[start:end]
amp_V_Geo1Data= amp_V_Geo1Data[start:end]
amp_E_Geo2Data= amp_E_Geo2Data[start:end]



##6. 開始做steepest descent method

#輸入learningRate、猜測次數(guessTime)、初始猜測值
learningRate=100000
guessTime=100

G=Go
k=1100
M=m
h0=0.03
R=Rt
start=[G,k,M,h0,R]

#紀錄每步的猜測(G,k,M,h0,R、z(查看目標函數是否越來越小))
recordingBox= np.zeros((guessTime,6))


z=0
for n in range(guessTime):
    dzdG=0
    dzdk=0
    dzdM=0
    dzdh0=0
    dzdR=0
    
    for i in range(len(wd)):    
        temp_commonTerm= -((amp_V_Geo1Data[i])*amp_E_Geo2Data[i]/(wd[i])**2)*((G**(-2)*k**2*M**(-2)-2*wd[i]**2*G**(-2)*k*M**(-1)+wd[i]**4*G**(-2)+4*wd[i]**2*G**(-2)*k*M**(-1)*h0**2+4*wd[i]**2*k**(1/2)*M**(-3/2)*h0*R**(-1)+wd[i]**2*G**2*M**(-2)*R**(-2)))**(-1/2)
        
        temp_G= (-2)*G**(-3)*k**2*M**(-2)+4*wd[i]**2*G**(-3)*k*M**(-1)-2*wd[i]**4*G**(-3)-8*wd[i]**2*G**(-3)*k*M**(-1)*h0**2+2*wd[i]**2*G*M**(-2)*R**(-2)
        dzdG += (amp_E_Geo2Data[i]**2/wd[i]**4)*temp_G+temp_commonTerm*temp_G
        
        temp_k= 2*G**(-2)*k*M**(-2)-2*wd[i]**2*G**(-2)*M**(-1)+4*wd[i]**2*G**(-2)*M**(-1)*h0**2+2*wd[i]**2*k**(-1/2)*M**(-3/2)*h0*R**(-1)
        dzdk += (amp_E_Geo2Data[i]**2/wd[i]**4)*temp_k+temp_commonTerm*temp_k
        
        temp_M= (-2)*G**(-2)*k**2*M**(-3)+2*wd[i]**2*G**(-2)*k*M**(-2)-4*wd[i]**2*G**(-2)*k*M**(-2)*h0**2-6*wd[i]**2*k**(1/2)*M**(-5/2)*h0*R**(-1)-2*wd[i]**2*G**2*M**(-3)*R**(-2)
        dzdM += (amp_E_Geo2Data[i]**2/wd[i]**4)*temp_M+temp_commonTerm*temp_M
        
        temp_h0= 8*wd[i]**2*G**(-2)*k*M**(-1)*h0+4*wd[i]**2*k**(1/2)*M**(-3/2)*R**(-1)
        dzdh0 += (amp_E_Geo2Data[i]**2/wd[i]**4)*temp_h0+temp_commonTerm*temp_h0
        
        temp_R= (-4)*wd[i]**2*k**(1/2)*M**(-3/2)*h0*R**(-2)-2*wd[i]**2*G**2*M**(-2)*R**(-3)
        dzdR += (amp_E_Geo2Data[i]**2/wd[i]**4)*temp_R+temp_commonTerm*temp_R
        
    G += learningRate*dzdG
    k += learningRate*dzdk
    M += learningRate*dzdM
    h0 += learningRate*dzdh0
    R += learningRate*dzdR
    start=[G,k,M,h0,R]
    
    for i in range(len(wd)):   
        z= ((amp_V_Geo1Data[i])**2)\
            -(2*amp_V_Geo1Data[i]*amp_E_Geo2Data[i]/(wd[i])**2)*(G**(-2)*k**2*M**(-2)-2*(wd[i]**2)*G**(-2)*k*M**(-1)+(wd[i])**4*G**(-2)+4*(wd[i])**2*G**(-2)*k*M**(-1)*(h0**2)+4*(wd[i])**2*k**(1/2)*M**(-3/2)*h0*R**(-1)+(wd[i])**2*G**2*M**(-2)*R**(-2))**0.5\
            +((amp_E_Geo2Data[i])**2/wd[i]**4*(G**(-2)*k**2*M**(-2)-2*wd[i]**2*G**(-2)*k*M**(-1)+wd[i]**4*G**(-2)+4*wd[i]**2*G**(-2)*k*M**(-1)*h0**2+4*wd[i]**2*k**(1/2)*M**(-3/2)*h0*R**(-1)+wd[i]**2*G**2*M**(-2)*R**(-2)))
    
    for i in range(5):
        recordingBox[n,i]=start[i]
    recordingBox[n,5]=z  
    
print('ok')

#把資料存入excel
# from openpyxl import Workbook
# wb= Workbook()
# wb.create_sheet('Data',index=0)

# row=["G","k","M","h0","R","z(目標函數)"]

# wb.save('RecordingData.xlsx')






    
        
            
    
    
    
        


    



















