# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# 畫某個函數的圖
x= np.linspace(-100,100,100)
y= np.linspace(-100,100,100)
X, Y= np.meshgrid(x, y)
Z= 3*(X-5)**2+2*(Y-20)**2


plt.figure(1)
C= plt.contour(X,Y,Z, 15)
plt.clabel(C, inline= True, fontsize= 10)

fig= plt.figure(2)
ax= Axes3D(fig)
ax.plot_surface(X,Y,Z, rstride=5, cstride=5, cmap='rainbow')
plt.show()


# 最陡下降法(SDM)
learningRate= -0.01
guessTime= 2000
x0= 100
y0= 100
t0= 0


start=[x0,y0,t0]
recordingBox= np.zeros((guessTime,7))

f=0
for n in range(guessTime):
    dfdx=0
    dfdy=0
    dfdt=0
    
    f= 3*((start[0]-5)**2)+2*(start[1]-20)**2+(20*start[2])
    dfdx=6*start[0]-30
    dfdy=4*start[1]-80
    dfdt=20
    
    for i in range(0,3):
        recordingBox[n,i]=start[i]
        
    recordingBox[n,3]=f
    recordingBox[n,4]=dfdx*learningRate
    recordingBox[n,5]=dfdy*learningRate
    recordingBox[n,6]=dfdt*learningRate
    
    start[0] += learningRate*dfdx
    start[1] += learningRate*dfdy
    start[2] += learningRate*dfdt  
print('ok')




# 記錄到excel
import xlwings as xw
wb= xw.Book("SDM紀錄檔案.xlsx")
sheet= wb.sheets["範例1"]
sheet.clear()

title= ["X","Y","T","Z","下一步(X)","下一步(Y)","下一步(T)"]
sheet.range("A1").value= title
sheet.range("A2").value= recordingBox

# wb.save("SDM紀錄檔案")













