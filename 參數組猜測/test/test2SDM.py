# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 18:47:29 2021

@author: r8891
"""

import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import xlwings as xw

# 畫某個函數的圖
x= np.linspace(-100,100,100)
y= np.linspace(-100,100,100)
X, Y= np.meshgrid(x, y)
Z= 6000*np.sin(np.pi*X/24)+3*(X-36)**2+2*(Y-20)**2

plt.figure(1)
C= plt.contour(X,Y,Z, 20)
plt.clabel(C, inline= True, fontsize= 10)

fig= plt.figure(2)
ax= Axes3D(fig)
ax.plot_surface(X,Y,Z, rstride=5, cstride=5, cmap='rainbow')
plt.show()

# 最陡下降法(SDM)
learningRate= -0.01
guessTime= 2500
x0= 100
y0= 100
t0= 0

##1. 不加momentum

# start=[x0,y0,t0]
# recordingBox= np.zeros((guessTime,7))

# f=0
# for n in range(guessTime):
#     dfdx=0
#     dfdy=0
#     dfdt=0
    
#     f= 6000*np.sin(np.pi*start[0]/24)+3*(start[0]-36)**2+2*(start[1]-20)**2+20*start[2]
#     dfdx=250*np.pi*np.cos(np.pi*start[0]/24)+6*start[0]-216
#     dfdy=4*start[1]-80
#     dfdt=20
    
#     for i in range(0,3):
#         recordingBox[n,i]=start[i]
        
#     recordingBox[n,3]= f
#     recordingBox[n,4]= dfdx*learningRate
#     recordingBox[n,5]= dfdy*learningRate
#     recordingBox[n,6]= dfdt*learningRate
    
#     start[0] += learningRate*dfdx
#     start[1] += learningRate*dfdy
#     start[2] += learningRate*dfdt  
# print('ok')

# wb= xw.Book("SDM紀錄檔案.xlsx")
# sheet= wb.sheets["範例2(不加momentum)"]
# sheet.clear()
# title= ["X","Y","T","Z","下一步(X)","下一步(Y)","下一步(T)"]
# sheet.range("A1").value= title
# sheet.range("A2").value= recordingBox


##2. 加momentum

start=[x0,y0,t0]
recordingBox= np.zeros((guessTime,10))
m_const= 0.99


for n in range(guessTime):
    f=0
    dfdx=0
    dfdy=0
    dfdt=0
    
    f= 6000*np.sin(np.pi*start[0]/24)+3*(start[0]-36)**2+2*(start[1]-20)**2+20*start[2]
    dfdx=250*np.pi*np.cos(np.pi*start[0]/24)+6*start[0]-216
    dfdy=4*start[1]-80
    dfdt=20
    
    for i in range(0,3):
        recordingBox[n,i]=start[i]
        
    recordingBox[n,3]= f
    
    recordingBox[n,4]= dfdx*learningRate
    recordingBox[n,5]= dfdy*learningRate
    recordingBox[n,6]= dfdt*learningRate
    
    
    for i in range(7,10):  
        if n==0:
            recordingBox[n,i]= (1-m_const)*recordingBox[n,i-3]
        else:
            recordingBox[n,i]=(1-m_const)*recordingBox[n,i-3]+m_const*recordingBox[n-1,i]
            
    
    start[0] += recordingBox[n,7]
    start[1] += recordingBox[n,8]
    start[2] += recordingBox[n,9] 
    
print('ok')

wb= xw.Book("SDM紀錄檔案.xlsx")
sheet= wb.sheets["範例2(加momentum)"]
sheet.clear()
title= ["X","Y","T","Z","原本下一步(X)","原本下一步(Y)","原本下一步(T)","加權後的下一步(X)","加權後的下一步(Y)","加權後的下一步(T)"]
sheet.range("A1").value= title
sheet.range("A2").value= recordingBox




