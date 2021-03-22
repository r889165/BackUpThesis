# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import xlrd
import xlwt
import numpy as np

sample_num= 25000
total_time= 5

file= xlrd.open_workbook('correctData.xlsx')
table= file.sheet_by_index(0)
print("總行數:"+str(table.nrows))
print("總列數:"+str(table.ncols))

data= table.col_values(0,0,sample_num)
data= np.array(data)

#使用原始資料加上一個隨機值(標準常態分布，平均:0、標準差:0.0001)
mean= 0
sigma= 0.0001
randomNum= mean+sigma* np.random.randn(len(data))

newData= data+randomNum
newData= newData.tolist()


#用xlwt將創造的資料傳回excel
testData= xlwt.Workbook()
table2= testData.add_sheet('data')

j=0
for i in newData:
    table2.write(j,0,i)
    j=j+1
    
testData.save('testData.xls')

