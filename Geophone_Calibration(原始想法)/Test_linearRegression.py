# -*- coding: utf-8 -*-
"""
Created on Sun Mar 21 13:41:39 2021

@author: r8891
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import linear_model
from sklearn.preprocessing import PolynomialFeatures
from sklearn.pipeline import make_pipeline
from sklearn.datasets import make_regression
from sklearn.model_selection import train_test_split

X,y=make_regression(n_samples=100, n_features=1, noise=50)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=0)

regr=linear_model.LinearRegression()
regr.fit(X_train, y_train)

plt.scatter(X_train, y_train, color='black')
plt.plot(X_train, regr.predict(X_train),color='blue',linewidth=1)
plt.show()

w_0=regr.intercept_
w_1=regr.coef_
R= regr.score(X_train,y_train)