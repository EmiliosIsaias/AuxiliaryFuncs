# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 11:20:47 2022

@author: Emilio Isaias-Camacho
"""

import pickle
from sklearn.model_selection import train_test_split

from scipy import stats
# from scipy.signal import convolve, windows # Plan to smooth the histograms

import pathlib as pl
import numpy as np
import matplotlib.pyplot as plt

from Neural_Decoding.decoders import KalmanFilterDecoder
from Neural_Decoding import metrics

from sklearn import model_selection

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

bin_size = 0.05

# Loading data
folder = pl.Path(
    r'Z:\Emilio\SuperiorColliculusExperiments\Roller' + 
    r'\Batch2_ephys\MC\GAD18\211205_C\ephys_F')
folder = folder.as_posix()

file_name = folder + r'/GADi18_SpkTms+vels.pickle'

with open(file_name, 'rb') as f:
    neural_data, vels_binned, vel_times, fr = pickle.load(f)

# ======================= Preprocessing =======================================
# For the Kalman filter, we use the position, velocity, and acceleration as
# outputs. Ultimately, we are only concerned with the goodness of fit of 
# velocity (for this dataset) but using them all as covariates helps 
# performance
#We will now determine position
pos_binned = np.zeros(vels_binned.shape) #Initialize 
pos_binned[0] = vel_times[0] #Assume starting position is at [0,0]
#Loop through time bins and determine positions based on the velocities
for i in range(pos_binned.shape[0]-1): 
    pos_binned[i+1] = pos_binned[i]+vels_binned[i]*bin_size
    
#We will now determine acceleration    
temp=np.diff(vels_binned,axis=0) 
#Assume acceleration at last time point is same as 2nd to last
acc_binned=np.concatenate((temp,temp[-1:,:]),axis=0) 
    
#The final output covariates include position, velocity, and acceleration
y_kf_o=np.concatenate((pos_binned,vels_binned,acc_binned),axis=1)
    
    
#The covariate is simply the matrix of firing rates for all neurons over time
X_kf_o = neural_data
# conv_kernel = windows.gaussian(3, std)
# for ccl in range(0, X_kf.shape[1]):
#     X_kf[:,ccl] = convolve

num_examples = X_kf_o.shape[0]

# Kalman filter history lag = -1 means 1 bin before the output
lags = -10
Nlags = np.abs(lags)
Nc = 20
Cs = np.logspace(-1, 3, num=Nc)
tss_it = model_selection.TimeSeriesSplit(n_splits=5)
r2_lag = np.zeros(Nlags)
for l, lag in enumerate(range(lags,0)):
    print("L: ", lag)
    X_kf = X_kf_o
    y_kf = y_kf_o
    
    
    #Re-align data to take lag into account
    if lag<0:
        y_kf=y_kf[-lag:,:]
        X_kf=X_kf[0:num_examples+lag,:]
    if lag>0:
        y_kf=y_kf[0:num_examples-lag,:]
        X_kf=X_kf[lag:num_examples,:]
    
    X_kf = stats.zscore(X_kf, axis=0)
    vel_mean = y_kf.mean(axis=0)
    y_kf -= vel_mean

    # ================================ Splitting =============================
    
    r2_mat = np.zeros((Nc, np.abs(lags)))
    for y, C in enumerate(Cs):
        print("C: ", C)
        r2s = []
        for x, idxs in enumerate(tss_it.split(X_kf)):
            train, test = idxs
            # Model definition (maybe not necessary to re-instanciate)
            kf_model = KalmanFilterDecoder(C=C)
            kf_model.fit(X_kf[train,:], y_kf[train,:])
            y_test = kf_model.predict(X_kf[test,:], y_kf[test,:])
            r2s.append(metrics.get_R2(y_kf[test,:], y_test))
            #print('L: {}, C: {}, R2: {}'.format(lag, C, r2s[1,:].mean(axis=0)))
        r2s = np.array(r2s)
        r2_mat[y,l]=r2s[1,:].mean(axis=0)
        
    r2_mean = r2_mat
    opt_C = np.argmax(r2_mean)
    print("Max R²: {}, optimal C: {}".format(r2_mean[opt_C], Cs[opt_C]))
    r2_lag[l] = r2_mean[opt_C]
plt.figure(figsize=(10,5))
plt.plot(np.arange(lags,0, 1), r2_lag)


#r2_mean = r2_mat.mean(axis=1)
#opt_C = np.argmax(r2_mean)
#kf_model = KalmanFilterDecoder(C=Cs[opt_c])
#kf_model.fit(X_kf[train,:], y_)
#y_test = kf_model.predict(X_kf[test,:], y_kf[test,:])
#plt.figure(figsize=(10,5))
#plt.xlabel('Constraint weight')
#plt.ylabel('R2')
#plt.semilogx(Cs, r2_mat.mean(axis=0))

