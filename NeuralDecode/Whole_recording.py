# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 11:20:47 2022

@author: Emilio Isaias-Camacho
"""

import pickle

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

# Kalman filter history lag = -1 means 1 bin before the output
lag = -4

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
y_kf=np.concatenate((pos_binned,vels_binned,acc_binned),axis=1)

#The covariate is simply the matrix of firing rates for all neurons over time
X_kf = neural_data

# conv_kernel = windows.gaussian(3, std)
# for ccl in range(0, X_kf.shape[1]):
#     X_kf[:,ccl] = convolve

num_examples = X_kf.shape[0]

#Re-align data to take lag into account
if lag<0:
    y_kf=y_kf[-lag:,:]
    X_kf=X_kf[0:num_examples+lag,:]
if lag>0:
    y_kf=y_kf[0:num_examples-lag,:]
    X_kf=X_kf[lag:num_examples,:]

X_kf = stats.zscore(X_kf, axis=0)
y_kf =- y_kf.mean(axis=0)

# ===================================== Splitting =============================
tss_it = model_selection.TimeSeriesSplit(n_splits=5)

r2_mean = []

for train, test in tss_it.split(X_kf):
    # Model definition (maybe not necessary to re-instanciate)
    kf_model = KalmanFilterDecoder(C=1)
    kf_model.fit(X_kf[train,:], y_kf[train,:])
    y_test = kf_model.predict(X_kf[test,:], y_kf[test,:])
    r2_mean.append(metrics.get_R2(y_kf[test,:], y_test))



    

