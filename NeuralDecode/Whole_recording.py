# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 11:20:47 2022

@author: Emilio Isaias-Camacho
"""

import pickle

from scipy import io
from scipy import stats

import pathlib as pl
import numpy as np
import matplotlib.pyplot as plt

from Neural_Decoding import preprocessing_funcs
from Neural_Decoding.decoders import KalmanFilterDecoder
from Neural_Decoding import metrics

from sklearn import linear_model
from sklearn import preprocessing

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


# Preprocessing ##
# Kalman filter history lag = -1 means 1 bin before the output
lag = -4

# For the Kalman filter, we use the position, velocity, and acceleration as outputs
# Ultimately, we are only concerned with the goodness of fit of velocity (for this dataset)
# But using them all as covariates helps performance

#We will now determine position
pos_binned = np.zeros(vels_binned.shape) #Initialize 
pos_binned[0] = vel_times[0] #Assume starting position is at [0,0]
#Loop through time bins and determine positions based on the velocities
for i in range(pos_binned.shape[0]-1): 
    pos_binned[i+1] = pos_binned[i]+vels_binned[i]*.05 #Note that .05 is the length of the time bin

#We will now determine acceleration    
temp=np.diff(vels_binned,axis=0) #The acceleration is the difference in velocities across time bins 
acc_binned=np.concatenate((temp,temp[-1:,:]),axis=0) #Assume acceleration at last time point is same as 2nd to last

#The final output covariates include position, velocity, and acceleration
y_kf=np.concatenate((pos_binned,vels_binned,acc_binned),axis=1)


