# -*- coding: utf-8 -*-
"""
Created on Sat Oct 15 17:57:23 2022

@author: jefe_
"""

#Import standard packages
import numpy as np
import matplotlib.pyplot as plt

from sklearn import linear_model as lm
from sklearn import model_selection as ms
from sklearn import discriminant_analysis as da
from sklearn import metrics as mc


from scipy import io, stats

import pathlib as pl

#Import metrics
from Neural_Decoding.metrics import get_rho, get_R2

# Import preprocessing functions
from Neural_Decoding.preprocessing_funcs import bin_output, bin_spikes

#Import decoder functions
from Neural_Decoding.decoders import WienerCascadeDecoder
from Neural_Decoding.decoders import WienerFilterDecoder
from Neural_Decoding.decoders import GRUDecoder
from Neural_Decoding.decoders import XGBoostDecoder
from Neural_Decoding.decoders import SVRDecoder

def my_zscore(a, axis=0, ddof=0, nan_policy='omit'):
    a_mu = np.mean(a, axis=axis)
    a_sig = np.std(a, axis=axis, ddof=ddof)
    a_z = stats.zscore(a, axis=axis, ddof=ddof, nan_policy=nan_policy)
    return a_mu, a_sig, a_z
    
folder=pl.Path(
    r'Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch2_ephys\MC\GAD18\211205_C')
neu_file = r'GADi18_SpkTms+vels.mat'
beh_file = r'BehaveSignals2021-12-05T14_28_34.mat'

neu_mat = folder.as_posix() + r'\\ephys_F\\'+ neu_file
beh_mat = folder.as_posix() + r'\\Behaviour\\'+ beh_file

neu_data = io.loadmat(neu_mat)
beh_data = io.loadmat(beh_mat)

# Neural data
variable_info = io.whosmat(neu_mat)
variable_names = []
for var_ in variable_info:
    variable_names.append(var_[0])
print(variable_names)
# Behaviour data
variable_info = io.whosmat(beh_mat)
variable_names = []
for var_ in variable_info:
    variable_names.append(var_[0])
print(variable_names)

fr = np.squeeze(neu_data['fr'])
fs = 3e4
spike_times = neu_data['spike_times']

nose = beh_data['nose']
lw = beh_data['lw']
rw = beh_data['rw']
time_axis = np.linspace(0, nose.shape[0], nose.shape[0], endpoint=True)/fr

spike_times = np.squeeze(spike_times)
for cst, sp in enumerate(spike_times):
    spike_times[cst] = np.squeeze(sp)

dt = 0.05
t_start = 0
t_end = nose.shape[0]/fr

neural_data = bin_spikes(spike_times, dt, t_start, t_end)
lw_binned = bin_output(lw, time_axis, dt, t_start, t_end)
rw_binned = bin_output(rw, time_axis, dt, t_start, t_end)
nose_binned = bin_output(nose, time_axis, dt, t_start, t_end)

# Format for Wiener Filter, Wiener Cascade, XGBoost, and Dense Neural Network
#Put in "flat" format, so each "neuron / time" is a single feature
neu_mu, neu_sg, X = my_zscore(neural_data, axis=0, ddof=1)
lw_mu, lw_sg, lw_z = my_zscore(lw_binned, ddof=1)
rw_mu, rw_sg, rw_z = my_zscore(rw_binned, ddof=1)
ns_mu, ns_sg, ns_z = my_zscore(nose_binned, ddof=1)

train_pc = 0.7
train = int(np.round(train_pc*X.shape[0]))
time_CV = ms.TimeSeriesSplit(n_splits=10)
rdge_ns = lm.RidgeCV(cv=time_CV, alphas=np.logspace(-3,3,num=7))
rdge_ns.fit(X[:train,:], ns_z[:train])

score, _, pval = ms.permutation_test_score(rdge_ns, X, ns_z, n_permutations=500, 
                                           cv=time_CV)

r2_wf = np.zeros((10, 3))
for cit, tt in enumerate(time_CV.split(nose_binned)):
    train, test = tt
    wf_nose = lm.LinearRegression()
    wf_lw = lm.LinearRegression()
    wf_rw = lm.LinearRegression()
    
    wf_nose.fit(X[train,:], ns_z[train])
    wf_lw.fit(X[train,:], lw_z[train])
    wf_rw.fit(X[train,:], rw_z[train])
    
    r2_wf[cit, 0] = get_R2(ns_z[test], wf_nose.predict(X[test,:]))
    r2_wf[cit, 1] = get_R2(lw_z[test], wf_lw.predict(X[test,:]))
    r2_wf[cit, 2] = get_R2(rw_z[test], wf_rw.predict(X[test,:]))
    
    
    
    



"""
# X contains the PSTH per trial of all units concatenated as if in continuous 
# time
X = data['X']
y = data['y']

Xn = stats.zscore(X, axis=1, nan_policy='omit')
Xn[np.isnan(Xn)] = 0
y_mu = y.mean()
yc = y - y_mu
y_svm = stats.zscore(y)

X_train, X_test, y_train, y_test = ms.train_test_split(Xn, yc, test_size=0.15,
                                                       shuffle=False)
_, _, y_svm_train, y_svm_test = ms.train_test_split(Xn, y_svm, test_size=0.15,
                                                    shuffle=False)

cv = 5
shuf_splitter = ms.ShuffleSplit(cv, test_size=0.2)
alphas = np.logspace(-4, 2)
scores = np.zeros((alphas.shape[0], cv))
for a, alpha in enumerate(alphas):
    ridge_mdl = lm.TweedieRegressor(power=0, alpha=alpha, link='auto')
    for t, idxs in enumerate(shuf_splitter.split(X_train)):
        train, vali = idxs
        ridge_mdl.fit(X_train[train,:], y_train[train])
        y_pred = ridge_mdl.predict(X_train[vali,:])
        #scores[a, t] = ridge_mdl.score(X_train[vali,:], y_train[vali])
        scores[a, t] = mc.r2_score(y_train[vali], y_pred=y_pred)
        #scores[a, t] = get_R2(y_train[vali], y_pred)
        

plt.figure()
plt.semilogx(alphas, scores.mean(axis=1))
plt.xlabel("Constraint weight")
plt.ylabel("R2 score in cross-validation")
plt.figure(figsize=(10,10))
optimal_C = np.argmax(scores.mean(axis=1))
ridge_mdl = lm.TweedieRegressor(power=0, alpha=alphas[optimal_C], link='auto')
ridge_mdl.fit(X_train[train,:], y_train[train])
r2_test = mc.r2_score(y_test, ridge_mdl.predict(X_test))
print("R2 on test: {}".format(r2_test))
coefficients = ridge_mdl.coef_.reshape(59,140)
plt.xlabel("Time [ms]")
plt.ylabel("Units")
plt.imshow(coefficients, cmap='coolwarm', vmin=-np.max(coefficients), 
           vmax=np.max(coefficients), extent=(-300,400,0,58), aspect='auto')
plt.colorbar(orientation = 'horizontal', shrink = .6, 
               label="Coefficient value.")
plt.tight_layout()
plt.show()
plt.figure(figsize=(10,5))
plt.plot(np.arange(0,X.shape[0]), y, ridge_mdl.predict(Xn)+y_mu)
plt.xlabel("Trial")
plt.ylabel("Roller speed [cm/s]")
plt.legend(["{}".format(lbl) for lbl in ("Ground truth", "Predicted {}".format(r2_test))])
plt.show()
"""