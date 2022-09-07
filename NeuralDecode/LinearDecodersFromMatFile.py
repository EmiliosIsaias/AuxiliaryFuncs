# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 17:29:44 2022

@author: Emilio Isaias-Camacho
"""
#Import standard packages
import numpy as np
import matplotlib.pyplot as plt

from sklearn import linear_model as lm
from sklearn import model_selection as ms
from sklearn import discriminant_analysis as da
from sklearn import metrics as mc


from scipy import io
from scipy import stats

import pathlib as pl

#Import metrics
from Neural_Decoding.metrics import get_R2
from Neural_Decoding.metrics import get_rho

#Import decoder functions
from Neural_Decoding.decoders import WienerCascadeDecoder
from Neural_Decoding.decoders import WienerFilterDecoder
from Neural_Decoding.decoders import DenseNNDecoder
from Neural_Decoding.decoders import SimpleRNNDecoder
from Neural_Decoding.decoders import GRUDecoder
from Neural_Decoding.decoders import LSTMDecoder
from Neural_Decoding.decoders import XGBoostDecoder
from Neural_Decoding.decoders import SVRDecoder


folder=pl.Path(
    r'Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch2_ephys\MC\GAD18\211205_C\ephys_F') #ENTER THE FOLDER THAT YOUR DATA IS IN
folder = folder.as_posix()
file_name = r'/GADi18_PSTH_trial x neuron-time + mvpt'

file_base = folder + file_name
mat_name = file_base + '.mat'

data=io.loadmat(mat_name)
variable_info = io.whosmat(mat_name)
variable_names = []
for var_ in variable_info:
    variable_names.append(var_[0])
print(variable_names)

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