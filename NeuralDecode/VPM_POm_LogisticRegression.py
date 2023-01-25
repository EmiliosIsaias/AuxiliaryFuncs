# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 14:13:30 2022

@author: Emilio Isaias-Camacho
"""

from sklearn import linear_model, model_selection, metrics

import matplotlib.pyplot as plt
import matplotlib.figure as Fig
import numpy as np
import scipy.stats as sts
from scipy import io
import pathlib as pl

data_path = pl.Path(r"E:\Database")
file_name = "vpm_pom_cluster_data_for_Emilio.mat"
file_path = data_path.as_posix() + pl.os.sep + file_name
data = io.loadmat(file_path)


variable_info = io.whosmat(file_path)
variable_names = []
for var_ in variable_info:
    variable_names.append(var_[0])
print(variable_names)

bC = data['Bins_PSTH']
puffH = data['Puff_PSTH']
touchH = data['Touch_PSTH']
Yv = np.squeeze(data['is_vpm'])
Yp = np.squeeze(data['is_pom'])


"""
# Puff and touch
X = np.concatenate((sts.zscore(puffH, axis=1),
                    sts.zscore(touchH, axis=1)), axis=1)
legStrs = (("Puff", "Touch"), "Puff", "Touch")
title = ("Puff + Touch","Puff","Touch")

# Puff
X = sts.zscore(puffH, axis=1)
legStr = "Puff"
title = "Puff"
"""
"""
# Touch
X = sts.zscore(touchH, axis=1)
legStr = "Touch"
title = "Touch"
"""
legStrs = (("Puff", "Touch"), "Puff", "Touch")
title = ("Puff + Touch","Puff","Touch")

def logReg(X, Yv, title_string, leg_string, test_size=0.3):
    X_train, X_test, y_train, y_test = model_selection.train_test_split(
        X, Yv, test_size=test_size)
    
    log_reg_model = linear_model.LogisticRegressionCV(
        Cs=np.logspace(-3, 3, num=50), penalty='l2', solver='lbfgs', 
        n_jobs=-1)
    
    #log_reg_model.fit(X_train, y_train)
    #print("Chosen C:{}".format(log_reg_model.C_))
    tot_score, perm_scores, p = model_selection.permutation_test_score(
        log_reg_model, X_train, y_train, n_permutations=256, verbose=True, n_jobs=-1,
        scoring='accuracy')
    
    log_reg_model.fit(X_train, y_train)
    
    print("Model accuracy: {} | Mean score: {} | P: {} | Precision: {}".format(
        metrics.accuracy_score(y_test, log_reg_model.predict(X_test)), 
        perm_scores.mean(), p, 
        metrics.precision_score(y_test, log_reg_model.predict(X_test))))
    
    # Plotting and saving results
    psth_tx = np.arange(-24.5, 75)
    coefs = log_reg_model.coef_
    if coefs.shape[1] > 100:
        coefs = coefs.reshape([2,100])
    plt.figure()
    plt.plot(psth_tx, coefs.transpose(), label=leg_string)
    plt.legend()
    plt.xlabel("Time [ms]"); plt.ylabel("Coefficient magnitude")
    
    
    fig, ax = plt.subplots()
    
    ax.hist(perm_scores, bins=20, density=True)
    ax.axvline(tot_score, ls="--", color="r")
    score_label = f"Score on original\ndata: {perm_scores:.2f}\n(p-value: {p:.3f})"
    ax.text(0.7, 10, score_label, fontsize=12)
    ax.set_xlabel("Accuracy score")
    _ = ax.set_ylabel("Probability")

for cx, X in enumerate((np.concatenate((sts.zscore(puffH, axis=1),
                    sts.zscore(touchH, axis=1)), axis=1),
                        sts.zscore(puffH, axis=1),
                        sts.zscore(touchH, axis=1))):
    print("Using {}".format(title[cx]))
    logReg(X, Yv, title[cx], legStrs[cx])

        