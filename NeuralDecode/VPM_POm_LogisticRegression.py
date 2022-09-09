# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 14:13:30 2022

@author: Emilio Isaias-Camacho
"""

from sklearn import linear_model, model_selection

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as sts
from scipy import io
import pathlib as pl

data_path = pl.Path(r"E:\Database\vpm_pom_cluster_data_for_Emilio.mat")
data = io.loadmat(data_path.as_posix())

variable_info = io.whosmat(data_path.as_posix())
variable_names = []
for var_ in variable_info:
    variable_names.append(var_[0])
print(variable_names)

bC = data['Bins_PSTH']
puffH = data['Puff_PSTH']
touchH = data['Touch_PSTH']
Yv = np.squeeze(data['is_vpm'])
Yp = np.squeeze(data['is_pom'])


X = np.concatenate((sts.zscore(puffH, axis=1),
                    sts.zscore(touchH, axis=1)), axis=1)

X_train, X_test, y_train, y_test = model_selection.train_test_split(
    X, Yv, 
    test_size=0.15, 
    random_state=71)

log_reg_model = linear_model.LogisticRegressionCV(
    Cs=np.logspace(-3, 3, num=20), penalty='l2', solver='lbfgs')

log_reg_model.fit(X_train, y_train)

tot_score, perm_scores, p = model_selection.permutation_test_score(
    log_reg_model, X, Yv)