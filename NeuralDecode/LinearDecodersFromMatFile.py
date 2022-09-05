# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 17:29:44 2022

@author: Emilio Isaias-Camacho
"""
#Import standard packages
import numpy as np
import matplotlib.pyplot as plt

from scipy import io
from scipy import stats

import pathlib as pl

from Neural_Decoding.preprocessing_funcs import bin_spikes
from Neural_Decoding.preprocessing_funcs import bin_output


#Import function to get the covariate matrix that includes spike history from previous bins
from Neural_Decoding.preprocessing_funcs import get_spikes_with_history

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
file_name = r'/GADi18_SpkTms+vels'

file_base = folder + file_name
mat_name = file_base + '.mat'

data=io.loadmat(mat_name)
variable_info = io.whosmat(mat_name)
variable_names = []
for var_ in variable_info:
    variable_names.append(var_[0])
print(variable_names)

