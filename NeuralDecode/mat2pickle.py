# -*- coding: utf-8 -*-
"""
Created on Fri Sep  2 11:51:32 2022

@author: Emilio Isaias-Camacho
"""

import numpy as np
from scipy import io
import pathlib as pl
import pickle

from Neural_Decoding.preprocessing_funcs import bin_spikes
from Neural_Decoding.preprocessing_funcs import bin_output


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

spike_times = data['spike_times']
vels = data['vels']
fr = data['fr']
vel_times = np.arange(0, vels.shape[1])/fr[0][0]

bin_size = 0.05 # in seconds
t_i = vel_times[0] # start of the experiment
t_e = vel_times[-1] # end of the experiment
downsampling_factor = 1 # no downsampling

spike_times = np.squeeze(spike_times)
for ccl in range(spike_times.shape[0]):
    spike_times[ccl]=np.squeeze(spike_times[ccl])

neural_data = bin_spikes(spike_times, bin_size, t_i, t_e)
vels_binned = bin_output(vels.transpose(), vel_times, bin_size, t_i, t_e)

pickle_name = file_base + '.pickle'

with open(pickle_name, 'wb') as f:
    pickle.dump([neural_data, vels_binned, vel_times, fr], f)

print('Done!')