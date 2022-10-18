# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 10:54:07 2021

@author: Emilio Isaias-Camacho

"""
import math
from scipy import io
import numpy as np
import quantities as pq
import neo
import elephant
#import viziphant
import pathlib as pl

#Windows
folder = pl.Path(
    r'Z:\Emilio\SuperiorColliculusExperiments\Roller' +
    r'\Batch2_ephys\MC\GAD18\211205_C\ephys_F')

"""
 #Linux
folder = pl.Path(
    r'/mnt/sds-hd/sd19b001/Emilio/SuperiorColliculusExperiments/'+
    r'Roller/Batch2_ephys/MC/GAD18/211205_C/ephys_F')
"""
folder = folder.as_posix()
file_name = folder + r'/GADi18_SpkTms+vels.mat'


data = io.loadmat(file_name)

variable_info = io.whosmat(file_name)
variable_names = []
for var_ in variable_info:
    variable_names.append(var_[0])
print(variable_names)

spike_times = data['spike_times']
uID = data['gclID']

spike_times = np.squeeze(spike_times)
spiketrains = []
for cst, spkTms in enumerate(spike_times):
    spiketrains.append(neo.core.SpikeTrain(
        np.squeeze(spike_times[cst])*pq.s,
        t_stop=1320*pq.s))
print("Spike train formatted. Starting SPADE:")

n_spikes = [len(spiketrain) for spiketrain in spiketrains]

bin_size = 5e-3 * pq.s

firing_rate_two_std = (np.mean(n_spikes) + 2 * np.std(n_spikes))/spiketrains[0].t_stop
firing_probability_two_std = (firing_rate_two_std * bin_size).simplified.item()

n_bins = int(spiketrains[0].t_stop/bin_size)

patterns = []
for min_spikes in range(2, 6):
    min_occ = math.ceil(firing_probability_two_std**min_spikes * n_bins)
    print(f'{min_spikes=}', f'{min_occ=}')
    spade_output = elephant.spade.spade(
        spiketrains=spiketrains,
        bin_size=bin_size,
        winlen=1,                   # 1 bin means synchronous patterns
        min_spikes=min_spikes,
        max_spikes=None,  #if min_spikes < 4 else None,
        min_occ=min_occ,
        min_neu=2,
        n_surrogates=1000,
        dither=15e-3*pq.s,
        psr_param=[0,0,3],
        spectrum='3d#')
    patterns.append(spade_output['patterns'])
    print(f'Patterns found for {min_spikes} spikes: {len(spade_output["patterns"])}')
