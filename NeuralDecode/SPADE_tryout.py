# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 10:54:07 2021

@author: Emilio Isaias-Camacho

"""
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
print("Spike train formated. Starting SPADE:")

spade_output = elephant.spade.spade(
    spiketrains=spiketrains,
    bin_size=5*pq.ms,
    winlen=1,                   # 1 bin means synchronous patterns
    min_spikes=3,
    n_surrogates=1,
    dither=3*pq.ms,
    psr_param=[0, 0, 3])

patterns = spade_output['patterns']

#viziphant.patterns.plot_patterns(spike_times, patterns)



"""
spiketrains = elephant.spike_train_generation.compound_poisson_process(
    rate=5*pq.Hz,  A=[0]+[0.98]+[0]*8+[0.02], t_stop=10*pq.s)
print(spiketrains)
for i in range(90):
    spiketrains.append(elephant.spike_train_generation.homogeneous_poisson_process(
        rate=5*pq.Hz, t_stop=10*pq.s))

patterns = elephant.spade.spade(
    spiketrains=spiketrains, binsize=1*pq.ms, winlen=1, min_spikes=3,
    n_surr=100,dither=5*pq.ms,
    psr_param=[0,0,0],
    output_format='patterns')['patterns']
    
viziphant.patterns.plot_patterns(spiketrains, patterns)
"""