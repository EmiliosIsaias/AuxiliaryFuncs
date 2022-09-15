# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 10:54:07 2021

@author: Emilio Isaias-Camacho

#import os
import numpy as np
import quantities as pq
import matplotlib.pyplot as plt
import neo
from viziphant.rasterplot import rasterplot
from viziphant.events import add_event
from elephant.spade import spade

def get_channel_units_for_block(block):
    all_channels_units = []
    for i, train in enumerate(block.segments[0].spiketrains):
        channel_id = train.annotations["channel_id"]
        unit_id = train.annotations["unit_id"]
        all_channels_units.append((channel_id, unit_id))
    return all_channels_units

# Relative path to data (change to where you saved them)
path = './data/'

# load data
purple = np.load(path + 'purple.npy', allow_pickle=True).item()
red = np.load(path + 'red.npy', allow_pickle=True).item()
blue = np.load(path + 'blue.npy', allow_pickle=True).item()
black = np.load(path + 'black.npy', allow_pickle=True).item()
orange = np.load(path + 'orange.npy', allow_pickle=True).item()
green = np.load(path + 'green.npy', allow_pickle=True).item()

colors = ['purple', 'red', 'blue', 'black', 'orange', 'green']



#for block_idx, block in enumerate([red, blue, black]):
for block_idx, block in enumerate([red]):
    all_sts = []
    all_channel_units = get_channel_units_for_block(block)
    for seg_id, segment in enumerate(block.segments):
        for i,(channel,units) in enumerate(all_channel_units):
          dict_list = [{'channel_id':channel},
                        {'unit_id': units}]
          sts = segment.filter(targdict=dict_list, objects=neo.SpikeTrain)[0]
          if seg_id == 0:
            all_sts.append(np.array(sts)[(0 <= sts) & (sts <= 0.5)])
          else:
            sts_to_append = np.array(sts)[(0 <= sts) & (sts <= 0.5)]+ seg_id*0.5
            all_sts[i] = np.hstack((all_sts[i], sts_to_append))
            neo_all_sts = [neo.SpikeTrain(all_spikes, t_start=0, 
                              t_stop=len(block.segments)*0.5, units=pq.s) 
              for all_spikes in all_sts]

neo_all_sts = [neo.SpikeTrain(all_spikes, t_start=0, 
                              t_stop=len(block.segments)*0.5, units=pq.s) 
            for all_spikes in all_sts]
spade_result = spade(neo_all_sts, 0.005 * pq.s, winlen=1, # synchronous events
                     min_spikes=5, 
                     min_occ=5, min_neu=5, n_surr=10)
print(len(spade_result['patterns']))
spade_result
"""
from scipy import io
import numpy as np
import quantities as pq
import neo
import elephant
import viziphant
import pathlib as pl

folder = pl.Path(
    r'Z:\Emilio\SuperiorColliculusExperiments\Roller' + 
    r'\Batch2_ephys\MC\GAD18\211205_C\ephys_F')
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
print("Spike train formated")

spade_output = elephant.spade.spade(
    spiketrains=spiketrains,
    bin_size=1*pq.s,
    winlen=1,                   # 1 bin means synchronous patterns
    min_spikes=3,
    n_surrogates=100,
    dither=5*pq.ms,
    psr_param=[0, 0, 0])

patterns = spade_output['patterns']

viziphant.patterns.plot_patterns(spike_times, patterns)



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