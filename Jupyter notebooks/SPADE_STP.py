# -*- coding: utf-8 -*-
"""
Created on Mon May 13 14:38:34 2024

@author: neuro
"""

import numpy as np
import quantities as pq
import neo
import elephant
import viziphant
from scipy import io
import pathlib as pl
import matplotlib.pyplot as plt
import math

"""
from mpi4py import MPI
assert MPI.COMM_WORLD.Get_size() > 1
rank = MPI.COMM_WORLD.Get_rank()
if rank == 0:
    1/0
    MPI.COMM_WORLD.send(None, dest=1, tag=42)
elif rank == 1:
    MPI.COMM_WORLD.recv(source=0, tag=42)
"""

rsp_mat_path = pl.Path(r"Z:\Jesus\Jittering\FULLYCurated\14_190716-2_jittering_3725_1620_1500 yes VPM good\KS2_newChanMap\ASSET.mat")
file_data = io.loadmat( rsp_mat_path.as_posix() )
file_data['rel_spike_times_4_ASSET'].shape

spike_times_c1 = np.squeeze( file_data['rel_spike_times_4_ASSET'] )
Ncl = spike_times_c1.shape[0]
spiketrains = []
for cluster in range(Ncl):
    spiketrains.append( neo.SpikeTrain( np.squeeze( 
        spike_times_c1[cluster] )*pq.s, t_stop=38.25*pq.s,
        sampling_rate=3e4 ) )
print("Spike train formatted. Starting SPADE:")

n_spikes = [len(spiketrain) for spiketrain in spiketrains]

bin_size = 0.5 * pq.ms

firing_rate_two_std = (np.mean(n_spikes) + 2 * np.std(n_spikes))/spiketrains[0].t_stop
firing_probability_two_std = (firing_rate_two_std * bin_size).simplified.item()

n_bins = int(spiketrains[0].t_stop/bin_size)

patterns = []
for min_spikes in range(2, 10):
    min_occ = math.ceil(firing_probability_two_std**min_spikes * n_bins)
    print(f'{min_spikes=}', f'{min_occ=}')
    spade_output = elephant.spade.spade(
        spiketrains=spiketrains,
        bin_size=bin_size,
        winlen=3,                   # 1 bin means synchronous patterns
        min_spikes=min_spikes,
        min_occ=min_occ,
        min_neu=min_occ,
        n_surr=100,
        dither=5*pq.ms)
    patterns.extend(spade_output['patterns'])
    print(f'Patterns found for {min_spikes} spikes: {len(spade_output["patterns"])}')

