# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 11:56:29 2024

@author: neuro
"""
import spikeinterface.full as si
import pathlib as pl
import probeinterface as pint
from probeinterface.plotting import plot_probe

import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import os


eph_dir = pl.Path(r"Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch18_ephys\iRNs\GADi40\240208_C+F_2350\ephys_E1")
                 
rec_files = [x.as_posix() for x in list( eph_dir.glob('./GADi*.bin') )]
fs = 3e4;

recs = si.read_binary(rec_files, fs, dtype='int16', num_channels=64, 
                      gain_to_uV=0.195, offset_to_uV=0., 
                      is_filtered=False)

# Probe E1
probe = pint.get_probe(manufacturer='cambridgeneurotech', 
                       probe_name='ASSY-77-E-1')
# ASSY-77>ADPT.A64-Om32x2-CN>RHD2164
connection = pint.get_available_pathways()[-1]
probe.wiring_to_device(connection)

recs = recs.set_probe(probe=probe)

recs_f = si.bandpass_filter(recs, freq_min=600, freq_max=10e3, 
                                    dtype='float32')

recs_cmr = si.common_reference(recs, reference='global', operator='median')

si.plot_traces({"filtered": recs_f, "common": recs_cmr}, 
                   time_range=[60.15, 60.4])

bad_channel_ids, channel_labels = si.detect_bad_channels(
    recs_f, method='coherence+psd')

