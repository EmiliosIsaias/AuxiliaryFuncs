# -*- coding: utf-8 -*-
"""
Created on Fri May  3 11:10:19 2024

@author: neuro
"""
import deeplabcut as dlc
import pathlib as pl

config_path = pl.Path(r'Z:/Emilio\\AwakenSC-EIC-2021-07-20\\config.yaml')
config_path = config_path.as_posix()

batch_path = pl.Path(r'Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch10_ephys.e\RNs\WTg63\221109_PTX100microM_2000\Behaviour')

video_paths = list( batch_path.glob(r'**/*.avi') )
video_paths = [v.as_posix() for v in video_paths]

dlc.extract_outlier_frames( config_path, video_paths, 
                           outlieralgorithm='uncertain', p_bound=0.7 )