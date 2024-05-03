# -*- coding: utf-8 -*-
"""
Created on Fri May  3 11:10:19 2024

@author: neuro
"""
import platform
import deeplabcut as dlc
import pathlib as pl

os_name = platform.system()

if os_name == 'Linux':
    path_base = '/mnt/sds-hd/sd19b001/'
else:
    path_base = 'Z:\\';

config_path = pl.Path(path_base + 'Emilio/AwakenSC-EIC-2021-07-20/config.yaml')
config_path = config_path.as_posix()

batch_path = pl.Path(r'./')

video_paths = list( batch_path.glob(r'**/*.avi') )
video_paths = [v.as_posix() for v in video_paths]

dlc.extract_outlier_frames( config_path, video_paths, 
                           outlieralgorithm='uncertain', p_bound=0.5 )