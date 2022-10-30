# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 12:46:06 2022

@author: neuro
"""
import deeplabcut as dlc
import pathlib as pl

config_path = pl.Path(r"E:\DLC\AwakenSC-EIC-2021-07-20\config.yaml")
config_path = config_path.as_posix()

batch_path = pl.Path(r"Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch10_ephys.e")
videos_path = list(batch_path.glob('**/*.avi'))
videos_path = [v.as_posix() for v in videos_path]

dlc.analyze_videos(config_path, videos_path, save_as_csv=True)
[dlc.filterpredictions(config_path, v, save_as_csv=True) for v in videos_path]
