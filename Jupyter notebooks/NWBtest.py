# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 17:27:24 2023

@author: neuro
"""
from pynwb import NWBHDF5IO, validate
from nwbwidgets import nwb2widget

fileName = r"Z:\\SC Anatomy paper data\\Roller\\pyNWB\\GADe59_221101.nwb"
io = NWBHDF5IO(fileName, mode='r')
validate(io)

nwb = io.read()

nwb2widget(nwb)

io.close()


