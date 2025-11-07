# -*- coding: utf-8 -*-
"""
Created on Tue Oct 14 17:03:29 2025

@author: Puercos
"""

import cv2 as cv
import pathlib as pl
import os.path as op
from concurrent.futures import ThreadPoolExecutor
# import functools

imgs_dir = pl.Path(r"C:\Users\Puercos\Pictures")
mb_den = 1024**2
img_paths = list(imgs_dir.glob(r'20*\**\*.png'))  # Convert to list for reuse

def reduce_bit_depth_of_images(img_path):
    print(img_path.parts[-2:])
    f_size = op.getsize(img_path) / mb_den
    if f_size > 40:
        print(f'File is {f_size:.2f} MB! Reducing size...')
        im = cv.imread(img_path, cv.IMREAD_COLOR)
        if im is None:
            print(f'Failed to read {img_path.name}')
            return
        if not cv.imwrite(img_path, im):
            print(f'Failed to write {img_path.name}')
        else:
            new_size = op.getsize(img_path) / mb_den
            print(f'Done, new weight: {new_size:.2f} MB')
    else:
        print(f'File is {f_size:.2f} MB. No reduction necessary...')

# Use ThreadPoolExecutor to parallelize
with ThreadPoolExecutor() as executor:
    executor.map(reduce_bit_depth_of_images, img_paths)
