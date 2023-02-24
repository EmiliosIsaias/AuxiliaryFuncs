# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import cv2 as cv
import numpy as np
import pathlib as pl
import matplotlib.pyplot as plt

imgs_dir = pl.Path(r"H:\2022\Mexico\\")
pSub = 134
rng = 2
imgs_nm = ["/_MG_0{}.png".format(x) for x in range(pSub,pSub+rng)]
#imgs_nm = ["/_MG_9209.png","/_MG_9210.png","/_MG_9211.png"]
imgs_fn = [imgs_dir.as_posix() + im for im in imgs_nm]

img_list = [cv.imread(fn) for fn in imgs_fn]

alignMTB = cv.createAlignMTB()
alignMTB.process(img_list, img_list)

merge_mertens = cv.createMergeMertens()
res_mertens = merge_mertens.process(img_list)

res_mertens_8bit = np.clip(res_mertens*255, 0, 255).astype('uint24')
plt.imshow(res_mertens_8bit)

out_name = imgs_nm[0].partition(".")
cv.imwrite(imgs_dir.as_posix() + out_name[0] + "_HDR.png", res_mertens_8bit)