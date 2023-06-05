# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 19:38:58 2023

@author: jefe_
"""

import cv2 as cv
import pathlib as pl
import matplotlib.pyplot as plt

imgs_dir = pl.Path(r"C:\Users\jefe_\06-Heidelberg, Schwarzwald")
pSub = 620
rng = 5
imgs_nm = ["/_MG_1{}.png".format(x) for x in range(pSub,pSub+rng)]
#imgs_nm = ["/_MG_9209.png","/_MG_9210.png","/_MG_9211.png"]
imgs_fn = [imgs_dir.as_posix() + im for im in imgs_nm]

img_list = [cv.imread(fn) for fn in imgs_fn]

alignMTB = cv.createAlignMTB()
alignMTB.process(img_list, img_list[0])

image_data = []
for img in images:
    this_image = cv2.imread(img, 1)
    image_data.append(this_image)

avg_image = image_data[0]
for i in range(len(image_data)):
    if i == 0:
        pass
    else:
        alpha = 1.0/(i + 1)
        beta = 1.0 - alpha
        avg_image = cv2.addWeighted(image_data[i], alpha, avg_image, beta, 0.0)

cv2.imwrite('avg_happy_face.png', avg_image)
avg_image = cv2.imread('avg_happy_face.png')
plt.imshow(avg_image)
plt.show()