import matplotlib.image as mpimg
import matplotlib.image
import numpy as np


img = mpimg.imread('figures/levioSAM_S_bw.png')
img2 = []
for irow, row in enumerate(img):
    img2.append([])
    for cell in row:
        if sum(cell) != 4:
            img2[-1].append([1, 1, 1, 1])
        else:
            img2[-1].append([1-cell[0], 1-cell[1],  1-cell[2],  1.])
matplotlib.image.imsave('figures/levioSAM_S_bw_dark.png', img2, dpi=300)

