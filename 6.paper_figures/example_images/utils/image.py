"""
Functions for preparing images for plotting
"""

import skimage
import numpy as np


def normalize_image(img):
    """
    Function copied from https://github.com/jr0th/segmentation/blob/master/visualization/CellArt.ipynb
    """
    # normalize to [0,1]
    percentile = 99
    high = np.percentile(img, percentile)
    low = np.percentile(img, 100 - percentile)

    img = np.minimum(high, img)
    img = np.maximum(low, img)

    # gives float64, thus cast to 8 bit later
    img = (img - low) / (high - low)

    img = skimage.img_as_ubyte(img)
    return img


def colorize_image(img, col):
    """
    Function copied from https://github.com/jr0th/segmentation/blob/master/visualization/CellArt.ipynb
    """
    # rescale image
    img_float = img.astype(np.float)
    img_float = img_float / 255

    # colorize
    img_col_float = np.reshape(img_float, img_float.shape + (1,)) * col
    img_col_byte = img_col_float.astype(np.uint8)

    return img_col_byte
