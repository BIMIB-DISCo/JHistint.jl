using PyCall

py"""
import numpy as np
import matplotlib.pyplot as plt

from tifffile import imread, imsave
from csbdeep.utils import Path, normalize
from csbdeep.utils.tf import keras_import
keras = keras_import()

from stardist import export_imagej_rois, random_label_cmap
from stardist.models import StarDist2D

np.random.seed(0)
cmap = random_label_cmap()

from zipfile import ZIP_DEFLATED


img = imread('C:/Users/nicom/Desktop/Extra/segmentation/TCGA-OR-A5J1-01A-01-TS1.CFE08710-54B8-45B0-86AE-500D6E36D8A5_004.svs')
model = StarDist2D.from_pretrained('2D_versatile_he')
from csbdeep.data import Normalizer, normalize_mi_ma

class MyNormalizer(Normalizer):
    def __init__(self, mi, ma):
            self.mi, self.ma = mi, ma
    def before(self, x, axes):
        return normalize_mi_ma(x, self.mi, self.ma, dtype=np.float32)
    def after(*args, **kwargs):
        assert False
    @property
    def do_after(self):
        return False

mi, ma = 0, 255
normalizer = MyNormalizer(mi, ma)

# labels, polys = model.predict_instances(img, normalizer=normalizer, n_tiles=(32,32,1))
labels, polys = model.predict_instances_big(img, axes='YXC', block_size=2048, min_overlap=128, context=128,
                                            normalizer=normalizer, n_tiles=(4,4,1))

imsave('labels.tif', labels, compress=ZIP_DEFLATED)
"""
