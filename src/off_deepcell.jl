using PyCall
using Images

link = joinpath(@__DIR__, "..", "input_example_demo", "slideExample1", "SlideExample_mini_1.tif")

py"""
from skimage.io import imread
from skimage import exposure, color
from deepcell.applications import Mesmer
import numpy as np
import cv2
import matplotlib.pyplot as plt

image = cv2.imread('SlideExample_mini_1.tif')

gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)

# Applica la segmentazione Otsu per ottenere una maschera binaria
_, mask = cv2.threshold(gray, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)

# Estrai i nuclei (regioni bianche) dalla maschera binaria
nuclei = cv2.bitwise_and(image, image, mask=mask)

# Estrai il citoplasma (regioni nere) dalla maschera binaria
cytoplasm = cv2.bitwise_and(image, image, mask=~mask)

# Imposta i canali dei nuclei su blu
nuclei[:, :, 0] = 0
nuclei[:, :, 2] = 0

# Imposta i canali del citoplasma su verde
cytoplasm[:, :, 0] = 0
cytoplasm[:, :, 2] = 0

# Combina nuclei e citoplasma in un'immagine multiplexed
multiplexed_image = cv2.add(nuclei, cytoplasm)



plt.imshow(multiplexed_image)
plt.axis('off')  # Rimuovi gli assi
plt.show()
# Create the application
#app = Mesmer()

# create the lab
#labeled_image = app.predict(im)
"""
