#import module
import cv2 as cv
import numpy as np
from matplotlib import pyplot as plt
import ipywidgets as widgets
from ipywidgets import interact, interact_manual

import argparse

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--sampleId", required=True, help="sampleId")
args = vars(ap.parse_args())

sampleId = args["sampleId"]

#read image
img = cv.imread(sampleId +'.jpg',0)

#filter non-tissue part
blur = cv.GaussianBlur(img,(5,5),0)
ret3,th3 = cv.threshold(blur,0,255,cv.THRESH_BINARY+cv.THRESH_OTSU)
numRows,numCols=th3.shape
pixel=50
pixel_count = 2*pixel-1;
pixel_w = numCols/pixel_count;
pixel_h = numRows/pixel_count;
useful_pixel ="";
filter_image=np.zeros((50,50),dtype=int)

#make position
for i in range(1,51):
    y = round(2*(i-1)*pixel_h+1);
    for j in range(1,51):
        x = round(2*(j-1)*pixel_w+1);
        pixel = th3[y:round(y+pixel_h-1),x:round(x+pixel_w-1)];
        C=sum(pixel);
        C=sum(C)
        if C > 0:
            useful_pixel=useful_pixel+","+str(j)+"x"+str(i)
            filter_image[i-1,j-1]=1;

#show figure after filtered
plt.figure(figsize=(10,8))
plt.subplot(1,2,1),plt.imshow(th3, "gray"),plt.title("Otsu's thresholding")
plt.subplot(1,2,2),plt.imshow(filter_image, "gray"),plt.title("bin")

# save position file
data = open("position_{}.txt".format(sampleId),'w',encoding="utf-8")
print(useful_pixel,file=data)
data.close()