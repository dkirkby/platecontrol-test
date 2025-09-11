"""Creates a set of simulated move requests for testing the anticollision code.
"""

import os
import sys
import csv
import pandas as pd
import numpy as np
import math
from scipy.ndimage.interpolation import rotate
import matplotlib.pyplot as plt

def rotmat2D_deg(angle):
    """Return the 2d rotation matrix for an angle given in degrees."""
    angle *= math.pi/180.
    return np.array([[math.cos(angle), -math.sin(angle)], [math.sin(angle), math.cos(angle)]])

## Read Targets ##
target_file='/home/msdos/fiberassign_collision_check/tile-060226.fits'
import fitsio
d=fitsio.FITS(target_file)
obsX_full=d[1]['FIBERASSIGN_X'][:]
obsY_full=d[1]['FIBERASSIGN_Y'][:]
petal_loc_full=d[1]['PETAL_LOC'][:]
device_loc_full=d[1]['DEVICE_LOC'][:]

for i in range(10):
    ind_this_petal=np.where(petal_loc_full == i)
    nn=len(ind_this_petal[0])
    obsX_this=obsX_full[ind_this_petal]
    obsY_this=obsY_full[ind_this_petal]
    obsXY_this=[]
    for j in range(nn):
        obsXY_this.append([obsX_this[j],obsY_this[j]])
    device_loc_this=device_loc_full[ind_this_petal]
    ind_sort=np.argsort(device_loc_this)
    ####### Sorting #########
    device_loc_this=device_loc_this[ind_sort]
    obsX_this=obsX_this[ind_sort]
    obsY_this=obsY_this[ind_sort]
    obsXY_this=np.array(obsXY_this)[ind_sort]

    ####### Rotation to petal3 position ########
    xy_np = np.transpose(obsXY_this)
    rot = rotmat2D_deg((3-i)*36.)
    xy_np = np.dot(rot, xy_np) 
    plt.scatter(xy_np[0],xy_np[1])
    for j in range(len(device_loc_this)):
        plt.text(xy_np[0][j],xy_np[1][j],str(device_loc_this[j]))
    plt.show()
    #plt.scatter(xy_np[0],xy_np[1])
    #plt.show()
    command=[]
    list_list=[]
    for j in range(nn):
        command.append('obsXY')
        list_list.append([device_loc_this[j],'obsXY',xy_np[0][j],xy_np[1][j]])
    df = pd.DataFrame(list_list, columns = ['DEVICE_LOC','command','u','v'])
    
    df.to_csv('move_request_sets/requests_0000'+str(i).strip()+'.csv',index=False)

