import numpy as np
import matplotlib.pyplot as plt
import os
import struct
import io
import pandas as pd


def read_GACOS_tiff(rootname, track, date, master):
    dirname = '/'.join([track, date, 'GACOS_data/', 'Tiff_data'])
    file = ''.join([rootname, '/', dirname, '/', master, '.ztd.tif'])
    header = ''.join([rootname, '/', dirname, '/', master, '.ztd.rsc'])

    with open(header, "r") as fid:

        FC = {}

        for ligne in fid:
            element = ligne.split()
            FC_name = element[0]
            FC_val = element[1]
            FC[FC_name] = FC_val

    xi = float(FC.get('X_FIRST'))
    yi = float(FC.get('Y_FIRST'))
    dx = float(FC.get('X_STEP'))
    dy = float(FC.get('Y_STEP'))
    j = float(FC.get('WIDTH'))
    i = float(FC.get('FILE_LENGTH'))

    i_arr = np.arange(1, i, 1)
    j_arr = np.arange(1, j, 1)

    Lon = xi + dx * j_arr - dx/2
    Lat = yi + dy * i_arr - dy/2

    with open(file, 'w', encoding='utf-16-le') as d:
        delay = io.TextIOWrapper(d)
        delay = io.TextIOWrapper.read(delay)
        print(delay)
        print(d)

    #with open(file, 'r', encoding='utf-16-le') as d:
        #delay = d.encode('utf-16-le')
    #print(delay)


    #delay = io.open(file)
    #delay = io.StringIO(delay)
    #struct.unpack('>HHHH', delay)
    #print(delay)

    return Lon, Lat, dx, dy


Lon, Lat, dx, dy = read_GACOS_tiff(rootname, track, date, master)

#plt.figure()
#plt.imshow(atmo)
#plt.show()




#def get_delay_geo(ztd_file, atr, cos_inc_angle):
    #"""calc single path tropo delay in line-of-sight direction
    #Parameters: ztd_file      - str, path of zenith delay file
               # atr           - dict, dictionary of attribute for output file
              #  cos_inc_angle - 2D np.ndarray in float32, cos(inc_angle)
    #Returns:    delay         - 2D np.ndarray in float32, LOS delay
   # """
    # get geo_box from ts_file attributes
    #length, width = int(atr['LENGTH']), int(atr['WIDTH'])
    #geo_box = ut.coordinate(atr).box_pixel2geo(pixel_box=(0, 0, width, length))
    # geo_box --> pix_box in ztd file
    #atr_ztd = readfile.read_attribute(ztd_file)
    #pix_box = ut.coordinate(atr_ztd).box_geo2pixel(geo_box)
    # read ztd file
    #delay = readfile.read(ztd_file, box=pix_box)[0]
    # interpolate/resample into the same resolution as ts_file
    #delay = resize(delay, (length, width),
    #order=1,
    #mode='constant',
    #anti_aliasing=True,
    #preserve_range=True)
    # project from zenith to line-of-sight
    #delay /= cos_inc_angle
    # reverse the sign for consistency between different phase correction steps/methods
    #delay *= -1
    #return delay

