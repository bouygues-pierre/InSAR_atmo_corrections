import numpy as np
from osgeo import gdal, gdal_array
import exifread


def Read_GACOS_tiff(root_name, track, date, master):

    # cath the necessary path to load datas
    dir_name_tiff = '/'.join([track, date, 'Gacos_data', 'Tiff_data'])
    file = ''.join([root_name, '/', dir_name_tiff, '/', master, '.ztd.tif'])

    # read data with the gdal library to extract data array of delay map
    delay_file = gdal.Open(file)
    delay_band = delay_file.GetRasterBand(1)
    delay_array = gdal_array.LoadFile(file)
    nodata = delay_band.GetNoDataValue()
    delay_array = np.ma.masked_equal(delay_array, nodata)

    # reading data to extract exif of tif delay map
    delay_file_exif = open(file, 'rb')

    # Return Exif tags
    tags = exifread.process_file(delay_file_exif)

    # metadata storage in a dictionary call FC
    FC = dict()

    for tag in tags.keys():
        FC[tag] = tags[tag].values

    # convert tuple to list
    FC['Image Tag 0x830E'] = list(sum(FC['Image Tag 0x830E'], ()))
    FC['Image Tag 0x8482'] = list(sum(FC['Image Tag 0x8482'], ()))
    FC['Image Tag 0x87B0'] = list(sum(FC['Image Tag 0x87B0'], ()))

    # get values of interest
    j = FC.get('Image ImageWidth')[0]
    i = FC.get('Image ImageLength')[0]
    xi = FC.get('Image Tag 0x8482')[3]  # get X_first data
    yi = FC.get('Image Tag 0x8482')[4]  # get Y_first data
    dx = FC.get('Image Tag 0x830E')[0]  # get X step data
    dy = FC.get('Image Tag 0x830E')[1]  # get Y step data

    # Lon and Lat vector creation
    i_arr = np.arange(1, i, 1)
    j_arr = np.arange(1, j, 1)

    Lon = xi + dx * j_arr - dx/2
    Lat = yi - dy * i_arr - dy/2

    return Lon, Lat, dx, dy, delay_array

