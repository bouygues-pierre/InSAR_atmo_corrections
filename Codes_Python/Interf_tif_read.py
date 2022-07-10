import numpy as np
from osgeo import gdal, gdal_array
import exifread
import matplotlib.pyplot as plt

root_name = '/Volumes/Pierre_2TO/Stage_M1_Albino/Data/Piton_de_la_Fournaise/GACOS_LiCSAR_geotiff/'
track = '151D_11112_000300'


def Read_Interf_tiff(root_name, track, date, master, slave):
    # cath the necessary path to load datas
    dir_name_tiff = '/'.join([track, date, 'Uncorrected_interf', 'unwrapped'])
    file = ''.join([root_name, '/', dir_name_tiff, '/', master, '_', slave, '.geo.unw.tif'])

    # read data with the gdal library to extract data array of delay map

    inter_file = gdal.Open(file)
    inter_band = inter_file.GetRasterBand(1)
    inter_array = gdal_array.LoadFile(file)
    # nodata = inter_band.GetNoDataValue()
    # inter_array = np.ma.masked_equal(inter_array, nodata)

    # reading data to extract exif of tif delay map
    inter_file_exif = open(file, 'rb')

    # Return Exif tags
    tags = exifread.process_file(inter_file_exif)

    # metadata storage in a dictionary call FC
    FC = {tag: tags[tag].values for tag in tags.keys()}

    # convert tuple to list
    FC['Image Tag 0x830E'] = list(sum(FC['Image Tag 0x830E'], ()))
    FC['Image Tag 0x8482'] = list(sum(FC['Image Tag 0x8482'], ()))
    # FC['Image Tag 0x87B0'] = list(sum(FC['Image Tag 0x87B0'], ()))

    # get values of interest
    j = FC.get('Image ImageWidth')[0]
    i = FC.get('Image ImageLength')[0]
    xi = FC.get('Image Tag 0x8482')[3]  # get X_first data
    yi = FC.get('Image Tag 0x8482')[4]  # get Y_first data
    dx = FC.get('Image Tag 0x830E')[0]  # get X step data
    dy = FC.get('Image Tag 0x830E')[1]  # get Y step data

    i_arr = np.arange(1, i, 1)
    j_arr = np.arange(1, j, 1)

    Lon = xi + dx * j_arr - dx / 2
    Lat = yi - dy * i_arr - dy / 2

    return Lon, Lat, dx, dy, inter_array, i, j


date = '20170813_20220425'
master = '20170813'
slave = '20170825'

Lon, Lat, dx, dy, inter_array, i, j = Read_Interf_tiff(root_name, track, date, master, slave)

# print(Lon)
# print(Lat)
print(inter_array.shape)
