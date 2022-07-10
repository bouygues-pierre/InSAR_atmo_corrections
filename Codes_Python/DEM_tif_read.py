import numpy as np
from osgeo import gdal, gdal_array
import exifread

root_name = '/Volumes/Pierre_2TO/Stage_M1_Albino/Data/Piton_de_la_Fournaise/GACOS_LiCSAR_geotiff/'
track = '151D_11112_000300'

def Read_DEM_tiff(root_name, track):

    dir_name_tiff = '/'.join([track])
    file = ''.join([root_name, '/', dir_name_tiff, '/', track, '.geo.hgt.tif'])

    return gdal_array.LoadFile(file)

'''''
    incid_file = gdal.Open(file)
    incid_band = incid_file.GetRasterBand(1)
    incid_array = gdal_array.LoadFile(file)
    nodata = incid_band.GetNoDataValue()
    inter_array = np.ma.masked_equal(incid_array, nodata)

    incid = incid_file.GetMetadata()
    print(incid)
    inter_file_exif = open(file, 'rb')

    tags = exifread.process_file(inter_file_exif)

    FC = dict()

    for tag in tags.keys():
        FC[tag] = tags[tag]

    print(FC)
    '''''



