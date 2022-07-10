from scipy.stats import linregress
from osgeo import osr
import numpy as np
from skimage.transform import resize
from osgeo import gdal, gdal_array
import exifread



# ----------------------------------------------------------------------------------------------------------------------

#                                              Array to raster function

# ----------------------------------------------------------------------------------------------------------------------

def array2raster(newRasterfn, rasterOrigin, pixelWidth, pixelHeight, array):

    if array.dtype == np.float32:
        arr_type = gdal.GDT_Float32
    else:
        arr_type = gdal.GDT_Int32

    cols = array.shape[1]
    rows = array.shape[0]
    originX = rasterOrigin[0]
    originY = rasterOrigin[1]

    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(newRasterfn, cols, rows, 1, arr_type)
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(array)
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromEPSG(4326)
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()


# ----------------------------------------------------------------------------------------------------------------------

#                                              Crop function

# ----------------------------------------------------------------------------------------------------------------------

def crop_function(inter_array, atmo, Lat_inter, Lon_inter, Lat_atm, Lon_atm):
    # crop interfero to have the same size at atmo

    imin = np.absolute(Lat_inter - np.min(Lat_atm))
    index_imin = imin.argmin()

    imax = np.absolute(Lat_inter - np.max(Lat_atm))
    index_imax = imax.argmin()

    jmin = np.absolute(Lon_inter - np.min(Lon_atm))
    index_jmin = jmin.argmin()

    jmax = np.absolute(Lon_inter - np.max(Lon_atm))
    index_jmax = jmax.argmin()

    inter_array = inter_array[index_imax:index_imin, index_jmin:index_jmax]
    atmo = resize(atmo, inter_array.shape)
    #if atmo.size > inter_array.size:
        #inter_array = inter_array[index_imax:index_imin, index_jmin:index_jmax]
       # atmo = resize(atmo, inter_array.shape)

   # elif atmo.size < inter_array.size:
        #atmo = atmo[index_imax:index_imin, index_jmin:index_jmax]
        #inter_array = resize(inter_array, atmo.shape)

    return inter_array, atmo


# ----------------------------------------------------------------------------------------------------------------------

#                                              DEM reader function

# ----------------------------------------------------------------------------------------------------------------------

def Read_DEM_tiff(root_name, track):

    dir_name_tiff = '/'.join([track])
    file = ''.join([root_name, '/', dir_name_tiff, '/', track, '.geo.hgt.tif'])

    DEM_array = gdal_array.LoadFile(file)

    return DEM_array


# ----------------------------------------------------------------------------------------------------------------------

#                                              Gacos reader function

# ----------------------------------------------------------------------------------------------------------------------

def Read_GACOS_tiff(root_name, track, date, master):

    # cath the necessary path to load datas
    dir_name_tiff = '/'.join([track, date, 'GACOS_data', 'Tiff_data'])
    file = ''.join([root_name, '/', dir_name_tiff, '/', master, '.ztd.tif'])

    # read data with the gdal library to extract data array of delay map
    delay_file = gdal.Open(file)
    delay_band = delay_file.GetRasterBand(1)
    delay_array = gdal_array.LoadFile(file)
    #nodata = delay_band.GetNoDataValue()
    #delay_array = np.ma.masked_equal(delay_array, nodata)

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


# ----------------------------------------------------------------------------------------------------------------------

#                                      Incidence angle array reader function

# ----------------------------------------------------------------------------------------------------------------------

def Read_Incidence_tiff(root_name, track):

    dir_name_tiff = '/'.join([track])
    file = ''.join([root_name, '/', dir_name_tiff, '/', track, '.geo.U.tif'])

    incid_array = gdal_array.LoadFile(file)

    return incid_array

# ----------------------------------------------------------------------------------------------------------------------

#                                      interfero Tif reader function

# ----------------------------------------------------------------------------------------------------------------------

def Read_Interf_tiff(root_name, track, date, master, slave):

    # cath the necessary path to load datas
    dir_name_tiff = '/'.join([track, date, 'Unwrapped_interf'])
    file = ''.join([root_name, '/', dir_name_tiff, '/', master, '_', slave, '.geo.unw.tif'])

    # read data with the gdal library to extract data array of delay map
    inter_file = gdal.Open(file)
    inter_band = inter_file.GetRasterBand(1)
    inter_array = gdal_array.LoadFile(file)
    #nodata = inter_band.GetNoDataValue()
    #inter_array = np.ma.masked_equal(inter_array, nodata)

    # reading data to extract exif of tif delay map
    inter_file_exif = open(file, 'rb')

    # Return Exif tags
    tags = exifread.process_file(inter_file_exif)

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

    i_arr = np.arange(1, i, 1)
    j_arr = np.arange(1, j, 1)

    Lon = xi + dx * j_arr - dx/2
    Lat = yi - dy * i_arr - dy/2

    return Lon, Lat, dx, dy, inter_array, i, j


# ----------------------------------------------------------------------------------------------------------------------

#                                      Stat calculation function

# ----------------------------------------------------------------------------------------------------------------------

def stats(X, Y):

    X2 = X[np.logical_not(np.isnan(X)) & np.logical_not(np.isnan(Y))]
    Y2 = Y[np.logical_not(np.isnan(X)) & np.logical_not(np.isnan(Y))]

    # coefficietn directeur et ordonnée à l origine
    coef = linregress(X2, Y2)
    y = coef[0] * X2 + coef[1]

    return coef, y, X2
