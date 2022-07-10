import numpy as np
from skimage.transform import resize


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
