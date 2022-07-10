import numpy as np
import matplotlib.pyplot as plt
import skgstat as skg
import Interf_tif_read as ITR
import Interf_cor_tif_read as ICTR
import Crop_function as CF
import Gacos_tif_read as GTR
import os
import time
start_time = time.time()

# ----------------------------------------------------------------------------------------------------------------------

#                                                  Data acquisition

# ----------------------------------------------------------------------------------------------------------------------

root_name = '/Volumes/Pierre_2TO/Stage_M1_Albino/Data/Piton_de_la_Fournaise/GACOS_LiCSAR_geotiff/'
track = '151D_11112_000300'

files_inter = os.listdir("/Volumes/Pierre_2TO/Stage_M1_Albino/Data/Piton_de_la_Fournaise/GACOS_LiCSAR_geotiff"
                         "/151D_11112_000300/20190908_20211213/Uncorrected_interf/unwrapped")

# files_inter = os.listdir("/Volumes/Pierre_2TO/Stage_M1_Albino/Data/Piton_de_la_Fournaise/GACOS_LiCSAR_geotiff"
# "/151D_11112_000300/20190908_20211213/Corrected_interf")

# Gacos map reading to extract Lon_atm and Lat_atm to crop Inter_uncor withInter_cor
Lon_atm, Lat_atm, dx_atm, dy_atm, atmo1 = GTR.Read_GACOS_tiff(root_name, track, '20190908_20211213', '20190908')

files_inter = files_inter[1: -1]

for i in range(140, len(files_inter)):
    files_inter[i] = files_inter[i].split('.')[0]
    date = files_inter[i]
    master = files_inter[i].split('_')[0]
    slave = files_inter[i].split('_')[1]
    print(i, date)

    # Inter map reading
    Lon_inter, Lat_inter, dx_inter, dy_inter, inter_array, i_inter, j_inter = \
        ITR.Read_Interf_tiff(root_name, track, '20190908_20211213', master, slave)

    # Inter corrected map reading
    inter_cor_array = ICTR.Read_Interf_cor_tiff(root_name, track, '20190908_20211213', master, slave)

    variogram_filename = '/Volumes/Pierre_2TO/Stage_M1_Albino/Data/Piton_de_la_Fournaise/GACOS_LiCSAR_geotiff/' \
                         '151D_11112_000300/20190908_20211213/Figures/directionnal_variogram/variogram_dir_{}.png'\
                         .format(date)

    inter_array, inter_cor_array = CF.crop_function(inter_array, inter_cor_array, Lat_inter, Lon_inter, Lat_atm, Lon_atm)

# ----------------------------------------------------------------------------------------------------------------------

#                                             parameters setup

# ----------------------------------------------------------------------------------------------------------------------

# number of sample
    samples = 5000

# creation of the array of coordinates wich take random pairs of coordinates of the smalest side of the inter_array
    coords_uncor = np.random.randint((inter_array.shape[0]), size=(samples, 2))
    coords_cor = coords_uncor.copy()

# creation of the array which contain the values of the coordinates
    values_uncor = np.ones((len(coords_uncor)))
    values_cor = np.ones((len(coords_cor)))

    i = 0
    for c in coords_uncor:
        values_uncor[i] = inter_array[c[0], c[1]]
        values_cor[i] = inter_cor_array[c[0], c[1]]
        i += 1

# NaN suppression
    coords_cor = coords_cor[np.isnan(values_cor) == False]
    values_cor = values_cor[np.isnan(values_cor) == False]
    coords_uncor = coords_uncor[np.isnan(values_uncor) == False]
    values_uncor = values_uncor[np.isnan(values_uncor) == False]


# ----------------------------------------------------------------------------------------------------------------------

#                              Interf_uncor directionnal variogram creation

# ----------------------------------------------------------------------------------------------------------------------

    Vnorth_uncor = skg.DirectionalVariogram(coords_uncor, values_uncor, azimuth=90, tolerance=90, n_lags=50, maxlag=300,
                                            model='stable')
    bins_Vnorth_uncor = Vnorth_uncor.bins
    coord_vario_Vnorth_uncor = Vnorth_uncor.coordinates
    data_Vnorth_uncor = Vnorth_uncor.data
    expe_data_Vnorth_uncor = Vnorth_uncor.experimental
    mean_res_Vnorth_uncor = Vnorth_uncor.mean_residual
    model_Vnorth_uncor = Vnorth_uncor.model_deviations()
    param_Vnorth_uncor = Vnorth_uncor.parameters
    res_Vnorth_uncor = Vnorth_uncor.residuals
    val_Vnorth_uncor = Vnorth_uncor.values

    Veast_uncor = skg.DirectionalVariogram(coords_uncor,  values_uncor, azimuth=0, tolerance=90, n_lags=50, maxlag=300,
                                       model='stable')
    bins_Veast_uncor = Veast_uncor.bins
    coord_vario_Veast_uncor = Veast_uncor.coordinates
    data_Veast_uncor = Veast_uncor.data
    expe_data_Veast_uncor = Veast_uncor.experimental
    mean_res_Veast_uncor = Veast_uncor.mean_residual
    model_Veast_uncor = Veast_uncor.model_deviations()
    param_Veast_uncor = Veast_uncor.parameters
    res_Veast_uncor = Veast_uncor.residuals
    val_Veast_uncor = Veast_uncor.values

# ----------------------------------------------------------------------------------------------------------------------

#                                       Interf_cor variogram creation

# ----------------------------------------------------------------------------------------------------------------------

    Vnorth_cor = skg.DirectionalVariogram(coords_cor, values_cor, azimuth=90, tolerance=90, n_lags=50, maxlag=300,
                                      model='stable')
    bins_Vnorth_cor = Vnorth_cor.bins
    coord_vario_Vnorth_cor = Vnorth_cor.coordinates
    data_Vnorth_cor = Vnorth_cor.data
    expe_data_Vnorth_cor = Vnorth_cor.experimental
    mean_res_Vnorth_cor = Vnorth_cor.mean_residual
    model_Vnorth_cor = Vnorth_cor.model_deviations()
    param_Vnorth_cor = Vnorth_cor.parameters
    res_Vnorth_cor = Vnorth_cor.residuals
    val_Vnorth_cor = Vnorth_cor.values

    Veast_cor = skg.DirectionalVariogram(coords_cor,  values_cor, azimuth=0, tolerance=90, n_lags=50, maxlag=300,
                                         model='stable')
    bins_Veast_cor = Veast_cor.bins
    coord_vario_Veast_cor = Veast_cor.coordinates
    data_Veast_cor = Veast_cor.data
    expe_data_Veast_cor = Veast_cor.experimental
    mean_res_Veast_cor = Veast_cor.mean_residual
    model_Veast_cor = Veast_cor.model_deviations()
    param_Veast_cor = Veast_cor.parameters
    res_Veast_cor = Veast_cor.residuals
    val_Veast_cor = Veast_cor.values

# ----------------------------------------------------------------------------------------------------------------------

#                                             Visualisation

# ----------------------------------------------------------------------------------------------------------------------

    fix1, axes1 = plt.subplots(1, 1, figsize=(15, 10))

    axes1.plot(bins_Vnorth_uncor * 0.09, np.sqrt(expe_data_Vnorth_uncor), '.k', markersize=8,
               label='North-South uncorrected data : • Sill: {} • range = {} • Nugget: {}'
               .format(round(np.sqrt(param_Vnorth_uncor[1]), 2), round(param_Vnorth_uncor[0]*0.09, 1),
                       round(np.sqrt(param_Vnorth_uncor[2]), 2)))

    axes1.plot(bins_Vnorth_uncor * 0.09, np.sqrt(expe_data_Vnorth_uncor + res_Vnorth_uncor), '-k',
               label='North-South fit model')

    axes1.plot(bins_Veast_uncor * 0.09, np.sqrt(expe_data_Veast_uncor), '.g', markersize=8,
               label='East-West uncorrected data : • Sill: {} • range = {} • Nugget: {}'
               .format(round(np.sqrt(param_Veast_uncor[1]), 2), round(param_Veast_uncor[0]*0.09, 1),
                       round(np.sqrt(param_Veast_uncor[2]), 2)))

    axes1.plot(bins_Veast_uncor * 0.09, np.sqrt(expe_data_Veast_uncor + res_Veast_uncor), '-g',
               label='East-West fit model')

    axes1.plot(bins_Vnorth_cor * 0.09, np.sqrt(expe_data_Vnorth_cor), '.k', markersize=8, alpha=0.3,
               label='North-South corrected data : • Sill: {} • range = {} • Nugget: {}'
               .format(round(np.sqrt(param_Vnorth_cor[1]), 2), round(param_Vnorth_cor[0]*0.09, 1),
                   round(np.sqrt(param_Vnorth_cor[2]), 2)))

    axes1.plot(bins_Vnorth_cor * 0.09, np.sqrt(expe_data_Vnorth_cor + res_Vnorth_cor), '-k', alpha=0.4,
               label='North-South fit model')

    axes1.plot(bins_Veast_cor * 0.09, np.sqrt(expe_data_Veast_cor), '.g', markersize=8, alpha=0.3,
               label='East-West corrected data : • Sill: {} • range = {} • Nugget: {}'
               .format(round(np.sqrt(param_Veast_cor[1]), 2), round(param_Veast_cor[0]*0.09, 1),
                       round(np.sqrt(param_Veast_cor[2]), 2)))

    axes1.plot(bins_Veast_cor * 0.09, np.sqrt(expe_data_Veast_cor + res_Veast_cor), '-g', alpha=0.4,
               label='East-West fit model')

    axes1.set_title('Directionnal variogram of {}'.format(date))
    axes1.set_ylabel('Variance square root (matheron)')
    axes1.set_xlabel('length in kilometers')
    plt.legend()

    plt.savefig(variogram_filename, format='png', dpi=100)
    plt.close()
    print("--- %s seconds ---" % (time.time() - start_time))

