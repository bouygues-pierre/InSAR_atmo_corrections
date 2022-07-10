import numpy as np
import matplotlib.pyplot as plt
import skgstat as skg
import Interf_tif_read as ITR
import Interf_cor_tif_read as ICTR
import Crop_function as CF
import Gacos_tif_read as GTR
import os
import itertools
import time
import csv
# from netCDF4 import Dataset
import concurrent.futures
from multiprocessing import Pool

start_time = time.time()

# ----------------------------------------------------------------------------------------------------------------------

#                                                  Data acquisition

# ----------------------------------------------------------------------------------------------------------------------

root_name = 'D:/Stage_M1_Albino/Data/Piton_de_la_Fournaise/GACOS_LiCSAR_geotiff/'
track = "151D_11112_000300"  # '144A_11114_000304' 151D_11112_000300
date_folder = '20170813_20220425'

files_inter = os.listdir("D:/Stage_M1_Albino/Data/Piton_de_la_Fournaise/GACOS_LiCSAR_geotiff"
                         "/151D_11112_000300/20170813_20220425/Uncorrected_interf/unwrapped")

# files_inter = os.listdir("/Volumes/Pierre_2TO/Stage_M1_Albino/Data/Piton_de_la_Fournaise/GACOS_LiCSAR_geotiff"
# "/151D_11112_000300/20190908_20211213/Corrected_interf")

# Gacos map reading to extract Lon_atm and Lat_atm to crop Inter_uncor withInter_cor
Lon_atm, Lat_atm, dx_atm, dy_atm, atmo1 = GTR.Read_GACOS_tiff(root_name, track, date_folder, '20171129')

files_inter = files_inter[1: -1]

with open('stats_variogram.csv', 'w', newline='') as csvfile:
    fieldnames = ["date", "range_cor", "sill_cor", "shape_cor", "nugget_cor", "range_uncor", "sill_uncor",
                  "shape_uncor", "nugget_uncor"]

    writer = csv.DictWriter(csvfile, fieldnames=fieldnames, quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writer.writeheader()

    # ------------------------------------------------------------------------------------------------------------------

    #                                                  data loop

    # ------------------------------------------------------------------------------------------------------------------

    for i in range(len(files_inter)):
        files_inter[i] = files_inter[i].split('.')[0]
        date = files_inter[i]
        master = files_inter[i].split('_')[0]
        slave = files_inter[i].split('_')[1]
        print(i, date)

        # Inter map reading
        Lon_inter, Lat_inter, dx_inter, dy_inter, inter_array, i_inter, j_inter = \
            ITR.Read_Interf_tiff(root_name, track, date_folder, master, slave)

        # Inter corrected map reading
        inter_cor_array = ICTR.Read_Interf_cor_tiff(root_name, track, date_folder, master, slave)

        inter_array, inter_cor_array = CF.crop_function(inter_array, inter_cor_array, Lat_inter, Lon_inter, Lat_atm,
                                                        Lon_atm)

        variogram_filename_SVG = 'D:/Stage_M1_Albino/Data/Piton_de_la_Fournaise/GACOS_LiCSAR_geotiff/' \
                                 '151D_11112_000300/20170813_20220425/Figures/variogram/SVG/variogram_{}.svg'.format(
            date)
        variogram_filename_PNG = 'D:/Stage_M1_Albino/Data/Piton_de_la_Fournaise/GACOS_LiCSAR_geotiff/' \
                                 '151D_11112_000300/20170813_20220425/Figures/variogram/PNG/variogram_{}.png'.format(
            date)

        # --------------------------------------------------------------------------------------------------------------

        #                                             parameter setup

        # --------------------------------------------------------------------------------------------------------------

        coordinates = a = np.array(list(itertools.product(range(inter_array.shape[0]), range(inter_array.shape[1]))))
        # setup of the coordinates array.Array of shape (m, n). Will be used as m observation points of n-dimensions.

        values_uncor = np.ones((coordinates.shape[0], 2))
        values_cor = np.ones((coordinates.shape[0], 2))
        # arrays in which the value (here in rad) of the pixel will be stored

        # ----------------------------------------------------------------------------------------------------------------------

        #                                       Interf_uncor variogram creation

        # ----------------------------------------------------------------------------------------------------------------------

        # for loop that associates a value, for each coordinates, in an array
        i = 0
        for c in coordinates:
            values_uncor[i] = inter_array[c[0], c[1]]
            i += 1

        # removal of NAN values
        coordinates_uncor = coordinates[np.isnan(values_uncor) == False]
        values_uncor = values_uncor[np.isnan(values_uncor) == False]

        # Isnan function change the size of array, we need to reshape them
        coordinates_uncor = np.reshape(coordinates_uncor, (int(coordinates_uncor.shape[0] / 2), 2))
        values_uncor = np.reshape(values_uncor, (int(values_uncor.shape[0] / 2), 2))
        values_uncor = values_uncor[:, 0]

        # variogramm creation with skgstat package : n_lags = number of bins, model = method used to fit experimental datas,
        # maxlag = the maximum value of distance that will be classify
        V_uncor = skg.Variogram(coordinates_uncor, values_uncor, estimator='matheron', use_nugget=True, n_lags=50,
                                model='stable', maxlag=300, samples=5000)

        # extraction of values of interest
        bins_uncor = V_uncor.bins
        coord_vario_uncor = V_uncor.coordinates
        data_uncor = V_uncor.data
        expe_data_uncor = V_uncor.experimental
        mean_res_uncor = V_uncor.mean_residual
        model_uncor = V_uncor.model_deviations()
        param_uncor = V_uncor.parameters
        res_uncor = V_uncor.residuals
        val_uncor = V_uncor.values

        # ----------------------------------------------------------------------------------------------------------------------

        #                                       Interf_cor variogram creation

        # ----------------------------------------------------------------------------------------------------------------------

        # for loop that associates a value, for each coordinates, in an array
        i = 0
        for c in coordinates:
            values_cor[i] = inter_cor_array[c[0], c[1]]
            i += 1

        # removal of NAN values
        coordinates_cor = coordinates[np.isnan(values_cor) == False]
        values_cor = values_cor[np.isnan(values_cor) == False]

        # Isnan function change the size of array, we need to reshape them
        coordinates_cor = np.reshape(coordinates_cor, (int(coordinates_cor.shape[0] / 2), 2))
        values_cor = np.reshape(values_cor, (int(values_cor.shape[0] / 2), 2))
        values_cor = values_cor[:, 0]

        # variogramm creation with skgstat package : n_lags = number of bins, model = method used to fit experimental datas,
        # maxlag = the maximum value of distance that will be classify
        V_cor = skg.Variogram(coordinates_cor, values_cor, estimator='matheron', use_nugget=True, n_lags=50,
                              model='stable', maxlag=300, samples=5000)

        bins_cor = V_cor.bins
        coord_vario_cor = V_cor.coordinates
        data_cor = V_cor.data
        expe_data_cor = V_cor.experimental
        mean_res_cor = V_cor.mean_residual
        model_cor = V_cor.model_deviations()
        param_cor = V_cor.parameters
        res_cor = V_cor.residuals
        val_cor = V_cor.values
        print("data created")
        # ----------------------------------------------------------------------------------------------------------------------

        #                                             Visualisation

        # ----------------------------------------------------------------------------------------------------------------------

        # variogram plot creation

        fig1, axes1 = plt.subplots(1, 1, figsize=(8, 6))
        axes1.plot(bins_uncor * 0.09, np.sqrt(expe_data_uncor), '.r', markersize=8,
                   label='Uncorrected data:\n• Sill = {}\n• range = {}\n• Nugget = {}'
                   .format(round(np.sqrt(param_uncor[1]), 2), round(param_uncor[0] * 0.09, 1),
                           round(np.sqrt(param_uncor[2]),
                                 2)))
        axes1.plot(bins_cor * 0.09, np.sqrt(expe_data_cor), '.b', markersize=8,
                   label='Corrected data\n• Sill: {}\n• range = {}\n• Nugget: {}'
                   .format(round(np.sqrt(param_cor[1]), 2), round(param_cor[0] * 0.09, 1),
                           round(np.sqrt(param_cor[2]), 2)))
        axes1.plot(bins_uncor * 0.09, np.sqrt((expe_data_uncor + res_uncor)), '-r', label='Uncorrected fit model')
        axes1.plot(bins_cor * 0.09, np.sqrt((expe_data_cor + res_cor)), '-b', label='Corrected fit model')
        axes1.legend(fontsize=8.5)
        axes1.set_title('Variogram of {}'.format(date))
        axes1.set_xlabel('Length in kilometers')
        axes1.set_ylabel('Variance square root (Matheron)')
        print("figures created")

        plt.savefig(variogram_filename_PNG, format='png', dpi=100)
        plt.savefig(variogram_filename_SVG, format='svg', dpi=100)

        print("figure saved")
        plt.close()

        writer.writerow({"date": date, "range_cor": param_cor[0], "sill_cor": param_cor[1], "shape_cor": param_cor[2],
                         "nugget_cor": param_cor[3], "range_uncor": param_uncor[0], "sill_uncor": param_uncor[1],
                         "shape_uncor": param_uncor[2], "nugget_uncor": param_uncor[3]})

        print("line wrote")

        print("--- %s seconds ---" % (time.time() - start_time))
