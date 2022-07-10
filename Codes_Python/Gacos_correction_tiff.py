import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import Gacos_tif_read as GTR
import Interf_tif_read as ITR
import array2raster as a2r
import Crop_function as CF
import Incidence_tif_read as INTR
import DEM_tif_read as DEM
import Stats as st
import os
import seaborn as sns
import matplotlib.ticker as mticker
import csv
from skimage.transform import resize
import time

start_time = time.time()

# ----------------------------------------------------------------------------------------------------------------------

#                                               Parameters definition

# ----------------------------------------------------------------------------------------------------------------------
# function tutoriel
# Gacos_tif_read imput : root_name, track, date, master and return Lon, Lat, dx, dy, delay_array, nodata
# Interf_tif_read imput : root_name, track, date, master and return Lon, Lat, dx, dy, inter_array, nodata
# Crop_function input : inter_array, atmo, Lat_inter, Lon_inter, Lat_atm, Lon_atm
# Incidence_tif_read input : root_name, track

root_name = '/Volumes/Pierre_2TO/Stage_M1_Albino/Data/Piton_de_la_Fournaise/GACOS_LiCSAR_geotiff/'
track = '144A_11114_000304'  # '151D_11112_000300'
# date = '20191118_20220307'
# master = '20191118'
# slave = '20191130'
textfile = open("stats_144A_11114_000304.txt", 'w', newline='')
csvfile = open('stats_144A_11114_000304.csv', 'w', newline='')
fieldnames = ['Inter', 'Quality_coef', 'stdi', 'stdf', 'Quality_coef_Pit', 'stdi_Pit', 'stdf_Pit', 'r_value_inter_atmo',
              'r_value_DEM_atmo', 'r_value_DEM_inter', 'r_value_DEM_corr', 'a_inter_atmo', 'a_DEM_atmo', 'a_DEM_inter',
              'a_DEM_intcorr']
writer = csv.DictWriter(csvfile, fieldnames=fieldnames, quotechar='|', quoting=csv.QUOTE_MINIMAL)

writer.writeheader()

files_inter = os.listdir("/Volumes/Pierre_2TO/Stage_M1_Albino/Data/Piton_de_la_Fournaise/GACOS_LiCSAR_geotiff"
                         "/144A_11114_000304/20191118_20220307/Uncorrected_interf/unwrapped")

files_inter = files_inter[1: -1]

for i in range(len(files_inter)):
    files_inter[i] = files_inter[i].split('.')[0]
    date = files_inter[i]
    master = files_inter[i].split('_')[0]
    slave = files_inter[i].split('_')[1]
    print(date)
    print(i)

    # int_corr_filename = '/Volumes/Pierre_2TO/Stage_M1_Albino/Data/Piton_de_la_Fournaise/GACOS_LiCSAR_geotiff' \
    # '/144A_11114_000304/20191118_20220307/Corrected_interf/{}_corrected_interf.tif'.format(date)

    # ------------------------------------------------------------------------------------------------------------------

    #                                               Data acquisition

    # --------------------------------------------------------------------------------------------------------------------

    # delay difference map reading
    Lon_atm, Lat_atm, dx_atm, dy_atm, atmo1 = GTR.Read_GACOS_tiff(root_name, track, '20191118_20220307', master)
    Lon_atm, Lat_atm, dx_atm, dy_atm, atmo2 = GTR.Read_GACOS_tiff(root_name, track, '20191118_20220307', slave)

    atmo = (atmo2 - atmo1)

    # Inter map reading
    Lon_inter, Lat_inter, dx_inter, dy_inter, inter_array, i_inter, j_inter = \
        ITR.Read_Interf_tiff(root_name, track, '20191118_20220307', master, slave)

    # Incidence array reading
    incidence_angle = np.arccos(INTR.Read_Incidence_tiff(root_name, track))

    # DEM array reading
    DEM_array = DEM.Read_DEM_tiff(root_name, track)
    DEM_array[DEM_array == 0] = np.nan

    # Crop Image
    Lon_DEM = 55.00050000000001 + np.arange(1, 943, 1) * 0.0010000000000000002 - 0.0010000000000000002 / 2
    Lat_DEM = -20.7601111 - np.arange(1, 727, 1) * 0.0010000000000000037 - 0.0010000000000000037 / 2
    # Pour le descendante remplacer lat/lon inter par lat/lon DEM dans incidence anle et DEM
    inter_array, atmo_array = CF.crop_function(inter_array, atmo, Lat_inter, Lon_inter, Lat_atm, Lon_atm)
    incidence_angle, atmo_array = CF.crop_function(incidence_angle, atmo, Lat_inter, Lon_inter, Lat_atm, Lon_atm)
    DEM_array, atmo_array = CF.crop_function(DEM_array, atmo, Lat_inter, Lon_inter, Lat_atm, Lon_atm)
    inter_array = resize(inter_array, atmo_array.shape)

    # difference from ref and reprojection of zenith to LOS

    wavelength = 0.055465
    m2rad = 4 * np.pi / wavelength
    rad2m = (wavelength / (4 * np.pi))
    incidence_angle[incidence_angle == 1.5707964] = np.nan
    zen2los = 1 / np.cos(incidence_angle / 180 * np.pi)

    inter_array[inter_array == 0] = np.nan
    inter_array = inter_array - np.nanmean(inter_array)
    atmo_array = atmo_array - np.nanmean(atmo)
    atmo_array = atmo_array * zen2los * m2rad

    # Inter_corr calculation

    int_corr = inter_array - atmo_array
    stdi = np.nanstd(inter_array)
    stdf = np.nanstd(int_corr)
    quality_coef = (stdi - stdf) / stdi

    # Inter_corr tif writing
    # a2r.array2raster(int_corr_filename, [i_inter, j_inter], dx_inter, dy_inter, int_corr)

    # Piton coordonates

    Piton_uncor = inter_array[300: 460, 480: 640]
    Piton_cor = int_corr[300: 460, 480: 640]

    stdi_Pit = np.nanstd(Piton_uncor)
    stdf_Pit = np.nanstd(Piton_cor)
    quality_coef_Pit = (stdi_Pit - stdf_Pit) / stdi_Pit

    # csv and txt files writing

    coef_inter_atmo, y_inter_atmo, X2_1 = st.stats(inter_array, atmo_array)
    coef_DEM_atmo, y_DEM_atmo, X2_2 = st.stats(DEM_array, atmo_array)
    coef_DEM_inter, y_DEM_inter, X2_3 = st.stats(DEM_array, inter_array)
    coef_DEM_corr, y_DEM_corr, X2_4 = st.stats(DEM_array, int_corr)

    textfile.write('Inter : {}, Quality_coef : {}, stdi : {}, stdf : {}, \n coef_inter_atmo : {}, coef_DEM_atmo : {},'
                   ' coef_DEM_inter : {}, coef_DEM_corr : {} \n'.format(date, quality_coef, stdi, stdf,
                                                                        coef_inter_atmo[2], coef_DEM_atmo[2],
                                                                        coef_DEM_inter[2], coef_DEM_corr[2]))

    writer.writerow({'Inter': date, 'Quality_coef': quality_coef, 'stdi': stdi, 'stdf': stdf,
                     'Quality_coef_Pit': quality_coef_Pit, 'stdi_Pit': stdi_Pit, 'stdf_Pit': stdf_Pit,
                     'r_value_inter_atmo': coef_inter_atmo[2], 'r_value_DEM_atmo': coef_DEM_atmo[2],
                     'r_value_DEM_inter': coef_DEM_inter[2], 'r_value_DEM_corr': coef_DEM_corr[2],
                     'a_inter_atmo': coef_inter_atmo[0], 'a_DEM_atmo': coef_DEM_atmo[0],
                     'a_DEM_inter': coef_DEM_inter[0], 'a_DEM_intcorr': coef_DEM_corr[0]})

    # ----------------------------------------------------------------------------------------------------------------------

    #                                                   Visualisation

    # ----------------------------------------------------------------------------------------------------------------------

    # axes parameters
    min_lat_inter = np.min(Lat_inter)
    max_lat_inter = np.max(Lat_inter)
    min_lon_inter = np.min(Lon_inter)
    max_lon_inter = np.max(Lon_inter)

    min_lat_atm = np.min(Lat_atm)
    max_lat_atm = np.max(Lat_atm)
    min_lon_atm = np.min(Lon_atm)
    max_lon_atm = np.max(Lon_atm)

    Q = np.round(quality_coef, 3)
    fontsize = 13
    figsize = (10, 10)

    # #color map
    #
    # cmap = mpl.cm.jet
    # norm_unwra = mpl.colors.Normalize(vmin=-20, vmax=20)
    # norm_wra = mpl.colors.Normalize(vmin=0, vmax=7)
    #
    # plt.figure(1, figsize=figsize)
    # plt.imshow(inter_array, cmap=cmap, norm=norm_unwra)
    # axes = plt.gca()
    # ticks_loc = axes.get_xticks().tolist()
    # axes.xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
    # axes.xaxis.set_ticklabels([np.round(Lat_inter[int(x)], 2) for x in ticks_loc], fontsize=fontsize, rotation=60)
    # ticks_loc = axes.get_yticks().tolist()
    # axes.yaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
    # axes.yaxis.set_ticklabels([np.round(Lon_inter[int(x)], 2) for x in ticks_loc], fontsize=fontsize)
    # plt.colorbar()
    # plt.title('Interferogram_{}_{}'.format(master, slave))
    # plt.savefig('/Volumes/Pierre_2TO/Stage_M1_Albino/Data/Piton_de_la_Fournaise/GACOS_LiCSAR_geotiff/144A_11114_000304/'
    #             '20191118_20220307/Figures/SVG/images/Inter/Unw_im_{}.svg'.format(date), format='svg', dpi=100)
    #
    # plt.figure(2, figsize=figsize)
    # plt.imshow(atmo_array, cmap=cmap, norm=norm_unwra)
    # axes = plt.gca()
    # ticks_loc = axes.get_xticks().tolist()
    # axes.xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
    # axes.xaxis.set_ticklabels([np.round(Lat_inter[int(x)], 2) for x in ticks_loc], fontsize=fontsize, rotation=60)
    # ticks_loc = axes.get_yticks().tolist()
    # axes.yaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
    # axes.yaxis.set_ticklabels([np.round(Lon_inter[int(x)], 2) for x in ticks_loc], fontsize=fontsize)
    # plt.title('GACOS_model_{}_{}'.format(master, slave))
    # plt.colorbar()
    # plt.savefig('/Volumes/Pierre_2TO/Stage_M1_Albino/Data/Piton_de_la_Fournaise/GACOS_LiCSAR_geotiff/144A_11114_000304/'
    #             '20191118_20220307/Figures/SVG/images/Gacos_inter/Gac_im_{}.svg'.format(date), format='svg', dpi=100)
    #
    # plt.figure(3, figsize=figsize)
    # plt.imshow(int_corr, cmap=cmap, norm=norm_unwra)
    # axes = plt.gca()
    # ticks_loc = axes.get_xticks().tolist()
    # axes.xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
    # axes.xaxis.set_ticklabels([np.round(Lat_inter[int(x)], 2) for x in ticks_loc], fontsize=fontsize, rotation=60)
    # ticks_loc = axes.get_yticks().tolist()
    # axes.yaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
    # axes.yaxis.set_ticklabels([np.round(Lon_inter[int(x)], 2) for x in ticks_loc], fontsize=fontsize)
    # plt.title('Corrected_Interferogram_{}_{}'.format(master, slave))
    # plt.text(12, 17, 'Q : {}'.format(str(Q), fontsize=15))
    # plt.colorbar()
    # plt.savefig('/Volumes/Pierre_2TO/Stage_M1_Albino/Data/Piton_de_la_Fournaise/GACOS_LiCSAR_geotiff/144A_11114_000304/'
    #             '20191118_20220307/Figures/SVG/images/Inter_corr/Inter_corr_im_{}.svg'.format(date), format='svg', dpi=100)
    #
    # plt.figure(4, figsize=figsize)
    # sns.histplot(inter_array[~np.isnan(inter_array)], color='r', pmax=0)
    # sns.histplot(atmo_array[~np.isnan(atmo_array)], color='b', pmax=0)
    # plt.legend(['Interferegram', 'GACOS'])
    # plt.savefig('/Volumes/Pierre_2TO/Stage_M1_Albino/Data/Piton_de_la_Fournaise/GACOS_LiCSAR_geotiff/144A_11114_000304/'
    #             '20191118_20220307/Figures/SVG/plot/histograms/histo_{}.svg'.format(date), format='svg', dpi=100)
    #
    # plt.figure(5, figsize=figsize)
    # plt.imshow(int_corr % (2 * np.pi), cmap=cmap, norm=norm_wra)
    # axes = plt.gca()
    # ticks_loc = axes.get_xticks().tolist()
    # axes.xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
    # axes.xaxis.set_ticklabels([np.round(Lat_inter[int(x)], 2) for x in ticks_loc], fontsize=fontsize, rotation=60)
    # ticks_loc = axes.get_yticks().tolist()
    # axes.yaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
    # axes.yaxis.set_ticklabels([np.round(Lon_inter[int(x)], 2) for x in ticks_loc], fontsize=fontsize)
    # plt.title('Corrected_wrapped_Interferogram_{}_{}'.format(master, slave))
    # plt.colorbar()
    # plt.savefig('/Volumes/Pierre_2TO/Stage_M1_Albino/Data/Piton_de_la_Fournaise/GACOS_LiCSAR_geotiff/144A_11114_000304/'
    #             '20191118_20220307/Figures/SVG/images/Inter_wra/Inter_wra_{}.svg'.format(date), format='svg', dpi=100)
    #
    # plt.figure(6, figsize=figsize)
    # plt.imshow(atmo_array % (2 * np.pi), cmap=cmap, norm=norm_wra)
    # axes = plt.gca()
    # ticks_loc = axes.get_xticks().tolist()
    # axes.xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
    # axes.xaxis.set_ticklabels([np.round(Lat_inter[int(x)], 2) for x in ticks_loc], fontsize=fontsize, rotation=60)
    # ticks_loc = axes.get_yticks().tolist()
    # axes.yaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
    # axes.yaxis.set_ticklabels([np.round(Lon_inter[int(x)], 2) for x in ticks_loc], fontsize=fontsize)
    # plt.title('Wrapped_GACOS_model_{}_{}'.format(master, slave))
    # plt.colorbar()
    # plt.savefig('/Volumes/Pierre_2TO/Stage_M1_Albino/Data/Piton_de_la_Fournaise/GACOS_LiCSAR_geotiff/144A_11114_000304/'
    #             '20191118_20220307/Figures/SVG/images/Gacos_inter_wra/Gacos_inter_wra_{}.svg'.format(date), format='svg'
    #             , dpi=100)
    #
    # plt.figure(7, figsize=figsize)
    # coef_inter_atmo, y_inter_atmo, X2_1 = st.stats(inter_array, atmo_array)
    # plt.plot(inter_array, atmo_array, '.r')
    # plt.plot(X2_1, y_inter_atmo, '-b')
    # plt.xlabel('LOS change (rad)')
    # plt.ylabel('GACOS model (rad)')
    # plt.ylim([-20, 20])
    # plt.text(10, 15, 'r value: {}'.format(np.round(coef_inter_atmo[2], 2)), ha='left', va='top', fontsize=fontsize)
    # plt.savefig('/Volumes/Pierre_2TO/Stage_M1_Albino/Data/Piton_de_la_Fournaise/GACOS_LiCSAR_geotiff/144A_11114_000304/'
    #             '20191118_20220307/Figures/SVG/plot/Inter_atmo/Inter_atmo_{}.svg'.format(date), format='svg', dpi=100)

    # plt.figure(8, figsize=figsize)
    # coef_DEM_atmo, y_DEM_atmo, X2_2 = st.stats(DEM_array, atmo_array)
    # plt.plot(DEM_array, atmo_array, '+g')
    # plt.plot(X2_2, y_DEM_atmo, '-b')
    # plt.xlabel('reunion island elevation (m)')
    # plt.ylabel('GACOS LOS changes (rad)')
    # plt.ylim([-15, 15])
    # plt.text(2000, 12, 'r value: {}'.format(np.round(coef_DEM_atmo[2], 2)), ha='left', va='top', fontsize=fontsize)
    # plt.savefig('/Volumes/Pierre_2TO/Stage_M1_Albino/Data/Piton_de_la_Fournaise/GACOS_LiCSAR_geotiff/144A_11114_000304/'
    #             '20191118_20220307/Figures/SVG/plot/DEM_atmo/DEM_atmo_{}.svg'.format(date), format='svg', dpi=100)

    # plt.figure(9, figsize=figsize)
    # coef_DEM_inter, y_DEM_inter, X2_3 = st.stats(DEM_array, inter_array)
    # plt.plot(DEM_array, inter_array, '*r')
    # plt.plot(X2_3, y_DEM_inter, '-b')
    # plt.xlabel('reunion island elevation (m)')
    # plt.ylabel('RAW Interferogramm LOS changes (rad)')
    # plt.ylim([-30, 30])
    # plt.text(2000, 25, 'r value: {}'.format(np.round(coef_DEM_inter[2], 2)), ha='left', va='top', fontsize=fontsize)
    # plt.savefig('/Volumes/Pierre_2TO/Stage_M1_Albino/Data/Piton_de_la_Fournaise/GACOS_LiCSAR_geotiff/144A_11114_000304/'
    #             '20191118_20220307/Figures/SVG/plot/DEM_inter/DEM_inter_{}.svg'.format(date), format='svg', dpi=100)

    # plt.figure(10, figsize=figsize)
    # coef_DEM_corr, y_DEM_corr, X2_4 = st.stats(DEM_array, int_corr)
    # plt.plot(DEM_array, int_corr, '.k')
    # plt.plot(X2_4, y_DEM_corr, '-b')
    # plt.xlabel('reunion island elevation (m)')
    # plt.ylabel('Corrected Interferogramm LOS changes (rad)')
    # plt.ylim([-30, 30])
    # plt.text(2000, 25, 'r value: {}'.format(np.round(coef_DEM_corr[2], 2)), ha='left', va='top', fontsize=fontsize)
    # plt.savefig('/Volumes/Pierre_2TO/Stage_M1_Albino/Data/Piton_de_la_Fournaise/GACOS_LiCSAR_geotiff/144A_11114_000304/'
    #             '20191118_20220307/Figures/SVG/plot/DEM_corr/DEM_corr_{}.svg'.format(date), format='svg', dpi=100)

    # plt.close('all')
    print("job {} done in --- %s seconds ---".format(date) % (time.time() - start_time))

textfile.close()
csvfile.close()
