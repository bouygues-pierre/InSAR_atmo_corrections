import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import matplotlib.ticker as mticker
import pandas as pd
import matplotlib.dates as mdates
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os
from matplotlib_venn import venn2
import time
import netCDF4

start_time = time.time()

# ----------------------------------------------------------------------------------------------------------------------

#                                               Data extraction

# ----------------------------------------------------------------------------------------------------------------------


Ascendant = pd.read_csv('stats_144A_11114_000304.csv')
Descendant = pd.read_csv('stats_151D_11112_000300_2.csv')

# Ascendant values

Inter_A = Ascendant.Inter
Quality_coef_A = Ascendant.Quality_coef
stdi_A = Ascendant.stdi
stdf_A = Ascendant.stdf
Quality_coef_Pit_A = Ascendant.Quality_coef_Pit
stdi_Pit_A = Ascendant.stdi_Pit
stdf_Pit_A = Ascendant.stdf_Pit
r_value_inter_atmo_A = Ascendant.r_value_inter_atmo
r_value_DEM_atmo_A = Ascendant.r_value_DEM_atmo
r_value_DEM_inter_A = Ascendant.r_value_DEM_inter
r_value_DEM_corr_A = Ascendant.r_value_DEM_corr
a_inter_atmo_A = Ascendant.a_inter_atmo
a_DEM_atmo_A = Ascendant.a_DEM_atmo
a_DEM_inter_A = Ascendant.a_DEM_inter
a_DEM_intcorr_A = Ascendant.a_DEM_intcorr

# Descendant values
Inter_D = Descendant.Inter
Quality_coef_D = Descendant.Quality_coef
stdi_D = Descendant.stdi
stdf_D = Descendant.stdf
Quality_coef_Pit_D = Descendant.Quality_coef_Pit
stdi_Pit_D = Descendant.stdi_Pit
stdf_Pit_D = Descendant.stdf_Pit
r_value_inter_atmo_D = Descendant.r_value_inter_atmo
r_value_DEM_atmo_D = Descendant.r_value_DEM_atmo
r_value_DEM_inter_D = Descendant.r_value_DEM_inter
r_value_DEM_corr_D = Descendant.r_value_DEM_corr
a_inter_atmo_D = Descendant.a_inter_atmo
a_DEM_atmo_D = Descendant.a_DEM_atmo
a_DEM_inter_D = Descendant.a_DEM_inter
a_DEM_intcorr_D = Descendant.a_DEM_intcorr

date = pd.date_range(start='2019-11-18', end='2022-03-19', freq='12.5D').strftime('%Y-%m-%d').tolist()
date2 = pd.date_range(start='2019-11-18', end='2022-03-19', freq='D').strftime('%Y-%m-%d').tolist()
date3 = pd.date_range(start='2019-11-19', end='2022-02-24', freq='11.65D').strftime('%Y-%m-%d').tolist()
date4 = pd.date_range(start='2019-11-19', end='2022-02-24', freq='D').strftime('%Y-%m-%d').tolist()
Inter_split_A = np.copy(Inter_A)
Inter_split_D = np.copy(Inter_D)
x = np.arange(0, 68)
x2 = np.arange(0, 65)
x3 = np.arange(0, 72)
x4 = np.arange(0, 68)
duration_A = np.arange(len(Inter_A))
duration_D = np.arange(len(Inter_D))

for i in range(len(Inter_A)):
    Inter_split_A[i] = Inter_A[i].split('_')
    Inter_split_A_date = pd.to_datetime(Inter_split_A[i])
    duration_A[i] = (Inter_split_A_date[1] - Inter_split_A_date[0]).days

for i in range(len(Inter_D)):
    Inter_split_D[i] = Inter_D[i].split('_')
    Inter_split_D_date = pd.to_datetime(Inter_split_D[i])
    duration_D[i] = (Inter_split_D_date[1] - Inter_split_D_date[0]).days
'''''
# ----------------------------------------------------------------------------------------------------------------------

#                                   plot of interferograms as a function of Q

# ----------------------------------------------------------------------------------------------------------------------

plt.figure()

for i in range(len(Inter_split_A)):
    ax = plt.gca()
    plt.plot([Inter_split_A[i][0], Inter_split_A[i][1]], [Quality_coef_A[i], Quality_coef_A[i]])


plt.fill_between(x, 0, np.max(Quality_coef_A), color='green', alpha=0.4)
plt.fill_between(x, 0, np.min(Quality_coef_A), color='red', alpha=0.4)
axes = plt.gca()
fig = plt.gcf()
plt.ylabel('Q')
plt.title('Ascendant data')
ticks_loc = axes.get_xticks()
axes.xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
axes.xaxis.set_ticklabels([date[x] for x in x], fontsize=10, rotation=60)


plt.figure()

for i in range(len(Inter_split_D)):
    ax = plt.gca()
    plt.plot([Inter_split_D[i][0], Inter_split_D[i][1]], [Quality_coef_D[i], Quality_coef_D[i]])

plt.fill_between(x3, 0, np.max(Quality_coef_D), color='green', alpha=0.4)
plt.fill_between(x3, 0, np.min(Quality_coef_D), color='red', alpha=0.4)
plt.ylabel('Q')
plt.title('Descendant data')
axes = plt.gca()
fig = plt.gcf()
ticks_loc = axes.get_xticks()
axes.xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
axes.xaxis.set_ticklabels([date3[x3] for x3 in x3], fontsize=10, rotation=60)


# ----------------------------------------------------------------------------------------------------------------------

#                             plot of duration of interferogramm as a function of time

# ----------------------------------------------------------------------------------------------------------------------

# the size of the circle is a function of stdi

plt.figure()

for i in range(len(Inter_split_A)):
    ax = plt.gca()
    plt.scatter(Inter_split_A[i][0], duration_A[i], s=np.exp(stdi_A[i])/100)

axes = plt.gca()
fig = plt.gcf()
plt.title('Ascendant data')
ticks_loc = axes.get_xticks()
axes.xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
axes.xaxis.set_ticklabels([date[x2] for x2 in x2], fontsize=10, rotation=60)
plt.ylim(0, 350)

plt.figure()

for i in range(len(Inter_split_D)):
    ax = plt.gca()
    plt.scatter(Inter_split_D[i][0], duration_D[i], s=np.exp(stdi_D[i])/100)

axes = plt.gca()
fig = plt.gcf()
plt.title('Descendant data')
ticks_loc = axes.get_xticks()
axes.xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
axes.xaxis.set_ticklabels([date4[x4] for x4 in x4], fontsize=10, rotation=60)

# ----------------------------------------------------------------------------------------------------------------------

#                        histogram of swath of the Piton compared to the all island

# ----------------------------------------------------------------------------------------------------------------------

plt.figure()

sns.histplot((Quality_coef_Pit_A, Quality_coef_A), hue_order=['Quality_coef_Pit', 'Quality_coef'], stat='percent',
             cumulative=True, fill=False, element='step')
plt.title('histogram of swath of the Piton compared to the all island (Ascendant data)')


plt.figure()

sns.histplot((Quality_coef_Pit_D, Quality_coef_D), hue_order=['Quality_coef_Pit', 'Quality_coef'], stat='percent',
             cumulative=True, fill=False, element='step')
plt.title('histogram of swath of the Piton compared to the all island (Descendant data)')

# ----------------------------------------------------------------------------------------------------------------------

#                                    plot of Q as a function of stdi

# ----------------------------------------------------------------------------------------------------------------------

plt.figure()

plt.scatter(stdi_D, Quality_coef_D, marker=".", color='blue', label='Descendant', s=50)
plt.scatter(stdi_A, Quality_coef_A, marker=".", color='red', label='Ascendant', s=50)
plt.axhline(y=np.mean(Quality_coef_A), linestyle="--", color="red")
plt.axhline(y=np.mean(Quality_coef_D), linestyle="--", color="blue")
plt.legend()
plt.title('plot of Q as a function of stdi for Ascendant and Descendant data')
plt.xlabel('initial standard deviation')
plt.ylabel('Quality coefficient')
'''
# ----------------------------------------------------------------------------------------------------------------------

#                                               Venn diagramms

# ----------------------------------------------------------------------------------------------------------------------

# For all positive Q, we want to calculate the percentage of infero where the value of the directing coefficient
# of the line of the signal according to the DEM, has decreased

Q_well_corr_A = np.arange(len(Quality_coef_A))
a_well_coor_A = np.arange(len(a_DEM_intcorr_A))
a_inter_Q_A = np.arange(len(a_DEM_intcorr_A))
print(type(Q_well_corr_A))

for i in range(len(Quality_coef_A)):
    if Quality_coef_A[i] > 0:
        Q_well_corr_A[i] = 1
    elif Quality_coef_A[i] < 0:
        Q_well_corr_A[i] = -1
    if np.abs(a_DEM_intcorr_A[i]) < np.abs(a_DEM_inter_A[i]):
        a_well_coor_A[i] = 1
    elif np.abs(a_DEM_intcorr_A[i]) > np.abs(a_DEM_inter_A[i]):
        a_well_coor_A[i] = -1
    if Quality_coef_A[i] > 0 and abs(a_DEM_intcorr_A[i]) < abs(a_DEM_inter_A[i]):
        a_inter_Q_A[i] = 1
    elif Quality_coef_A[i] < 0 and abs(a_DEM_intcorr_A[i]) > abs(a_DEM_inter_A[i]):
        a_inter_Q_A[i] = -1


percent_Q_well_corr_A = np.round((np.size(Q_well_corr_A[Q_well_corr_A == 1])/np.size(Q_well_corr_A))*100, 2)
percent_a_well_coor_A = np.round((np.size(a_well_coor_A[a_well_coor_A == 1])/np.size(a_well_coor_A))*100, 2)
percent_a_inter_Q_A = np.round((np.size(a_inter_Q_A[a_inter_Q_A == 1])/np.size(a_inter_Q_A))*100, 2)


Q_well_corr_D = np.arange(len(Quality_coef_D))
a_well_coor_D = np.arange(len(a_DEM_intcorr_D))
a_inter_Q_D = np.arange(len(a_DEM_intcorr_D))
print(type(Q_well_corr_D))

for i in range(len(Quality_coef_D)):
    if Quality_coef_D[i] > 0:
        Q_well_corr_D[i] = 1
    elif Quality_coef_D[i] < 0:
        Q_well_corr_D[i] = -1
    if np.abs(a_DEM_intcorr_D[i]) < np.abs(a_DEM_inter_D[i]):
        a_well_coor_D[i] = 1
    elif np.abs(a_DEM_intcorr_D[i]) > np.abs(a_DEM_inter_D[i]):
        a_well_coor_D[i] = -1
    if Quality_coef_D[i] > 0 and abs(a_DEM_intcorr_D[i]) < abs(a_DEM_inter_D[i]):
        a_inter_Q_D[i] = 1
    elif Quality_coef_D[i] < 0 and abs(a_DEM_intcorr_D[i]) > abs(a_DEM_inter_D[i]):
        a_inter_Q_D[i] = -1

percent_Q_well_corr_D = np.round((np.size(Q_well_corr_D[Q_well_corr_D == 1])/np.size(Q_well_corr_D))*100, 2)
percent_a_well_coor_D = np.round((np.size(a_well_coor_D[a_well_coor_D == 1])/np.size(a_well_coor_D))*100, 2)
percent_a_inter_Q_D = np.round((np.size(a_inter_Q_D[a_inter_Q_D == 1])/np.size(a_inter_Q_D))*100, 2)

plt.figure()
venn2(subsets=(percent_Q_well_corr_A, percent_a_well_coor_A, percent_a_inter_Q_A),
      set_labels=('Q > 0', 'a_cor < a_uncor'))
plt.title('Ascendant')
plt.figure()
venn2(subsets=(percent_Q_well_corr_D, percent_a_well_coor_D, percent_a_inter_Q_D),
      set_labels=('Q > 0', 'a_cor < a_uncor'))
plt.title('Descendant')


print("--- %s seconds ---" % (time.time() - start_time))

plt.show()



