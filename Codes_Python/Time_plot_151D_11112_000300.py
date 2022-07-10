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
import time
start_time = time.time()

# ----------------------------------------------------------------------------------------------------------------------

#                                               Data extraction

# ----------------------------------------------------------------------------------------------------------------------


csv_file = pd.read_csv('stats_151D_11112_000300.csv')

Inter = csv_file.Inter
Quality_coef = csv_file.Quality_coef
stdi = csv_file.stdi
stdf = csv_file.stdf
Quality_coef_Pit = csv_file.Quality_coef_Pit
stdi_Pit = csv_file.stdi_Pit
stdf_Pit = csv_file.stdf_Pit
r_value_inter_atmo = csv_file.r_value_inter_atmo
r_value_DEM_atmo = csv_file.r_value_DEM_atmo
r_value_DEM_inter = csv_file.r_value_DEM_inter
r_value_DEM_corr = csv_file.r_value_DEM_corr

print(max(Quality_coef))
'''''
date = pd.date_range(start='2019-09-08', end='2022-03-20', freq='6.46D').strftime('%Y-%m-%d').tolist()
date2 = pd.date_range(start='2019-09-08', end='2022-03-20', freq='D').strftime('%Y-%m-%d').tolist()
Inter_split = np.copy(Inter)
x = np.arange(0, 144)
x2 = np.arange(0, 140)
duration = np.arange(len(Inter))

for i in range(len(Inter)):
    Inter_split[i] = Inter[i].split('_')
    Inter_split_date = pd.to_datetime(Inter_split[i])
    duration[i] = (Inter_split_date[1] - Inter_split_date[0]).days

# ----------------------------------------------------------------------------------------------------------------------

#                                   plot of interferograms as a function of Q

# ----------------------------------------------------------------------------------------------------------------------

plt.figure()
color = np.arange(np.min(stdi), np.max(stdi), 0.1)


for i in range(len(Inter_split)):
    ax = plt.gca()
    #color = next(ax._get_lines.prop_cycler)['color']
    plt.plot([Inter_split[i][0], Inter_split[i][1]], [Quality_coef[i], Quality_coef[i]]) #color=color[i])


plt.fill_between(x, 0, 0.5, color='green', alpha=0.4)
plt.fill_between(x, 0, -1.25, color='red', alpha=0.4)
axes = plt.gca()
fig = plt.gcf()
ticks_loc = axes.get_xticks()
axes.xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
axes.xaxis.set_ticklabels([date[x] for x in x], fontsize=7, rotation=60)

#colors = plt.cm.Reds(np.linspace((np.min(Quality_coef)-np.min(Quality_coef)), (np.max(Quality_coef)-np.min(Quality_coef))))

# ----------------------------------------------------------------------------------------------------------------------

#                             plot of duration of interferogramm as a function of time

# ----------------------------------------------------------------------------------------------------------------------

# the size of the circle is a function of stdi
duration_150 = duration[duration < 150]
Inter_split_150 = Inter_split[duration < 150]

print(len(duration_150), Inter_split_150, Inter_split)

plt.figure()

for i in range(len(Inter_split_150)):
    ax = plt.gca()
    plt.scatter(Inter_split_150[i][0], duration_150[i], s=(stdi[i]) ** 2)  #color=color[i])

axes = plt.gca()
fig = plt.gcf()
ticks_loc = axes.get_xticks()
axes.xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
axes.xaxis.set_ticklabels([date[x2] for x2 in x2], fontsize=7, rotation=60)

# ----------------------------------------------------------------------------------------------------------------------

#                        histogram of swath of the Piton compared to the all island

# ----------------------------------------------------------------------------------------------------------------------

plt.figure()

sns.histplot((Quality_coef_Pit, Quality_coef), hue_order=['Quality_coef_Pit', 'Quality_coef'], stat='percent',
             cumulative=True, fill=False, element='step')

# ----------------------------------------------------------------------------------------------------------------------

#                                    plot of Q as a function of stdi

# ----------------------------------------------------------------------------------------------------------------------

plt.figure()

plt.scatter(stdi, Quality_coef, marker=".", color='red')
plt.title('plot of Q as a function of stdi')
plt.xlabel('initial standard deviation')
plt.ylabel('Quality coefficient')

print("--- %s seconds ---" % (time.time() - start_time))

plt.show()
'''''
