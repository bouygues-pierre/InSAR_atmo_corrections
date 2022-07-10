import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
import struct
import io
import pandas as pd
import re
from datetime import datetime
import datetime
import time
from ZTD_read_NASA import ZTD_read_NASA
from ZTD_read_OVPF import ZTD_read_OVPF
from ZTD_read_IGN import ZTD_read_IGN
from ZTD_read_OVPF import ZTD_read_OVPF
import re
import concurrent.futures
import warnings


start_time = time.time()
warnings.filterwarnings('ignore', '.*Could not.*', )

# ----------------------------------------------------------------------------------------------------------------------

#                                       NASA IGN OVPF plot of one day

# ----------------------------------------------------------------------------------------------------------------------
''''
# OVPF file
file_name_OVPF = 'met_reun.21005'
path_OVPF = '/Volumes/Pierre_2TO/Stage_M1_Albino/Data/Piton_de_la_Fournaise/ZTD_data/OVPF/2021'

# NASA file
file_name_NASA = 'reun0050.21zpd'
path_NASA = '/Volumes/Pierre_2TO/Stage_M1_Albino/Data/Piton_de_la_Fournaise/ZTD_data/NASA/ZTD_REUN00REU'

# IGN file
file_name_IGN = "reun_2021_005_V2-00.iwv"
path_IGN = "/Volumes/Pierre_2TO/Stage_M1_Albino/Data/Piton_de_la_Fournaise/ZTD_data/IGN/2021/2021_01_05"

values_OVPF = ZTD_read_OVPF(file_name_OVPF, path_OVPF)
values_NASA = ZTD_read_NASA(file_name_NASA, path_NASA)
values_IGN = ZTD_read_IGN(file_name_IGN, path_IGN)
values_IGN.insert(0, "time2", values_NASA["EPOCH"])
values_IGN["ZHD"] *= 1000
values_IGN["Wet_delay"] *= 1000
values_IGN["ZTD"] = values_IGN["ZHD"] + values_IGN["Wet_delay"]


values_OVPF_resamp = values_OVPF.resample('1min', on="time").mean()
values_NASA_resamp = values_NASA.resample('1min', on="EPOCH").mean()
values_IGN_resamp = values_IGN.resample('1min', on="time").mean()


values_NASA_interp = values_NASA_resamp.interpolate().reset_index()
values_OVPF_interp = values_OVPF_resamp.interpolate().reset_index()
values_IGN_interp = values_IGN_resamp.interpolate().reset_index()

fig, axs = plt.subplots(figsize=(14, 9))
values_NASA.plot.scatter(x="EPOCH", y="TROTOT", c="blue", s=5, label="NASA", ax=axs)
values_OVPF.plot.scatter(x="time", y="Total_Zen", ax=axs, c="red", s=5, label="OVPF")
values_IGN.plot.scatter(x="time", y="ZTD", ax=axs, c="green", s=5, label="IGN")
axs.set_title('Zenithal total delay (ZTD) of {}'.format(str(values_NASA["EPOCH"].iat[0])))
axs.set_ylabel("ZTD")
axs.set_xlabel("Time")
axs.axvline(pd.to_datetime('2021-01-05 14:53:00'),
            label="satellite acquisition : \n NASA : {} \n OVPF : {} \n IGN : {}"
            .format(values_NASA_interp["TROTOT"][values_NASA_interp.eq('2021-01-05 14:53:00').any(1)],
                    values_OVPF_interp["Total_Zen"][values_OVPF_interp.eq('2021-01-05 14:53:00').any(1)],
                    values_IGN_interp["ZTD"][values_IGN_interp.eq('2021-01-05 14:53:00').any(1)], index=False))
plt.legend()
plt.savefig(
    '/Volumes/Pierre_2TO/Stage_M1_Albino/Data/Piton_de_la_Fournaise/ZTD_data/figures/NASA_IGN_OVPF_one_day_compar.svg'
    , dpi=100, format='svg'
)
plt.show()
#axs.text(str(values_IGN_interp["ZTD"][values_IGN_interp.eq('2021-01-05 18:53:00').any(1)]))
'''''


# ----------------------------------------------------------------------------------------------------------------------

#                                      IGN time serie of 2021 for each station

# ----------------------------------------------------------------------------------------------------------------------
'''''
path_IGN = "/Volumes/Pierre_2TO//Stage_M1_Albino/Data/Piton_de_la_Fournaise/ZTD_data/IGN/2021" # D:

liste_directory = os.listdir(path_IGN)[1:]
lenght = []

# trouver la liste de station avec le nombre de station maximal parmi toutes les données
for dir in liste_directory:
    os.chdir('/'.join([path_IGN, dir]))
    liste_file = os.listdir(os.getcwd())
    lenght.append(len(liste_file))

station = os.listdir(os.path.join(path_IGN, str(liste_directory[np.max(lenght)])))

for i in range(len(station)):
    station[i] = station[i].split('_')[0]

# liste des noms des dossiers de chaque jours de l'année 2021
liste_directory = os.listdir(path_IGN)[1:]

# récupération des données de chaque jours de l'année 2021 pour chaque station

for i in range(len(station)):
    print("working on {} station".format(station[i]))
    motif = re.compile(r'^{}'.format(station[i]))

    data = []
    i2 = 0
    for i2 in range(len(liste_directory)):
        liste_file2 = []
        os.chdir('/'.join([path_IGN, liste_directory[i2]]))
        liste_file = os.listdir(os.getcwd())
        for i3 in range(len(liste_file)):
            liste_file2.append(liste_file[i3].split("_")[0])
        if str(station[i]) in liste_file2:
            liste_date = []
            [liste_date.append(liste_file[j]) for j in range(len(liste_file)) if motif.match(liste_file[j])]
            [data.append(ZTD_read_IGN(str(liste_date[h]), '/'.join([path_IGN, str(liste_directory[i2])]))) for h in
             range(len(liste_date))]

    time_serie = pd.concat(data)
    time_serie["ZTD"] = time_serie["ZHD"] * 1000 + time_serie["Wet_delay"] * 1000
    fig, ax = plt.subplots(figsize=(14, 9))
    time_serie.plot(x="Julian_day", y="ZTD", kind="line", ax=ax)
    ax.set_title("2021 ZTD time serie for {} station".format(station[i]))
    ax.set_xlabel("time in day")
    ax.set_ylabel("Zenital total delay (mm)")
    plt.savefig(
        '/Volumes/Pierre_2TO/Stage_M1_Albino/Data/Piton_de_la_Fournaise/ZTD_data/figures/IGN_station_time_serie/svg/{}.svg'
        .format(station[i])
        , format='svg'
    )
    plt.savefig(
        '/Volumes/Pierre_2TO/Stage_M1_Albino/Data/Piton_de_la_Fournaise/ZTD_data/figures/IGN_station_time_serie/png/{}.png'
        .format(station[i])
        , format='png'
    )
    plt.close('all')
'''
# ----------------------------------------------------------------------------------------------------------------------

#                                    offset between IGN and OVPF over a year

# ----------------------------------------------------------------------------------------------------------------------
'''''
path_IGN = "/Volumes/Pierre_2TO/Stage_M1_Albino/Data/Piton_de_la_Fournaise/ZTD_data/IGN/2021"
path_OVPF = "/Volumes/Pierre_2TO/Stage_M1_Albino/Data/Piton_de_la_Fournaise/ZTD_data/OVPF/2021"

liste_directory = os.listdir(path_IGN)[1:]
lenght = []
station = 'reun'
motif = re.compile(r'^{}'.format(station))

data_IGN = []
data_OVPF = []
i2 = 0


print('load IGN data')
for i2 in range(len(liste_directory)):
    liste_file2 = []
    os.chdir('/'.join([path_IGN, liste_directory[i2]]))
    liste_file = os.listdir(os.getcwd())
    for i3 in range(len(liste_file)):
        liste_file2.append(liste_file[i3].split("_")[0])
    if str(station) in liste_file2:
        liste_date = []
        [liste_date.append(liste_file[j]) for j in range(len(liste_file)) if motif.match(liste_file[j])]
        [data_IGN.append(ZTD_read_IGN(str(liste_date[h]), '/'.join([path_IGN, str(liste_directory[i2])]))) for h in
         range(len(liste_date))]

time_serie_IGN = pd.concat(data_IGN)
#time_serie_IGN["ZTD"] = time_serie_IGN["ZHD"] * 1000 + time_serie_IGN["Wet_delay"] * 1000


print("load OVPF data")
for i in range(1, 365):
    if i < 10:
        file_OVPF = 'met_' + station + '.21'+'00{}'.format(i)
    if 10 < i < 100:
        file_OVPF = 'met_' + station + '.21' + '0{}'.format(i)
    if i > 100:
        file_OVPF = 'met_' + station + '.21' + '{}'.format(i)

    data_OVPF.append(ZTD_read_OVPF(file_OVPF, path_OVPF))


time_serie_OVPF = pd.concat(data_OVPF)
time_serie_OVPF = time_serie_OVPF.resample('5min', on="time").mean()
time_serie_OVPF = time_serie_OVPF.interpolate().reset_index()

#time_serie_IGN["time"] = time_serie_OVPF["time"]

ax = plt.subplot()
time_serie_OVPF.plot(x="time", y="Total_Zen", ax=ax)
time_serie_IGN.plot(x="Julian_day", y="ZTD", ax=ax)
plt.plot()

#diff = pd.DataFrame()
#diff = diff.assign(diff=time_serie_OVPF["Total_Zen"] - time_serie_IGN["ZTD"])
#print(diff)
#diff = pd.DataFrame(data=[time_serie_IGN["time"], time_serie_OVPF["Total_Zen"], time_serie_IGN["ZTD"]],
                    #columns=["time", "ZTD_OVPF", "ZTD_IGN"])
#diff["time"] = time_serie_OVPF["time"]
#diff["ZTD"] = time_serie_OVPF["Total_Zen"] - time_serie_IGN["ZTD"]



#time_serie_OVPF.plot(x="time", y="Total_Zen")

'''''
# ----------------------------------------------------------------------------------------------------------------------

#                                    difference between IGN and OVPF for 1 year

# ----------------------------------------------------------------------------------------------------------------------

path_IGN = "/Volumes/Pierre_2TO/Stage_M1_Albino/Data/Piton_de_la_Fournaise/ZTD_data/IGN/2021"  # D:/
path_OVPF = "/Volumes/Pierre_2TO/Stage_M1_Albino/Data/Piton_de_la_Fournaise/ZTD_data/OVPF/2021"  # D:/

liste_directory = os.listdir(path_IGN)[1:]
lenght = []
station = 'pier'
motif = re.compile(f'^{station}')

data_IGN = []
data_OVPF = []
i2 = 0


print('load IGN data')
for i2 in range(len(liste_directory)):
    os.chdir('/'.join([path_IGN, liste_directory[i2]]))
    liste_file = os.listdir(os.getcwd())
    liste_file2 = [liste_file[i3].split("_")[0] for i3 in range(len(liste_file))]
    if station in liste_file2:
        liste_date = []
        [liste_date.append(liste_file[j]) for j in range(len(liste_file)) if motif.match(liste_file[j])]
        [data_IGN.append(ZTD_read_IGN(str(liste_date[h]), '/'.join([path_IGN, str(liste_directory[i2])]))) for h in
         range(len(liste_date))]

time_serie_IGN = pd.concat(data_IGN)

print("load OVPF data")
for i in range(1, 365):
    if i < 10:
        file_OVPF = f'met_{station}.21' + f'00{i}'
    if 10 < i < 100:
        file_OVPF = f'met_{station}.21' + f'0{i}'
    if i > 100:
        file_OVPF = f'met_{station}.21' + f'{i}'
    if file_OVPF in os.listdir('/Volumes/Pierre_2TO/Stage_M1_Albino/Data/Piton_de_la_Fournaise/ZTD_data/OVPF/2021'):
        data_OVPF.append(ZTD_read_OVPF(file_OVPF, path_OVPF))

# faire une condition sur l'offset de 30 min des données et mettre des nan à la place pour eviter les sauts de valeurs.

time_serie_OVPF = pd.concat(data_OVPF)
time_serie_OVPF = time_serie_OVPF.resample('5min', on="time").mean()
time_serie_OVPF = time_serie_OVPF.interpolate().reset_index()
time_serie_OVPF["Hours"] = time_serie_OVPF["time"].dt.time
time_serie_OVPF["Hours"] = time_serie_OVPF["Hours"].values.astype(str)


# acquisition des donnnées ascendantes et descendantes
data_IGN_acqui_desc = time_serie_IGN[time_serie_IGN["Hours"] == "01:45:00"]

data_OVPF_acqui_desc = time_serie_OVPF[time_serie_OVPF["Hours"] == "01:45:00"]

data_IGN_acqui_asc = time_serie_IGN[time_serie_IGN["Hours"] == "14:50:00"]

data_OVPF_acqui_asc = time_serie_OVPF[time_serie_OVPF["Hours"] == "14:50:00"]


ax = data_IGN_acqui_desc.plot(x="time", y="ZTD", c='red', label="IGN")
data_OVPF_acqui_desc.plot(x="time", y="Total_Zen", ax=ax, c='blue', label="OVPF")
ax.set_title("Time serie of ZTD at descendant acquisition hour for IGN and OVPF of reun station")
ax.set_ylabel("ZTD (mm)")
ax.set_xlabel("2021 year")

ax = data_IGN_acqui_asc.plot(x="time", y="ZTD", c='red', label="IGN")
data_OVPF_acqui_asc.plot(x="time", y="Total_Zen", ax=ax, c='blue', label="OVPF")
ax.set_title("Time serie of ZTD at ascendant acquisition hour for IGN and OVPF of reun station")
ax.set_ylabel("ZTD (mm)")
ax.set_xlabel("2021 year")

time_serie_IGN["ZHD"] = time_serie_IGN["ZHD"] * 1000
ax = time_serie_IGN.plot(x="time", y="ZHD", c='red', label="IGN")
time_serie_OVPF.plot(x="time", y="ZHD (mm)", ax=ax, c='blue', label="OVPF")
ax.set_title("Time serie of ZHD at ascendant acquisition hour for IGN and OVPF of reun station")
ax.set_ylabel("ZHD (mm)")
ax.set_xlabel("2021 year")

time_serie_IGN["Wet_delay"] = time_serie_IGN["Wet_delay"] * 1000
time_serie_OVPF["Wet_Zen"] = - time_serie_OVPF["Wet_Zen"]

ax = time_serie_IGN.plot(x="time", y="Wet_delay", c='red', label="IGN")
time_serie_OVPF.plot(x="time", y="Wet_Zen", ax=ax, c='blue', label="OVPF")
ax.set_title("Time serie of wet delay at ascendant acquisition hour for IGN and OVPF of reun station")
ax.set_ylabel("Wet delay (mm)")
ax.set_xlabel("2021 year")


# acquisition des données des différences de ztd entre IGN et OVPF pour une année.
indexname = data_OVPF_acqui_asc[data_OVPF_acqui_asc["Day"] == 235].index
indexname2 = data_OVPF_acqui_asc[data_OVPF_acqui_asc["Day"] == 236].index
data_OVPF_acqui_asc = data_OVPF_acqui_asc.drop(indexname)
data_OVPF_acqui_asc = data_OVPF_acqui_asc.drop(indexname2)
print(data_IGN_acqui_asc.to_string())

print(data_OVPF_acqui_asc.to_string())

# #diff = time_serie_OVPF["Total_Zen"] - time_serie_IGN["ZTD"]
# print(diff.columns, diff.index)
# print(diff)

plt.figure()
diff = []

for i in range(len(data_IGN_acqui_asc)-1):
    diff.append(data_OVPF_acqui_asc.iloc[i, 6] - data_IGN_acqui_asc.iloc[i, 8])
    plt.plot(
        data_OVPF_acqui_asc.iloc[i, 0], data_OVPF_acqui_asc.iloc[i, 6] - data_IGN_acqui_asc.iloc[i, 8], '.r'
    )
#
# for index, row, index2, row2 in data_OVPF_acqui_asc.iterrows, data_IGN_acqui_asc.iterrows:
#     plt.plot(row["Day"], row["Total_Zen"] - row2["ZTD"])
plt.title('difference between OVPF and IGN value for 2021')
plt.xlabel('time')
plt.ylabel('ZTS(OVPF) - ZTD(IGN)')

print(f"--- {time.time() - start_time} seconds ---")
plt.show()
