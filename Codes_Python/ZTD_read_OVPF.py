import matplotlib.pyplot as plt
import os
import pandas as pd
import re
from datetime import datetime
import datetime

file_name = 'met_reun.21005'
path = '/Volumes/Pierre_2TO/Stage_M1_Albino/Data/Piton_de_la_Fournaise/ZTD_data/OVPF/2021'


def ZTD_read_OVPF(file_name, path):
    file = '/'.join([path, file_name])
    file2 = '/'.join([path, file_name])
    files_ZTD = os.listdir(path)

    motif = re.compile(r'^(.*Yr| *\d)')
    motif2 = re.compile(r'^(.*estimate)')

    with open(file, 'r') as file:
        lines_values = [line.strip() for line in file if motif.match(line)]

    with open(file2, 'r') as file:
        height = [line.strip() for line in file if motif2.match(line)]

    height = re.findall(r"-?\d+\.?\d+", str(height))
        
    time = []

    for i in range(1, len(lines_values)):
        lines_values[i] = lines_values[i].split()
        time.append('-'.join((lines_values[i][0], lines_values[i][1], lines_values[i][2], lines_values[i][3],
                              lines_values[i][4][0:-1])))

    for i in range(len(time)):
        time[i] = datetime.datetime.strptime(time[i], '%Y-%j-%H-%M-%S')

    values = pd.DataFrame(data=lines_values[1:], columns=['Year', 'Day', 'Hour', 'Min', 'Sec', 'Total_Zen', 'Wet_Zen',
    'Sig_Zen', 'PW', 'Sig_PW (mm)', 'Pressure (hPa)',
    'Temperature', 'ZHD (mm)', 'Grad NS', 'Sig NS', 'Grad EW',
    'Sig EW'], dtype=float)
    values.style.hide(axis='index')
    time = pd.to_datetime(time, utc=True)
    time = pd.DataFrame(data=time, columns=['time'])
    values = pd.concat([values, time, pd.DataFrame(data=[pd.to_datetime('2020-01-05 18:53:00')], columns=['heure sat'])]
                       , axis=1)

    return values


#values = ZTD_read_OVPF(file_name, path)




