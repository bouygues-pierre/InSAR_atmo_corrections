import matplotlib.pyplot as plt
import os
import pandas as pd
import re
from datetime import datetime
import datetime
import numpy as np

path = "/Volumes/Pierre_2TO/Stage_M1_Albino/Data/Piton_de_la_Fournaise/ZTD_data/IGN/2021/2021_01_05"
file_name = "reun_2021_005_V2-00.iwv"

def ZTD_read_IGN(file_name, path):
    file = '/'.join([path, file_name])

    with open(file, 'r') as file:
        lines_values = [line.strip() for line in file]

    for line in range(len(lines_values)):
        lines_values[line] = lines_values[line].split()
    time = []

    for line in range(len(lines_values)):
        time.append([lines_values[line][0], lines_values[line][1].split('.')[0], lines_values[line][1].split('.')[1]])

    for i in range(len(time)):
        time[i][2] = '.'.join(['0', time[i][2]])
        time[i][2] = float(time[i][2])
        time[i][2] = round(time[i][2] * 24 * 60 * 60)
        time[i][2] = str(datetime.timedelta(seconds=time[i][2]))

    for i in range(len(lines_values)):
        lines_values[i].append(':'.join(time[i]))
        lines_values[i][7] = pd.to_datetime(datetime.datetime.strptime(lines_values[i][7], '%Y:%j:%H:%M:%S'))

    lines_values = pd.DataFrame(data=lines_values, columns=["Year", "Julian_day", "ZHD", "Wet_delay", "T(K)",
                                                            "IWV", "error", "time"], dtype=float)

    lines_values["ZTD"] = lines_values["ZHD"] * 1000 + lines_values["Wet_delay"] * 1000
    lines_values["Hours"] = pd.to_datetime(lines_values["time"]).dt.time
    lines_values["Hours"] = lines_values["Hours"].values.astype(str)

    return lines_values


