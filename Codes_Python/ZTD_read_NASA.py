import matplotlib.pyplot as plt
import pandas as pd
import re
from datetime import datetime
import datetime
import time


# file_name = 'reun0050.21zpd'
# path = '/Volumes/Pierre_2TO/Stage_M1_Albino/Data/Piton_de_la_Fournaise/ZTD_data/NASA/ZTD_REUN00REU'

def ZTD_read_NASA(file_name, path):
    file = '/'.join([path, file_name])
    motif = re.compile(r'^( *REUN +\d\d:|.*EPOCH)')

    with open(file, 'r') as f2:

        lines_values = [line.strip() for line in f2 if motif.match(line)]

    for i in range(0, len(lines_values)):
        lines_values[i] = lines_values[i].split()

    for i in range(1, len(lines_values)):
        lines_values[i][1] = lines_values[i][1].split(":")
        lines_values[i][1][2] = str(datetime.timedelta(seconds=float(lines_values[i][1][2])))

        lines_values[i][1] = ':'.join(lines_values[i][1])

        lines_values[i][1] = datetime.datetime.strptime(lines_values[i][1], '%y:%j:%H:%M:%S')
        # lines_values[i][1] = pd.to_datetime(lines_values[i][1], utc=True).tz_convert(tz='Indian/Reunion')

    lines_values = pd.DataFrame(lines_values[1:], columns=[str(lines_values[0][0]), 'EPOCH', lines_values[0][2],
                                                           lines_values[0][3], lines_values[0][4], lines_values[0][5],
                                                           lines_values[0][6], lines_values[0][7]], dtype=float)
    lines_values.style.hide(axis='index')

    return lines_values

