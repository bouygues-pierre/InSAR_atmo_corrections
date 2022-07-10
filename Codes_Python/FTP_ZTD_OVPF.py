from ftplib import FTP
import os
import concurrent.futures
import time

ftp_address = 'tramontane.univ-reunion.fr'
#ftp.cwd('/ovpf/gnss/zt!d(d/ovpf/2022_015')

#Open ftp connection
ftp = FTP(ftp_address)
print("logging ...")
ftp.login(user="falbino", passwd="Sejab33Ramey)")
print("connected")

year = '2021'

for day in range(1, 365, 1):
    if day < 10:
        numday = '00'+str(day)
    if 10 < day < 100:
        numday = '0'+str(day)
    if day > 100:
        numday = str(day)

    ftp.cwd('/ovpf/gnss/ztd/ovpf/'+year+'_' + numday)

    # Get All Files
    ztd_files = [filename for filename in ftp.nlst() if 'met' in filename]

    for file in ztd_files:
        local_filename = os.path.join(
            "/Volumes/Pierre_2TO/Stage_M1_Albino/Data/Piton_de_la_Fournaise/ZTD_data/OVPF/2021", file)
        print("{} checked".format(file))
        if not os.path.isfile(local_filename):
            with open(local_filename, 'wb') as f:
                print("Downloading {}...".format(file))
                ftp.retrbinary('RETR %s' % file, f.write)







'''''
ftp_address = 'tramontane.univ-reunion.fr'
year = '2021'
#ftp.cwd('/ovpf/gnss/zt!d(d/ovpf/2022_015')
# Open ftp connection


def download_process(ztd_files):  # for file in ztd_files:
    local_filename = os.path.join(
        "/Volumes/Pierre_2TO/Stage_M1_Albino/Data/Piton_de_la_Fournaise/ZTD_data/OVPF/2021", ztd_files)

    print("connecting to...")
    ftp_2 = FTP(ftp_address)
    ftp_2.login(user="falbino", passwd="Sejab33Ramey)")

    print(f'checking files {local_filename}')

    if not os.path.isfile(local_filename):
        with open(local_filename, 'wb') as f:
            print("Downloading {}...".format(ztd_files))
            ftp_2.retrbinary('RETR %s' % ztd_files, f.write)

    ftp_2.close()

def searching_ztd_files(day):
    time.sleep(1)
    print("connecting to...")
    ftp = FTP(ftp_address)
    ftp.login(user="falbino", passwd="Sejab33Ramey)")
    print("connected")

    if day < 10:
        numday = '00'+str(day)
    if 10 < day < 100:
        numday = '0'+str(day)
    if day > 100:
        numday = str(day)

    ftp.cwd('/ovpf/gnss/ztd/ovpf/' + year + '_' + numday)

    # Get All Files
    ztd_files = [filename for filename in ftp.nlst() if 'met' in filename]
    ftp.close()

    return ztd_files


with concurrent.futures.ProcessPoolExecutor(max_workers=10) as executor:
    res = executor.map(searching_ztd_files, range(1, 4, 1))

print(list(res))

#with concurrent.futures.ProcessPoolExecutor() as executor:
    #executor.map(download_process, ztd_files)
'''
