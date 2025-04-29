import os
with open("file_list.txt") as nfile:
    for line in nfile.readlines():
        os.rename(line[0:-1],line[0:13] + '.nc')
