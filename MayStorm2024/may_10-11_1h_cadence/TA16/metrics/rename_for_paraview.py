import os
with open("file_names.txt") as nfile:
    for line in nfile.readlines():
        os.rename(line[0:-1],line[-9:-4] + "_" + line[0:13] + '.nc')
