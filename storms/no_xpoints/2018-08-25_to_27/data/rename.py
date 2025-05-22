import os

with open("names.txt","r") as namefile:
    for i,line in enumerate(namefile.readlines()):
        os.rename(line[0:-1],"{:04d}.nc".format(i+1))
