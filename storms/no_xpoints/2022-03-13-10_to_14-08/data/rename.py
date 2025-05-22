import os

with open("names.txt","r") as namefile:
    for i,line in enumerate(namefile.readlines()):
        os.rename(line[0:-1],f"2021_data_{i+1}.nc")
