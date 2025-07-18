import xarray as xr
import numpy as np
import pandas as pd
import re
import os

if __name__ == '__main__':

    for dirr in ["july_scans/s01"]:
        directory = "/users/xnb22215/magnetosphere-T-models/data/" + dirr

        filenames = []
        pattern = re.compile(r"output_.*\.nc")
        for fname in os.listdir(directory):
            if pattern.match(fname):
                filenames.append(str(directory + "/" + fname))

        print(filenames)
        times = []
        with open(directory + "/input_data.lst","r") as omni_data:
            for line in omni_data.readlines():
                data = line.split()
                times.append(data[0] + '    ' + data[1] + '    ' + data[2]) 

        field_strength = []
        for files in filenames:
            ds = xr.open_dataset(files,chunks="auto")
            magB = np.sqrt(ds.bx**2 + ds.by**2 + ds.bz**2).where(ds.x <= -2)
            field_strength.append(float(magB.sum()))
            

        header = '''Format:
           1) Year
           2) Day
           3) Min
           4) SUM(|B|) for x<=-2
           \n'''
        with open( directory + "/total_nightside_field.lst", "w") as f:
            f.write(header)
            for time,field in zip(times,field_strength):
                f.write( times + '   ' )
                f.write( f"{field:.4f}" + '\n' )
