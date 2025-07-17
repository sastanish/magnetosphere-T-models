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

        field_strength = []
        times = []
        for i,file in enumerate(filenames):
            ds = xr.open_dataset(file,chuncks="auto")
            magB = np.sqrt(ds.bx**2 + ds.by**2 + ds.bz**2).where(ds.x <= -2)
            field_strength[i] = float(magB.sum())
            times[i] = str(ds.attrs["time"])

        header = '''Format:
           1) Time
           2) SUM(|B|) for x<=-2
           \n'''
        with open( directory + "/total_nightside_field.lst", "w") as f:
            f.write(header)
            for time,field in zip(times,field_strength):
                f.write( times + '   ' )
                f.write( f"{field:.4f}" + '\n' )
