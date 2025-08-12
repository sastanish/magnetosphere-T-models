import xarray as xr
import numpy as np
import pandas as pd

def read_dscover_data():
    ds = xr.open_mfdataset([
        "data/dscover_data/oe_m1m_dscovr_s20240509000000_e20240509235959_p20240511021629_pub.nc.gz",
        "data/dscover_data/oe_m1m_dscovr_s20240510000000_e20240510235959_p20240511034031_pub.nc.gz",
        "data/dscover_data/oe_m1m_dscovr_s20240511000000_e20240511235959_p20240512021629_pub.nc.gz",
        "data/dscover_data/oe_m1m_dscovr_s20240512000000_e20240512235959_p20240513021624_pub.nc.gz"
        ])

    bxgse = ds.bx_gse.values
    bygse = ds.by_gse.values
    bzgse = ds.bz_gse.values
    times = ds.time.values

    return (times, bxgse, bygse, bzgse)

infiles = ["omni_data_s06.lst"]

for infile in infiles:

    file_data = np.genfromtxt('data/' + infile,dtype=None)
    (times, bx, by, bz) = read_dscover_data()

    dscover_ind = 0
    for ind,line in enumerate(file_data):
        ## IMF Fill
        if line[4] > 9999 or line[5] > 9999 or line[6] > 9999: # Data gap
            found = False
            # Find the right field data
            date = pd.to_datetime(str(line[0]) + "/" + str(line[1]) + "/" + str(line[2]) + " " + str(line[3]), format="%Y/%j/%H %M")
            while found is False:
                if date == times[dscover_ind]:
                    found = True
                    line[4] = bx[dscover_ind]
                    line[5] = by[dscover_ind]
                    line[6] = bz[dscover_ind]
                dscover_ind += 1
            

    F_fmt =  ["{: >4d}", "{: >4d}", "{: >3d}", "{: >3d}", "{: > 8.2f}", "{: > 8.2f}", "{: > 8.2f}", "{: > 8.1f}", "{: > 8.1f}", "{: > 8.1f}", "{: > 7.2f}", "{: > 9.0f}", "{: > 6.2f}", "{: >6d}", "{: >6d}"]

    with open('data/' + infile,"w") as out_file:
        for i,line in enumerate(file_data):
            for data,fmt in zip(line,F_fmt):
                out_file.write(fmt.format(data))
            out_file.write("\n")

