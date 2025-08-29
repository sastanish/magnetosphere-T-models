import numpy as np
import json

## Get the correct storm id's from storms.json
with open("../../storms.json", "r") as file:
    storm_dict = json.load(file)
    storm_ids = list(storm_dict.keys())

infiles = [f'./data/omni_data_{storm_dict[id]["name"]}_filled_gaps_with_tilt.lst' for id in storm_ids]
outfiles = [f'./data/{storm_dict[id]["name"]}_TS05_parameters.lst' for id in storm_ids]

for infile,outfile in zip(infiles,outfiles):

    indata = np.genfromtxt(infile,dtype=None)

    bx = np.zeros( len(indata) )
    by = np.zeros( len(indata) )
    bz = np.zeros( len(indata) )
    vx = np.zeros( len(indata) )
    vy = np.zeros( len(indata) )
    vz = np.zeros( len(indata) )
    rho = np.zeros( len(indata) )
    symh = np.zeros( len(indata) )
    pydn = np.zeros( len(indata) )

    for i,line in enumerate(indata):
        bx[i] = line[4]
        by[i] = line[5]
        bz[i] = line[6]
        vx[i] = line[7]
        vy[i] = line[8] + 29.78
        vz[i] = line[9]
        rho[i] = line[10]
        symh[i] = line[14]
        pydn[i] = line[18]

    # Formula taken from TS05
    # BIndex = np.sqrt(rho/5) * (np.sqrt(vx**2 + vy**2 + vz**2)/400)**(5/2) * (by/5) * np.sin(np.abs(np.atan(by/bz)/2))**6
    theta = np.arctan2(by,bz)
    NIndex = 10**(-4) * np.sqrt(vx**2 + vy**2 + vz**2)**(4/3) * np.sqrt(by**2+bz**2)**(2/3) * np.power(np.abs(np.sin(theta/2)),8/3)
    SymHc = 0.8*symh-13*np.sqrt(pydn)

    F_fmt =  ["{: >4d}", "{: >4d}", "{: >3d}", "{: >3d}", "{: > 8.2f}", "{: > 8.2f}", "{: > 8.2f}", "{: > 8.1f}", "{: > 8.1f}", "{: > 8.1f}", "{: > 7.2f}", "{: > 9.0f}", "{: > 6.2f}", "{: >6d}", "{: >6d}", "{: >3d}", "{: >3d}", "{: > 8.4f}", "{: > 7.2f}"]
    with open(outfile,"w") as out_file:
        for i,line in enumerate(indata):
            if i>=6 and i<len(indata)-3:
                # Average field
                abx = np.average(bx[i-6:i+1])
                aby = np.average(by[i-6:i+1])
                abz = np.average(bz[i-6:i+1])
                # Average Bindex
                # Abind = np.average(BIndex[i-30:i])
                # Average Nindex
                Anind = np.average(NIndex[i-6:i+1])
                # Average SymHc
                ASymHc = np.average(SymHc[i-3:i+3])

                for data,fmt in zip(line,F_fmt):
                    out_file.write(fmt.format(data))
                out_file.write("{: > 8.2f}".format(abx))
                out_file.write("{: > 8.2f}".format(aby))
                out_file.write("{: > 8.2f}".format(abz))
                # out_file.write("{: > 8.4f}".format(Abind))
                out_file.write("{: > 8.4f}".format(Anind))
                out_file.write("{: > 7.1f}".format(ASymHc))
                out_file.write("\n")
