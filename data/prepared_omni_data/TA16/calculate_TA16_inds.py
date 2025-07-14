import numpy as np

dates = ['2018-08-25', '2022-03-13', '2023-03-22', '2024-03-03', '2024-08-11', '2021-11-03', '2022-10-22', '2023-04-23', '2024-03-24', '2024-10-10', '2022-01-14', '2023-02-26', '2023-11-06', '2024-05-10']
for date in dates:
    infile = f'data/{date}_filled_gaps_with_tilt.lst'
    outfile = f'data/{date}_TA16_parameters.lst'

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

    # Formula taken from TA15
    # BIndex = np.sqrt(rho/5) * (np.sqrt(vx**2 + vy**2 + vz**2)/400)**(5/2) * (by/5) * np.sin(np.abs(np.atan(by/bz)/2))**6
    NIndex = 10**(-4) * np.sqrt(vx**2 + vy**2 + vz**2)**(4/3) * np.abs(by)**(2/3) * np.power(np.abs(np.sin(np.atan2(by,bz)/2)),(8/3))
    SymHc = 0.8*symh-13*np.sqrt(pydn)

    F_fmt =  ["{: >4d}", "{: >4d}", "{: >3d}", "{: >3d}", "{: > 8.2f}", "{: > 8.2f}", "{: > 8.2f}", "{: > 8.1f}", "{: > 8.1f}", "{: > 8.1f}", "{: > 7.2f}", "{: > 9.0f}", "{: > 6.2f}", "{: >6d}", "{: >6d}", "{: >3d}", "{: >3d}", "{: > 8.4f}", "{: > 7.2f}"]
    with open(outfile,"w") as out_file:
        for i,line in enumerate(indata):
            if i>=30 and i<len(indata)-15:
                # Average field
                abx = np.average(bx[i-30:i])
                aby = np.average(by[i-30:i])
                abz = np.average(bz[i-30:i])
                # Average Bindex
                # Abind = np.average(BIndex[i-30:i])
                # Average Nindex
                Anind = np.average(NIndex[i-30:i])
                # Average SymHc
                ASymHc = np.average(SymHc[i-15:i+15])

                for data,fmt in zip(line,F_fmt):
                    out_file.write(fmt.format(data))
                out_file.write("{: > 8.2f}".format(abx))
                out_file.write("{: > 8.2f}".format(aby))
                out_file.write("{: > 8.2f}".format(abz))
                # out_file.write("{: > 8.4f}".format(Abind))
                out_file.write("{: > 8.4f}".format(Anind))
                out_file.write("{: > 7.1f}".format(ASymHc))
                out_file.write("\n")
