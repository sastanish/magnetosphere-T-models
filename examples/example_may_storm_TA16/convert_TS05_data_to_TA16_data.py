import numpy as np

infile = 'may_2024_storm_with_TS05_vars.dat'
outfile = 'may_2024_storm_with_TS16_vars.dat'

indata = np.genfromtxt(infile)
outdata = np.genfromtxt(infile)

bx = indata[:,4]
by = indata[:,5]
bz = indata[:,6]
vx = indata[:,7]
vy = indata[:,8]
vz = indata[:,9]
rho = indata[:,10]
symh = indata[:,12]
pydn = indata[:,16]

# Formula taken from TA15
BIndex = np.sqrt(rho/5) * (np.sqrt(vx**2 + vy**2 + vz**2)/400)**(5/2) * (by/5) * np.sin(np.atan(by/bz)/2)**6
NIndex = 10**(-4) * np.sqrt(vx**2 + vy**2 + vz**2)**(4/3) * np.abs(by)**(2/3) * np.power(np.abs(np.sin(np.atan2(by,bz)/2)),(8/3))
SymHc = 0.8*symh-13*np.sqrt(pydn)

for i in range(30,len(indata)-15):

    # Average field
    outdata[i,4] = np.average(bx[i-30:i])
    outdata[i,5] = np.average(by[i-30:i])
    outdata[i,6] = np.average(bz[i-30:i])

    # Average Bindex
    outdata[i,17] = np.average(BIndex[i-30:i])
    # Average Nindex
    outdata[i,18] = np.average(NIndex[i-30:i])
    # Average SymHc
    outdata[i,19] = np.average(SymHc[i-15:i+15])

fmt = ["%4.0f", "%3.0f", "%3.0f", "%2.0f", "%8.2f", "%8.2f", "%8.2f", "%8.1f", "%8.1f", "%8.1f", "%7.2f", "%9.0f", "%7.1f", "%1.0f", "%1.0f", "%8.4f", "%7.2f", "%8.4f", "%8.4f", "%7.1f"]

np.savetxt(outfile,outdata[:,0:20],fmt=fmt)
