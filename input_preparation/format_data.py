'''
This python file takes an OMNI web dataset "omni_data.lst" and outputs three files:

  omni_data_with_filled_gaps.lst <-- Missing data filled via linear interpolation.
  omni_data_with_filled_gaps_with_tilt.lst <-- Dataset with geo-tilt calculated from the geopack module
  omni_data_with_filled_gaps_with_tilt_with_ta16_inds.lst <-- Dataset with ta16 indicies

To use this script, point it to the correct input file by changing 
  omni_file = YOUR_FILE.lst
in the main section below. The fortran routine, calculate_tilt.f90, will need to be
compiled as calculate_tilt.out and accessible to this script (see the main section
below).

The omni dataset should be of the format:
    
    ITEMS                      FORMAT   
     
 1 Year                          I4        
 2 Day                           I4        
 3 Hour                          I3        
 4 Minute                        I3        
 5 BX, nT (GSE, GSM)             F8.2      
 6 BY, nT (GSE)                  F8.2      
 7 BZ, nT (GSE)                  F8.2      
 8 Vx Velocity,km/s              F8.1      
 9 Vy Velocity, km/s             F8.1      
10 Vz Velocity, km/s             F8.1      
11 Proton Density, n/cc          F7.2      
12 Proton Temperature, K         F9.0      
13 Flow pressure, nPa            F6.2      
14 AE-index, nT                  I6        
15 SYM/H, nT                     I6        
'''

import numpy as np
import os
import datetime

def calc_ta16_inds(infile):

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
  theta = np.arctan2(by,bz)
  NIndex = 10**(-4) * np.sqrt(vx**2 + vy**2 + vz**2)**(4/3) * np.sqrt(by**2+bz**2)**(2/3) * np.power(np.abs(np.sin(theta/2)),8/3)
  SymHc = 0.8*symh-13*np.sqrt(pydn)

  F_fmt =  ["{: >4d}", "{: >4d}", "{: >3d}", "{: >3d}", "{: > 8.2f}", "{: > 8.2f}", "{: > 8.2f}", "{: > 8.1f}", "{: > 8.1f}", "{: > 8.1f}", "{: > 7.2f}", "{: > 9.0f}", "{: > 6.2f}", "{: >6d}", "{: >6d}", "{: >3d}", "{: >3d}", "{: > 8.4f}", "{: > 7.2f}"]

  outfile = infile.replace('.lst','') + '_with_ta16_inds.lst'
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

  return

def fill_gaps(infile):

  file_data = np.genfromtxt(infile,dtype=None)
  IMF_flag = np.zeros(file_data.shape[0],dtype=int)
  SW_flag = np.zeros(file_data.shape[0],dtype=int)

  gap_start = - 1
  s_gap_start = - 1
  for ind,line in enumerate(file_data):
    ## IMF Fill
    if line[4] > 9999 or line[5] > 9999 or line[6] > 9999: # Data gap starts
      IMF_flag[ind] = -1
      if gap_start == -1:
        gap_start = ind - 1
    if IMF_flag[ind] == 0 and gap_start + 1 > 0: # Data gap is over
      gap_end = ind
      # Do linear extrap
      for i in range(gap_start+1,gap_end):
        file_data[i][4] = file_data[gap_start][4]*( (gap_end - i)/(gap_end - gap_start) )\
                   + file_data[gap_end][4]*( (i - gap_start)/(gap_end - gap_start) )
        file_data[i][5] = file_data[gap_start][5]*( (gap_end - i)/(gap_end - gap_start) )\
                   + file_data[gap_end][5]*( (i - gap_start)/(gap_end - gap_start) )
        file_data[i][6] = file_data[gap_start][6]*( (gap_end - i)/(gap_end - gap_start) )\
                   + file_data[gap_end][6]*( (i - gap_start)/(gap_end - gap_start) )
      gap_start = -1 # reset
    ## SW Fill
    if line[7] > 9999 or line[8] > 9999 or line[9] > 9999: # Data gap starts
      SW_flag[ind] = -1
      if s_gap_start == -1:
        s_gap_start = ind - 1
    if SW_flag[ind] == 0 and s_gap_start + 1 > 0: # Data gap is over
      s_gap_end = ind
      # Do linear extrap
      for i in range(s_gap_start+1,s_gap_end):
        file_data[i][7] = file_data[s_gap_start][7]*( (s_gap_end - i)/(s_gap_end - s_gap_start) )\
                   + file_data[s_gap_end][7]*( (i - s_gap_start)/(s_gap_end - s_gap_start) )
        file_data[i][8] = file_data[s_gap_start][8]*( (s_gap_end - i)/(s_gap_end - s_gap_start) )\
                   + file_data[s_gap_end][8]*( (i - s_gap_start)/(s_gap_end - s_gap_start) )
        file_data[i][9] = file_data[s_gap_start][9]*( (s_gap_end - i)/(s_gap_end - s_gap_start) )\
                   + file_data[s_gap_end][9]*( (i - s_gap_start)/(s_gap_end - s_gap_start) )
        file_data[i][10] = file_data[s_gap_start][10]*( (s_gap_end - i)/(s_gap_end - s_gap_start) )\
                   + file_data[s_gap_end][10]*( (i - s_gap_start)/(s_gap_end - s_gap_start) )
        file_data[i][11] = file_data[s_gap_start][11]*( (s_gap_end - i)/(s_gap_end - s_gap_start) )\
                   + file_data[s_gap_end][11]*( (i - s_gap_start)/(s_gap_end - s_gap_start) )
        file_data[i][12] = file_data[s_gap_start][12]*( (s_gap_end - i)/(s_gap_end - s_gap_start) )\
                   + file_data[s_gap_end][12]*( (i - s_gap_start)/(s_gap_end - s_gap_start) )
      s_gap_start = -1 # reset

  F_fmt =  ["{: >4d}", "{: >4d}", "{: >3d}", "{: >3d}", "{: > 8.2f}", "{: > 8.2f}", "{: > 8.2f}", "{: > 8.1f}", "{: > 8.1f}", "{: > 8.1f}", "{: > 7.2f}", "{: > 9.0f}", "{: > 6.2f}", "{: >6d}", "{: >6d}"]

  outfile = infile.replace('.lst','') + '_filled_gaps.lst'

  with open(outfile,"w") as out_file:
    for i,line in enumerate(file_data):
      for data,fmt in zip(line,F_fmt):
        out_file.write(fmt.format(data))
      out_file.write("{:>3d}".format(IMF_flag[i]))
      out_file.write("{:>3d}".format(SW_flag[i]))
      out_file.write("\n")

  return outfile

if __name__=="__main__":

  omni_file = "MY_FILE_DIR/MY_OMNI_DATA.lst"
  filled_name = fill_gaps(omni_file)

  # Point this script to the compiled calculate_tilt.out executable
  os.system(f"{"./calculate_tilt.out"} {filled_name.replace(".lst","")}")
  tilt_name = filled_name.replace(".lst","") + "_with_tilt.lst"
  calc_ta16_inds(tilt_name)
