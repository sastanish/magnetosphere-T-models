import numpy as np

infile = 'data/omni_5min_2022-03-13-10_to_14-08.lst'

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

with open(infile[0:-4] + "_filled_SW_and_IMF_gaps.lst","w") as out_file:
    for i,line in enumerate(file_data):
        for data,fmt in zip(line,F_fmt):
            out_file.write(fmt.format(data))
        out_file.write("{:>3d}".format(IMF_flag[i]))
        out_file.write("{:>3d}".format(SW_flag[i]))
        out_file.write("\n")
