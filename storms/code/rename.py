import os
dates = ["2022-01-14_to_16", "2023-11-06_to_06", "2022-03-13-10_to_14-08", "2024-03-03-0000_to_04-1200", "2022-10-22_to_24", "2024-03-24_to_25", "2023-02-26-1200_to_27-2400", "2024-05-10_to_11", "2023-04-23-1200_to_24-2400", "2024-08-11_to_13"]
for date in dates:
    dir_str = f"./{date}/data/"
    directory = os.fsencode(dir_str)
 
    for i,file in enumerate(os.listdir(directory)):
        filename = os.fsdecode(file)
        try:
            os.rename(dir_str + f"{i+1}.nc",dir_str + f"{date[0:7]}_data_{i+1}.nc")
#            print(dir_str + "{:04d}.nc".format(i+1) + "-->" + dir_str + f"{i+1}.nc")
        except:
            continue





