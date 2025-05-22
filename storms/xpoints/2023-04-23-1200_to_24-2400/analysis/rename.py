import os

date = "2023-04-24"
dir_str = f"../data/"
directory = os.fsencode(dir_str)

for i,file in enumerate(os.listdir(directory)):
    filename = os.fsdecode(file)
    try:
         os.rename(dir_str + filename,dir_str + f"{date[0:7]}_data_{i+1}.nc")
         print(dir_str + filename + "-->" + dir_str + f"{i+1}.nc")
    except:
        continue
