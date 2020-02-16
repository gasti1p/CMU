import numpy as np
import os, glob, sys, h5py

run_number = sys.argv[-1]


path='dsk_data/dsk{}/'.format(run_number)      #change path to reflect data location


# ------------------
# A) FILE CONVERTER
# ------------------

sys.stdout.write("\nConverting raw files for run {}...".format(run_number))
sys.stdout.flush()
num = 0
for filename in sorted(glob.glob(os.path.join(path, 'tag*'))):	#use commands below for all ev_e.xxx files in path
    with open(path+'tag_{}.dat'.format(num),'w') as fw:        #open file to write the final result
        with open(filename, 'r') as f:
          header1 = f.readline()		#getting rid of header...
          n=8         			        #we need to read every n=8 characters
          for line in f:
             line=line.strip('\n')   		#remove end of line character
             line=line.strip('\t')   		#remove begin of line character if it exists
             for i in range(0,len(line),n):	#loop over each line
                fw.write(line[i:i+n])		#write 8 characters to file
                fw.write("\n ") 		#type end of line character \n to make columns
        f.close()
    fw.close()
    num += 1
sys.stdout.write("\rConverting raw files for run {}... DONE".format(run_number))





# -----------------------
# B) CREATE BINARY FILE
# -----------------------
sys.stdout.write("\nCreate binary file for run {}...".format(run_number))
sys.stdout.flush()

tag_lenda = ["tag_8.dat","tag_9.dat","tag_10.dat","tag_11.dat"]
tag_ls = ["tag_24.dat","tag_25.dat","tag_26.dat","tag_27.dat","tag_30.dat","tag_31.dat"]
tag_monitor =  ["tag_32.dat","tag_33.dat","tag_34.dat","tag_35.dat","tag_38.dat"]
tag_general =  ["tag_63.dat","tag_64.dat"]

hist_basic = ["ea","es","psd","tof","tof_g","tof_n"]
hist_general = ["bci","dtime"]

h5data=h5py.File("dsk{}.h5".format(run_number),'w')

for filename in tag_lenda:
    dset = np.loadtxt(path+filename)
    leaf_name = hist_basic[tag_lenda.index(filename)]
    h5data.create_dataset("/raw/lenda/{}".format(leaf_name),data=dset)

for filename in tag_ls:
    dset = np.loadtxt(path+filename)
    leaf_name = hist_basic[tag_ls.index(filename)]
    h5data.create_dataset("/raw/ls/{}".format(leaf_name),data=dset)

for filename in tag_monitor:
    dset = np.loadtxt(path+filename)
    leaf_name = hist_basic[tag_monitor.index(filename)]
    h5data.create_dataset("/raw/monitor/{}".format(leaf_name),data=dset)

for filename in tag_general:
    dset = np.loadtxt(path+filename)
    leaf_name = hist_general[tag_general.index(filename)]
    h5data.create_dataset("/raw/{}".format(leaf_name),data=dset)

h5data.close()
sys.stdout.write("\rCreate binary file for run {}... DONE\n\n".format(run_number))


# -----------------------
# C) CLEAN TEMP FILES
# -----------------------
for filename in sorted(glob.glob(os.path.join(path, 'tag_*'))):	os.remove(filename)
