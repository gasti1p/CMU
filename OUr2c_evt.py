import numpy as np
import os, glob, sys, h5py

run_number = sys.argv[-1]


path='evt_data/evt{}/'.format(run_number)					#change path to reflect data location


evt_filenames = ["bci","dtime","ev_route","ev_tof","ev_psd","ev_es","ev_ea","ev_eea","ev_e"]


# ------------------
# A) FILE CONVERTER
# ------------------

sys.stdout.write("\nConverting raw files for run {}...".format(run_number))
sys.stdout.flush()

stop_times = []
for name in evt_filenames:
    segment = 0
    for filename in sorted(glob.glob(os.path.join(path, '{}.*'.format(name)))):	#use commands below for all ev_e.xxx files in path
        with open(path + '{}_{}.dat'.format(name,segment),'w') as fw:				#open file to write the final result
            with open(filename, 'r') as f:
                header1 = f.readline()				#getting rid of header...
                if name==evt_filenames[0]: stop_times.append(header1[-7:])
                n=8         						#we need to read every n=8 characters
                for line in f:
                    line=line.strip('\n')   				#remove end of line character
                    line=line.strip('\t')   				#remove begin of line character if it exists
                    for i in range(0,len(line),n):			#loop over each line
                        fw.write(line[i:i+n])				#write 8 characters to file
                        fw.write("\n ") 				#type end of line character \n to make columns
            f.close()
        fw.close()
        segment +=1
sys.stdout.write("\rConverting raw files for run {}... DONE".format(run_number))



# ------------------
# B) SEGMENTS INFO
# ------------------
sys.stdout.write("\nCollecting segments info for run {}...".format(run_number))
sys.stdout.flush()

durations = []
for i in xrange(1,len(np.array(stop_times)),1):
    start = int(stop_times[i-1][:2])*3600 + int(stop_times[i-1][2:4])*60 + int(stop_times[i-1][-3:])
    stop = int(stop_times[i][:2])*3600 + int(stop_times[i][2:4])*60 + int(stop_times[i][-3:])

    if stop-start>0: durations.append(stop-start)
    elif stop-start<0: durations.append((24*3600-stop)+start)
    else: durations.append(0.0)
durations.insert(0, np.mean(durations[:10]))
sys.stdout.write("\rCollecting segments info for run {}... DONE".format(run_number))




# -----------------------
# C) CREATE BINARY FILE
# -----------------------
sys.stdout.write("\nCreate binary file for run {}...".format(run_number))
sys.stdout.flush()

h5data=h5py.File("evt{}.h5".format(run_number),'w')

for name in evt_filenames:
    for filename in sorted(glob.glob(os.path.join(path, '{}_*'.format(name)))):
        dset = np.loadtxt(filename)
        leaf_name = os.path.basename(filename)[:-4]
        h5data.create_dataset("/raw/all_files/{}".format(leaf_name),data=dset,compression="lzf")
        os.remove(filename)


h5data.create_dataset("/raw/all_files/durations",data=durations,compression="lzf")

h5data.close()

sys.stdout.write("\rCreate binary file for run {}... DONE".format(run_number))




# -------------------------------
# D) CREATE BASIC ARRAYS
# -------------------------------
sys.stdout.write("\nCreate basic arrays for run {}...".format(run_number))
sys.stdout.flush()

_filenames_ = ["ev_tof","ev_psd","ev_es","ev_ea","ev_eea","ev_e"]

segments = np.array(durations)

h5data = h5py.File("evt{}.h5".format(run_number),'a')

for name in _filenames_:
    list_lenda = []
    list_ls = []
    list_monitor = []
    if name==_filenames_[0]: list_route = []
    for i in xrange(0,len(np.array(segments)),1):
        route_raw = h5data["/raw/all_files/ev_route_{}".format(i)]
        route = route_raw[...]
        if name==_filenames_[0]: list_route.append(route)

        raw_file = h5data["/raw/all_files/{}_{}".format(name,i)]
        file = raw_file[...]

        GateLENDA = np.flatnonzero(route==1.0)
        GateLS = np.flatnonzero(route==4.0)
        GateMONITOR = np.flatnonzero(route==8.0)

        list_lenda.append(file[GateLENDA])
        list_ls.append(file[GateLS])
        list_monitor.append(file[GateMONITOR])
    if name==_filenames_[0]: h5data.create_dataset("raw/route",data=[j for i in list_route for j in i],compression="lzf")
    h5data.create_dataset("raw/lenda/{}".format(name),data=[j for i in list_lenda for j in i],compression="lzf")
    h5data.create_dataset("raw/ls/{}".format(name),data=[j for i in list_ls for j in i],compression="lzf")
    h5data.create_dataset("raw/monitor/{}".format(name),data=[j for i in list_monitor for j in i],compression="lzf")




bci_tot = np.array(h5data["/raw/all_files/bci_0"][...])
for i in xrange(1,len(np.array(segments)),1):
    bci_raw = h5data["/raw/all_files/bci_{}".format(i)]
    bci_tot = bci_tot + np.array(bci_raw[...])
h5data.create_dataset("raw/bci",data=bci_tot)



dtime_tot = np.array(h5data["/raw/all_files/dtime_0"][...])
for i in xrange(1,len(np.array(segments)),1):
    dtime_raw = h5data["/raw/all_files/dtime_{}".format(i)]
    dtime_tot = dtime_tot + np.array(dtime_raw[...])
h5data.create_dataset("raw/dtime",data=dtime_tot)

h5data.close()

sys.stdout.write("\rCreate basic arrays for run {}... DONE\n\n".format(run_number))
