import numpy as np
from numpy import histogram
import pylab as pl
from matplotlib.colors import LogNorm
from scipy.stats import norm
import numpy.polynomial.polynomial as poly
import sys
import warnings
import h5py
warnings.simplefilter("ignore")

path_to_data = "../data/evt/"


run_number = 4150

print2file = False

bgr_min = 1.8  #MeV
bgr_max = 4.4  #MeV

n0_min = 0.89 #MeV
n0_max = 1.25  #MeV

n2_min = 0.26  #MeV
n2_max = 0.42  #MeV

threshold = 30.0 #keV
correction = False
calibration = "B"
Buncher = 1  # 1->1600ns  or  2->800ns

intervals = 8
printIntervals = False


#Info for calibration and photopeak limits
L = 5.10            # m
c = 299792458       # m/s
photons_tof = 1E+9*L/c #nsec

poly_order = 3

if Buncher == 1: #@1600ns
    slope = -0.409838   #nsec/channel
    photo_min = 3410
    photo_max = 3450

elif Buncher == 2: #@800ns
    slope = -0.409838/1.9   #nsec/channel
    photo_min = 2810 #4127
    photo_max = 2870 #4127


# -------------------------
#     LENDA calibration
# -------------------------
if calibration=="A":
    offset_ea = -9.111629
    offset_es = -13.717808
    slope_ea = 1.015557
    slope_es = 1.002159

if calibration=="B":
    offset_ea = -7.073681
    offset_es = -5.838526
    slope_ea = 1.026253
    slope_es = 0.990459
'''
if calibration=="A":
    ea_ch_points = [67.6,479]     #channel CAL- A
    es_ch_points = [73.1,490]     #channel CAL- A
if calibration=="B":
    ea_ch_points = [65.2,472]     #channel CAL- B
    es_ch_points = [63.34,488]    #channel CAL- B
e_points = [59.54,477.34]     # [Am, Cs]keV
calibration_ea = np.vectorize(poly.Polynomial(poly.polyfit(ea_ch_points, e_points, 1)))
calibration_es = np.vectorize(poly.Polynomial(poly.polyfit(es_ch_points, e_points, 1)))
offset_ea = calibration_ea(0)
offset_es = calibration_es(0)
slope_ea = (calibration_ea(300)-offset_ea)/300
slope_es = (calibration_es(300)-offset_es)/300
'''

# -------------------------
#     Load data
# -------------------------


h5data = h5py.File(path_to_data + "evt{}.h5".format(run_number),'a')
sys.stdout.write("\nReading leafs ...")
sys.stdout.flush()
ev_tof_lenda = h5data["raw/lenda/ev_tof"][...] - 4*4096
ev_es_lenda = h5data["raw/lenda/ev_es"][...] - 3*4096
ev_ea_lenda = h5data["raw/lenda/ev_ea"][...] - 2*4096
rtime =  h5data["raw/all_files/durations"][...]
dtime = h5data["raw/dtime"][...] #  512 channels -> histogram
bci = h5data["raw/bci"][...]     # 4096 channels -> histogram
sys.stdout.write("\rReading leafs ... DONE\n")
h5data.close()


h5data = h5py.File("../data/dsk/" + "dsk{}.h5".format(run_number),'a')
sys.stdout.write("\nReading dsk leafs ...")
sys.stdout.flush()
dsk_tof_lenda = h5data["raw/lenda/tof"][...]
bci = h5data["raw/bci"][...]     # 4096 channels -> histogram
dtime = h5data["raw/dtime"][...] #  512 channels -> histogram
sys.stdout.write("\rReading dsk leafs ... DONE\n")
h5data.close()


if run_number == 4138:
    h5data = h5py.File(path_to_data + "evt4139.h5",'a')
    ev_tof_lenda_b = h5data["raw/lenda/ev_tof"][...] - 4*4096
    ev_es_lenda_b = h5data["raw/lenda/ev_es"][...] - 3*4096
    ev_ea_lenda_b = h5data["raw/lenda/ev_ea"][...] - 2*4096
    ev_tof_lenda = np.append(ev_tof_lenda,ev_tof_lenda_b)
    ev_es_lenda = np.append(ev_es_lenda,ev_es_lenda_b)
    ev_ea_lenda = np.append(ev_ea_lenda,ev_ea_lenda_b)

# -------------------------
# Calculate dead time & BCI
# -------------------------
ADC_counts = np.sum(dtime[70:96])
TAC_counts = np.sum(dtime[140:180])
CLOCK_counts = np.sum(dtime[280:310])
ADC_deadtime = ADC_counts/(10*CLOCK_counts)
TAC_deadtime = TAC_counts/(10*CLOCK_counts)
live_time = 1. - max([ADC_deadtime,TAC_deadtime])



# -------------------------
#   Apply threshold cuts
# -------------------------
ev_es_lenda = ev_es_lenda*slope_es + offset_es
ev_ea_lenda = ev_ea_lenda*slope_ea + offset_ea

if threshold>0.0:
    gate1 = np.flatnonzero(ev_es_lenda>threshold)
    gate2 = np.flatnonzero(ev_ea_lenda>threshold)
    #thresh = reduce(np.union1d, (gate1,gate2))
    thresh = reduce(np.intersect1d,(gate1,gate2))
    ev_es_lenda = ev_es_lenda.ravel()[thresh]
    ev_ea_lenda = ev_ea_lenda.ravel()[thresh]
    ev_tof_lenda = ev_tof_lenda.ravel()[thresh]


# --------------------------------------
#   Timing correction and calibration
# --------------------------------------
Gate_photopeak = np.flatnonzero((ev_tof_lenda>photo_min) & (ev_tof_lenda<photo_max))
photopeak_tof = ev_tof_lenda[Gate_photopeak]

sys.stdout.write("Collecting points for timing corrections ...")
sys.stdout.flush()
x_points = np.arange(65,900,5)
y_points = []

for i in x_points:
    Gate = np.intersect1d( np.flatnonzero((ev_ea_lenda>=i-2)&(ev_ea_lenda<=i+2)),Gate_photopeak)
    (mu, sigma) = norm.fit(ev_tof_lenda[Gate])
    y_points.append(mu)

np.array(y_points)
sys.stdout.write("\rCollecting points for timing corrections ... DONE\n")

Fit_Function = poly.Polynomial(poly.polyfit(1./np.sqrt(x_points), y_points, poly_order))
np.vectorize(Fit_Function)
offset = photons_tof - slope*np.mean(y_points)
if correction: ev_tof_lenda_corrected = ev_tof_lenda - Fit_Function(1./np.sqrt(ev_ea_lenda)) + np.mean(y_points)
else: ev_tof_lenda_corrected = ev_tof_lenda

ev_tof_lenda_calibrated = slope*ev_tof_lenda_corrected + offset




# -------------------------
#    Peak integration
# -------------------------

min = int(min(ev_tof_lenda_calibrated))
max = int(max(ev_tof_lenda_calibrated))
histo_tof_calibrated = np.histogram(ev_tof_lenda_calibrated,bins=np.arange(min,max,-slope))


if Buncher == 1: Gate_energy = np.flatnonzero((histo_tof_calibrated[1]>=150)&(histo_tof_calibrated[1]<=1300))
elif Buncher == 2: Gate_energy = np.flatnonzero((histo_tof_calibrated[1]>=150)&(histo_tof_calibrated[1]<=580))
tof_to_energy = histo_tof_calibrated[1].ravel()[Gate_energy]
energy = np.array(0.5*939.5654133*(np.power(1E+9*L/tof_to_energy,2)/np.power(c,2)))
counts_energy = histo_tof_calibrated[0].ravel()[Gate_energy]

Gate_bgr = np.flatnonzero((energy>=bgr_min)&(energy<=bgr_max))
Gate_peak_n0 = np.flatnonzero((energy>=n0_min)&(energy<=n0_max))
Gate_peak_n2 = np.flatnonzero((energy>=n2_min)&(energy<=n2_max))

counts_peak_n0 = counts_energy.ravel()[Gate_peak_n0]
counts_peak_n2 = counts_energy.ravel()[Gate_peak_n2]

counts_bgr = counts_energy.ravel()[Gate_bgr]
bgr_per_channel = np.mean(counts_bgr)

integral_peak_n0 = np.sum(counts_peak_n0) - len(counts_peak_n0)*bgr_per_channel
error_peak_n0 = np.sqrt( np.sum(counts_peak_n0) + len(counts_peak_n0)*bgr_per_channel)

integral_peak_n2 = np.sum(counts_peak_n2) - len(counts_peak_n2)*bgr_per_channel
error_peak_n2 = np.sqrt( np.sum(counts_peak_n2) + len(counts_peak_n2)*bgr_per_channel)




# -------------------------
#    Print run info
# -------------------------

print "\n\n   ========= INFO CARD FOR RUN {} =========".format(run_number)
print "\tADC dead-time: {:0.2f} %".format(100*ADC_deadtime)
print "\tTAC dead-time: {:0.2f} %\n".format(100*TAC_deadtime)
print "\tLive-time correction factor: {:.4f}".format(live_time)
print "\tRun time (sec): {}".format(int(np.sum(rtime)))
print "\tTotal BCI: {}\n".format(np.sum(bci))
print "\t------ peak integrals ------".format(run_number)
print "\t Net counts (p,n0): {} +- {}".format(int(integral_peak_n0),int(error_peak_n0))
print "\t Net counts (n,n2): {} +- {}".format(int(integral_peak_n2),int(error_peak_n2))
print "\t Background: {} per bin".format(int(bgr_per_channel))
print "   ===========================================\n"



# -------------------------
#     ANALYSE PEAK
# -------------------------
step = (n2_max-n2_min)/intervals
e_in = n2_min
e_out = n2_min+step

if printIntervals: f = open("Interv_{}.txt".format(run_number),"wb")
print "\n------ INTERVALS ------"
print " Step = {:0.2f} MeV ".format(step)
print " DE = {:0.3f} MeV ".format(step/2)
print "\n  E_n(MeV)   Wi   dWi"
for i in xrange(0,intervals):
    Gate_interval = np.flatnonzero((energy>=e_in)&(energy<=e_out))
    counts_inteval = counts_energy.ravel()[Gate_interval]
    integral_inteval = np.sum(counts_inteval) - len(counts_inteval)*bgr_per_channel
    energy_bin = np.mean(energy.ravel()[Gate_interval])

    dxi = np.sqrt(integral_inteval)/integral_inteval
    dx = error_peak_n2/integral_peak_n2
    Wi = integral_inteval/integral_peak_n2
    dWi = np.sqrt(  dxi**2 + dx**2 )

    print "   {:0.2f}     {:0.2f}   {:0.3f}  ".format(energy_bin,Wi,dWi*Wi)
    if printIntervals: f.write("{:0.2f}    {:0.2f}   {:0.3f} \n".format(energy_bin,Wi,dWi*Wi))
    e_in+=step
    e_out+=step
print "   ------  ------"




# ----------------------------
#       Print to file
# ---------------------------
histo_tof = np.histogram(ev_tof_lenda,bins=np.arange(0,4096,1))
if print2file:
    f = open("../TV_spectra/{}.dat".format(run_number),"wb")
    for i in xrange(0,4095,1):
        f.write("{}\n".format(histo_tof[0][i]))
    f.close()



# -------------------------
#    Plotting
# -------------------------

pl.figure(1,figsize=(16,7))

pl.subplot(231)
binx2=np.arange(0,4000,4)
biny2=np.arange(0,4000,2)
pl.hist2d(ev_ea_lenda,ev_tof_lenda, bins=[binx2,biny2],norm=LogNorm(),label='Raw spectrum')
pl.plot(x_points,y_points,'ro',label='Timing correction points')
pl.xlabel("PMT energy (channel)")
pl.ylabel("ToF channel")
pl.legend(loc='best', title="Raw spectrum")


pl.subplot(232)
pl.plot(1./np.sqrt(x_points),y_points,"bo",label='Timing correction points')
pl.plot(1./np.sqrt(x_points),Fit_Function(1./np.sqrt(x_points)),"r-",label='Fit (polynomial order: {})'.format(poly_order))
pl.xlabel("1/sqrt(E)")
pl.ylabel("ToF channel")
pl.legend(loc='best')
pl.title("RUN {}".format(run_number))


pl.subplot(233)
binx2=np.arange(0,4000,4)
biny2=np.arange(0,4000,2)
pl.hist2d(ev_ea_lenda,ev_tof_lenda_corrected, bins=[binx2,biny2],norm=LogNorm())
pl.xlabel("PMT energy (channel)")
pl.ylabel("ToF channel")
pl.legend(loc='best',title="Corrected raw spectrum")





pl.subplot(234)
histo_tof = np.histogram(ev_tof_lenda,bins=np.arange(0,4096,1))
pl.plot(histo_tof[1][:-1],histo_tof[0],linewidth=0.5,label='Raw TOF spectrum - evt')
pl.plot(np.arange(0,4096,1),dsk_tof_lenda,linewidth=0.5,label='Raw TOF spectrum - dsk')

pl.legend(loc='best',title='Threshhold = {} keV'.format(threshold))
pl.xlabel("channel")
pl.ylabel("Counts/Channel")
#pl.xlim([0,1400])



pl.subplot(235)
pl.plot(histo_tof_calibrated[1][:-1],histo_tof_calibrated[0],linewidth=0.5,label='TOF calibrated/corrected')
pl.legend(loc='best', title='Threshhold = {} keV'.format(threshold))
pl.xlabel("TOF (nsec)")
pl.ylabel("Counts/bin")
pl.xlim([0,1400])



pl.subplot(236)
pl.plot(energy,counts_energy,linewidth=0.5,label='Neutron energy spectrum')
pl.plot(energy.ravel()[Gate_bgr],counts_bgr,label='Background region')
pl.plot(energy.ravel()[Gate_peak_n0],counts_peak_n0,label='(p,n0) region')
pl.plot(energy.ravel()[Gate_peak_n2],counts_peak_n2,label='(p,n2) region')
pl.legend(loc='best')
pl.xlabel("Energy (MeV)")
pl.ylabel("Counts/bin")
#pl.xlim([0,1400])

pl.subplots_adjust(top=0.95, bottom=0.08, left=0.06, right=0.98, hspace=0.22, wspace=0.21)


#TOF = histo_tof_calibrated[1][:-1]
#Counts_TOF = histo_tof_calibrated[0]
#f = open("neutron.dat","wb")
#f.write("TOF(nsec)   Counts_TOF    En(MeV)    Counts_En\n")
#for i in xrange(0,len(TOF),1):
#    try: f.write("{}     {}      {}      {}\n".format(TOF[i],Counts_TOF[i],energy[i],counts_energy[i]) )
#    except: f.write("{}     {}\n".format(TOF[i],Counts_TOF[i]) )
#f.close()


pl.show()
