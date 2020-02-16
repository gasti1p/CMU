
# Panagiotis Gastis
# Last update: 20 Nov 2018
# This program is set to integrate the photopeak at 770keV


import numpy.polynomial.polynomial as poly
import numpy as np
import matplotlib.pyplot as pl
import sys


run_number = 4170     #Set run number
poly_order = 1          #Set polynomial degree for the fit



### CODE STARTS HERE #############################

#----------------
#For 770keV
#----------------
fit_start = 680
fit_stop  = 816
#peak1_start = 748
#peak1_stop  = 800


#----------------
#For 843keV
#----------------
#fit_start = 800
#fit_stop  = 930
#peak1_start = 815
#peak1_stop  = 875

#----------------
#For 1014keV
#----------------
#fit_start = 930
#fit_stop  = 1080
#peak1_start = 983
#peak1_stop  = 1047


#----------------
#For background
#----------------
fit_start = 1460
fit_stop  = 1556
peak1_start = 1480
peak1_stop  = 1536


peak2_start = 1580
peak2_stop  = 1590


spec_file = "../data/LaBr_data/{}.Spe".format(run_number)



# Calibration parameters -> 2nd order
a2 = 0.000030
a1 = 0.818936
a0 = 2.367848


# Read spectra from ascii files
num_of_counts_spec  = []
num_of_counts_bgr  = []
channel = []
energy = []



chnl = 0
with open(spec_file) as f:
	for line in f:
		if line.startswith(('   ', '\t')):
			num_of_counts_spec.append( float(line.split()[0]) )
			channel.append(chnl)
			energy.append(a2*chnl**2+a1*chnl+a0)
			chnl +=1
f.close()



np.array(num_of_counts_spec)
np.array(channel)
np.array(energy)
pl.plot(energy,num_of_counts_spec,label='Spectrum')
pl.xlabel("E$_{\gamma}$ (keV)",fontsize=12)
pl.ylabel("counts/bin",fontsize=12)
pl.ticklabel_format(axis='y',style='sci',scilimits=(0,3),useMathText=True)

pl.show()
sys.exit()



# Preparations to fit the background
Peak_Region = []
Peak_Counts = []
Fit_Region = []
x = []
y = []
for i in energy:
	if i>peak1_start and i<peak1_stop:
		Peak_Region.append(int (i))
		Peak_Counts.append(num_of_counts_spec[energy.index(i)])

	if i>=fit_start and i<=fit_stop:
		Fit_Region.append(i)
		x.append(i)
		y.append(num_of_counts_spec[energy.index(i)])
np.array(x)
np.array(y)



# remove 1st peak (@770keV)
index1 = []
for i in x:
	if i>peak1_start and i<peak1_stop:
		#print "%f "%i + "%i "%x.index(i) + "%f "%y[x.index(i)]
		index1.append(x.index(i))
np.array(index1)
del x[min(index1):max(index1)]
del y[min(index1):max(index1)]



# remove 2nd peak (@843keV)
try:
	index2 = []
	for i in x:
		if i>peak2_start and i<peak2_stop:
			index2.append(x.index(i))
	np.array(index2)
	del x[min(index2):max(index2)]
	del y[min(index2):max(index2)]
except: print "\nThe 2nd peak has been ignored!\n"



# fit background with polynomial function
fit_bgr = poly.Polynomial(poly.polyfit(x, y, poly_order))



# Integrate peak and remove background
print Peak_Region
Background = np.array(fit_bgr(Peak_Region))
Net_Integral = np.array(Peak_Counts)

#Peak_Integral[Peak_Integral<0.0] = 0.0 #Set all negative numbers equal zero

Total_Number_of_Counts = np.sum(Net_Integral)
Error_total = np.sqrt(Total_Number_of_Counts)

Background_Counts = np.sum(fit_bgr(Peak_Region))
Error_Bgr = np.sqrt(Background_Counts)


Integral = Total_Number_of_Counts - Background_Counts

if Integral<0: Integral=0

Error = np.sqrt( np.power(Error_total,2) + np.power(Error_Bgr,2) )




print "\n------------------------------"
print "### PHOTOPEAK ANALYSIS ###"
print "------------------------------"
print "Minimum peak energy (keV): {}".format(min(Peak_Region))
print "Maximum peak energy (keV): {}".format(max(Peak_Region))
print "Average peak energy (keV): {}".format(np.mean(Peak_Region))
print "------------------------------"
print "NET COUNTS: {}".format(Total_Number_of_Counts) + " +- " + "{}".format(Error_total)
print "BACKGROUND COUNTS: {}".format(Background_Counts) + " +- " + "{}".format(Error_Bgr)
print "------------------------------"
print "INTEGRAL: {}".format(Integral) + " +- " + "{}".format(Error)
print "------------------------------\n"



# Make plots
pl.figure(1,figsize=(15,9))


pl.plot(energy,num_of_counts_spec,label='Spectrum')
#pl.plot(energy,num_of_counts_bgr)
#pl.plot(energy,Norm_bgr)
#pl.plot(energy,Clean_Spectrum)

pl.plot(x,y,'yo',label='Background points')
pl.plot(Fit_Region, fit_bgr(Fit_Region), 'r-', label='Fit function')


pl.xlabel('Energy (keV)')
pl.ylabel('Counts')
pl.legend(title="Run number: {}".format(run_number) + "\nPolynomial degree: {}".format(poly_order),loc='best')

#pl.yscale('log')
#pl.savefig('figure.jpeg')
pl.show()
