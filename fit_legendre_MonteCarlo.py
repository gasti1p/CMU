import numpy as np
import pylab as pl
import numpy.polynomial.legendre as leg
from scipy.stats import norm
import sys, math, re
import scipy.integrate
import xlrd
import warnings
warnings.simplefilter("ignore")




E = 3.8
R = "n2" #options: "n0","n2","ls"
Poly_order = 4
MC_trials = 500



# Read excel file
excel_path = ("../Analysis_1.xlsx")
excel = xlrd.open_workbook(excel_path)
sheet = excel.sheet_by_index(1)
angle,dsdQ,dsdQ_err = [],[],[]


if E==4.2: # 4.2MeV  - thresh = 50keV
    for i in xrange(9,14,1):
        angle.append(np.radians(sheet.cell_value(i,7)))

        if R=="n0":
            dsdQ.append(sheet.cell_value(i,24))
            dsdQ_err.append(sheet.cell_value(i,25))
        if R=="n2":
            dsdQ.append(sheet.cell_value(i,26))
            dsdQ_err.append(sheet.cell_value(i,27))
        if R=="ls":
            dsdQ.append(sheet.cell_value(i,49))
            dsdQ_err.append(sheet.cell_value(i,50))

elif E==4.1: # 4.1MeV  - thresh = 50keV
    for i in xrange(17,22,1):
        angle.append(np.radians(sheet.cell_value(i,7)))
        if R=="n0":
            dsdQ.append(sheet.cell_value(i,24))
            dsdQ_err.append(sheet.cell_value(i,25))
        if R=="n2":
            dsdQ.append(sheet.cell_value(i,26))
            dsdQ_err.append(sheet.cell_value(i,27))
        if R=="ls":
            dsdQ.append(sheet.cell_value(i,49))
            dsdQ_err.append(sheet.cell_value(i,50))


elif E==4.0:  # 4.0MeV  - thresh = 50keV
    for i in xrange(25,33,1):
        angle.append(np.radians(sheet.cell_value(i,7)))
        if R=="n0":
            dsdQ.append(sheet.cell_value(i,24))
            dsdQ_err.append(sheet.cell_value(i,25))
        if R=="n2":
            dsdQ.append(sheet.cell_value(i,26))
            dsdQ_err.append(sheet.cell_value(i,27))
        if R=="ls":
            dsdQ.append(sheet.cell_value(i,49))
            dsdQ_err.append(sheet.cell_value(i,50))


elif E==3.9: # 3.9MeV  - thresh = 40keV
    for i in xrange(37,42,1):
        angle.append(np.radians(sheet.cell_value(i,7)))
        if R=="n0":
            dsdQ.append(sheet.cell_value(i,24))
            dsdQ_err.append(sheet.cell_value(i,25))
        if R=="n2":
            dsdQ.append(sheet.cell_value(i,26))
            dsdQ_err.append(sheet.cell_value(i,27))
        if R=="ls":
            dsdQ.append(sheet.cell_value(i,49))
            dsdQ_err.append(sheet.cell_value(i,50))


elif E==3.8: # 3.8MeV  - thresh = 40keV
    for i in xrange(46,51,1):
        angle.append(np.radians(sheet.cell_value(i,7)))
        if R=="n0":
            dsdQ.append(sheet.cell_value(i,24))
            dsdQ_err.append(sheet.cell_value(i,25))
        if R=="n2":
            dsdQ.append(sheet.cell_value(i,26))
            dsdQ_err.append(sheet.cell_value(i,27))
        if R=="ls":
            dsdQ.append(sheet.cell_value(i,49))
            dsdQ_err.append(sheet.cell_value(i,50))


elif E==3.7: # 3.7MeV  - thresh = 30keV
    for i in xrange(55,60,1):
        angle.append(np.radians(sheet.cell_value(i,7)))
        if R=="n0":
            dsdQ.append(sheet.cell_value(i,24))
            dsdQ_err.append(sheet.cell_value(i,25))
        if R=="n2":
            dsdQ.append(sheet.cell_value(i,26))
            dsdQ_err.append(sheet.cell_value(i,27))
        if R=="ls":
            dsdQ.append(sheet.cell_value(i,49))
            dsdQ_err.append(sheet.cell_value(i,50))


d = 4.5
L = 510.

if R=="ls":
    d = 12.7 # 6.55*2
    L = 531.
    deltaTheta = math.atan(d/L)
    angle = np.array(angle)-deltaTheta


dsigma_exp = dsdQ
dsigma_exp_err = dsdQ_err
dtheta = 2*math.atan(d/(2*L)) #rad
weights = 1/np.power(np.array(dsigma_exp_err),2)

coeff = leg.Legendre.fit(np.cos(angle), dsigma_exp,Poly_order, w=weights, full=False)
chisquare = np.sum(  np.power( (coeff(np.cos(angle))-dsigma_exp)/dsigma_exp ,2) )

theta = np.arange(0,math.pi,dtheta)
dsigma = coeff(np.cos(theta))*dtheta

def differential(x):
    return 2*math.pi*coeff(np.cos(x))*np.sin(x) # x in radians


String_Coeff = str(coeff)[3:].translate(None, "(){}<>[]")
Initial_Coeff = []
for i in xrange(0,Poly_order+1,1):
    Initial_Coeff.append(float(String_Coeff.split()[i]))
Initial_Coeff= np.array(Initial_Coeff)


# Trapezoidal integration
manual_integ = 2*math.pi*0.5*dtheta*( coeff(np.cos(theta[0]))*np.sin(theta[0]) + coeff(np.cos(theta[-1]))*np.sin(theta[-1]) + 2*np.sum(coeff(np.cos(theta[1:-1]))*np.sin(theta[1:-1])))

# Python integral
auto_integ = scipy.integrate.romberg(differential, 0.0 , math.pi)


pl.figure(1,figsize=(12,4))

pl.subplot(121)


max_integral = max_sq = min_sq = 0
min_integral = 1E+10
CSQ = []
Integrals = []
a0,a1,a2,a3,a4 = [],[],[],[],[]
for i in range(0,MC_trials,1):

    mc_matrix = []
    #for i in xrange(0,len(dsigma_exp)):
    #    cf = np.random.normal(dsigma_exp[i],dsigma_exp_err[i]/2)
    #    mc_matrix.append(cf)
    #dsigma_mc = np.array(mc_matrix)

    for item in dsigma_exp_err:
        cf = np.random.uniform(-item,item)
        mc_matrix.append(cf)
    dsigma_mc = dsigma_exp + np.array(mc_matrix)

    coeff_mc = leg.Legendre.fit(np.cos(angle), dsigma_mc, Poly_order)
    LC_mc = [float(i) for i in re.findall('[-+]?\d*\.\d+|\d+', str(coeff_mc))]
    a0.append(LC_mc[0])
    a1.append(LC_mc[1])
    a2.append(LC_mc[2])
    #a3.append(LC_mc[3])
    #a4.append(LC_mc[4])

    #pl.plot(np.cos(theta),coeff_mc(np.cos(theta)),"b-",linewidth=0.05,alpha=0.05)
    pl.plot(np.cos(theta),coeff_mc(np.cos(theta)),"b-",linewidth=0.05,alpha=0.05)

    def differential_i(x):
        return 2*math.pi*coeff_mc(np.cos(x))*np.sin(x)

    auto_integ_i = scipy.integrate.romberg(differential_i, 0.0 , math.pi)

    Integrals.append(auto_integ_i)
    chisquare_i = np.sum(  np.power( (coeff_mc(np.cos(angle))-dsigma_exp)/dsigma_exp ,2) )
    CSQ.append(chisquare_i)
    if auto_integ_i > max_integral:
        max_integral=auto_integ_i
        max_line = coeff_mc
        max_sq = chisquare_i
    if auto_integ_i < min_integral:
        min_integral=auto_integ_i
        min_line = coeff_mc
        min_sq = chisquare_i

Err = [np.std(a0),np.std(a1),np.std(a2),np.std(a3),np.std(a4)]
LC = [float(i) for i in re.findall('[-+]?\d*\.\d+|\d+', str(coeff))]
print "\n------------- INFO CARD for (p,{}) @ {} ---------------------\n".format(R,E)
print "    Theta(deg)\tdsigma(mb/sr)"
for i in xrange(0,len(angle),1):
        print "\t{:0.1f}\t{:0.2f} +- {:0.2f}".format(np.degrees(angle[i]),dsdQ[i],dsdQ_err[i])
print "\nLegendre coefficients:"
for i in xrange(0,len(LC),1):
    try: print "\ta({}): {:0.3f} +- {:0.3f}".format(i,LC[i],3*Err[i])
    except: print "\ta({}): {:0.3f} +- {:0.3f}".format(i,LC[i],3*Err[1])
print "\nchi-square: {:0.5f}".format(chisquare)
print "\nTotal cross section (NumPy): {:0.2f} +-  {:0.2f} mb".format(auto_integ, 0.5*(max_integral-min_integral))


pl.figure(1,figsize=(12,4))

pl.subplot(121)
pl.errorbar(np.cos(angle),dsigma_exp, yerr=dsigma_exp_err,label="Experimental points",
    fmt="ro",        # marker style
    markersize=3,
    linewidth=3,     # width of plot line
    elinewidth=0.5,  # width of error bar line
    ecolor='k',      # color of error bar
    capsize=0.5,       # cap length for error bar
    capthick=0.5    # cap thickness for error bar)
    )

pl.plot(np.cos(theta),coeff(np.cos(theta)),"b-",linewidth=0.5,alpha=0.5,label="Monte-Carlo fits")
pl.plot(np.cos(theta),coeff(np.cos(theta)),"r-",label="Best Legendre fit",linewidth=0.5)
#pl.plot(np.degrees(theta),coeff(np.cos(theta)),"b-",linewidth=0.5,alpha=0.5,label="Monte-Carlo fits")
#pl.plot(np.degrees(theta),coeff(np.cos(theta)),"r-",label="Best Legendre fit",linewidth=0.5)

pl.xlabel("cos$\Theta_{cm}$",fontsize=12)
#pl.xlabel("$\Theta_{CM}$",fontsize=12)

pl.ylabel("$\dfrac{d\sigma}{d\Omega}$ (mb/sr)",fontsize=12)
#pl.legend(loc="best") #,title="$\chi^2$ = {:0.5f}".format(chisquare))
pl.ylim((0,3.5))
pl.xlim((-1,1))



(mu, sigma) = norm.fit(Integrals)
print "\nTotal cross section (Monte-Carlo): {:0.2f} +-  {:0.2f} mb".format(mu, 3*sigma)
print "\n---------------------------------------------------------------\n"

pl.subplot(122)

step = mu/80
bins = np.arange(mu-0.3*mu,mu+0.3*mu,step)
pl.hist(Integrals,bins,alpha=0.7,label="Monte-Carlo integrals")
xmin, xmax = pl.xlim()
x = np.linspace(xmin, xmax, MC_trials)
pl.plot(x, norm.pdf(x,mu,sigma)*MC_trials*step, "r--",label="Gaussian fit")
pl.ylabel("Frequency",fontsize=12)
pl.xlabel("Integral (mb)",fontsize=12)
pl.legend(loc="best") #,title="$\chi^2$ = {:0.5f}".format(chisquare))


pl.subplots_adjust(top=0.90, bottom=0.13, left=0.07, right=0.94, hspace=0.20, wspace=0.23)
pl.show()
