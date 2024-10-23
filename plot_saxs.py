"""
Plot the output of saxs.py

A. Martin - andrew.martin@rmit.edu.au
Oct 2024
"""

import numpy as np
import matplotlib.pyplot as plt
from saxsPlottingTools import parameters, extract_required_parameters, gaussian, convolve_gaussian, convolve_window


#
# Set up the parameters 
#

# location of output files from saxs.py
path = "./output2"

# name of the sample
sample = "TBAB500EtOH1000_1g_100ns_ff_norf_TEST"

# min and max q bins (q pixels)
plotqmin, plotqmax = 0.25, 1.8

# min and max r bins (r pixels)
plotrmin, plotrmax  = 0, 200

# convolve with a Gaussian?
convolve_flag = True
gaussian_bin_width = 1 

# element labels
#elems = ["H", "C", "N", "B", "O"] 
#elems = ["C", "N", "B", "O"] 

#experimental data
efile = "N44Ethnl_data.txt"
edata = np.loadtxt( efile, skiprows=2, delimiter="\t")
print(edata.shape)

#
# Start plotting! 
#
fname = f"{path}/{sample}_saxs.npy"
paramfile = f"{path}/{sample}_diffraction_batch_parameter_log.txt"
intensity = np.load(fname)
sp = extract_required_parameters(paramfile)
print(intensity.shape, sp.qpoints.shape)
print("elements :", sp.elements)
elems = sp.elements

if convolve_flag == True:
    disp = convolve_gaussian(intensity, 2)
else:
    disp = intensity

iq = np.where( (sp.qpoints>plotqmin)*(sp.qpoints<plotqmax))
iqe = np.where( (edata[:,0]>plotqmin)*(edata[:,0]<plotqmax))
print( "Debug iqe", np.min(edata[iqe,1]), np.max(edata[iqe,1]))
print( "Debug iqe", len(iqe[0]))

#Plot S(q)
plt.figure()
plt.plot( edata[:,0], edata[:,1] )
plt.plot( sp.qpoints[iq], disp[iq] )
plt.xlabel(r"q ($\mathrm{\AA}^{-1}$)")
plt.ylabel("S(q)")
plt.xlim([plotqmin, plotqmax])
plt.savefig(fname[:-4]+"_saxs.png")

np.savetxt(fname[:-4]+"_data.txt", np.array([np.abs(sp.qpoints[iq]),np.abs(disp[iq])]).T)

#Plot the partial pair distributions for each element pair
fname = f"{path}/{sample}_pairdist.npy"
pd = np.load(fname)
rmax = sp.nx/sp.qmaxin #points[-1]
print( "rmax ", rmax, sp.nx, sp.qpoints[-1])
r = np.arange(pd.shape[2])*rmax/pd.shape[2]
ir = np.where( (r>plotrmin)*(r<plotrmax) )
plt.figure()
plist, labels = [], []
for ie in range(pd.shape[0]):
    for ie2 in range(pd.shape[1]):
        if ie2>ie: continue
        p, = plt.plot( r[ir], np.squeeze(pd[ie,ie2,ir])) 
        plist.append(p)
        labels.append( f"{elems[ie]} {elems[ie2]}")
plt.ylabel("no. of pairs")
plt.xlabel(r"r ($\mathrm{\AA}$)")
plt.legend(plist,labels)
plt.savefig(fname[:-4]+"_pair_distribution.png")


# Plot S(q) for each element pair
fname = f"{path}/{sample}_saxs_partial.npy"
saxs_partial = np.load(fname)
plt.figure()
plist, labels = [], []
for ie in range(pd.shape[0]):
    for ie2 in range(pd.shape[1]):
        if ie2>ie: continue
        p, = plt.plot( sp.qpoints[iq], np.squeeze(saxs_partial[ie,ie2,iq]) )
        plist.append(p)
        labels.append( f"{elems[ie]} {elems[ie2]}")
plt.ylabel("S(q)")
plt.xlabel(r"q ($\mathrm{\AA}^{-1}$)")
plt.legend(plist,labels)
plt.savefig(fname[:-4]+"_partial_saxs.png")

# theoretical box pd
xmax, ymax, zmax = 73.0, 73.0, 73.0
boxpd_theory = (np.pi/2)*xmax*ymax*zmax  - (np.pi/4)*(xmax*ymax + ymax*zmax + xmax*zmax)*r
boxpd_theory += (1.0/3.0)*(xmax + ymax + zmax)*r*r #- (1.0/8.0)*r*r*r

#plot the box pair distribution
fname = f"{path}/{sample}_boxpd.npy"
boxpd = np.load(fname) 
plt.figure()
plt.plot( r[ir], boxpd[ir])
#plt.plot( r[ir], (3/4.75)*1e1*boxpd_theory[ir])
plt.xlabel(r" ($\mathrm{\AA}$)")
plt.ylabel("np. of pairs")
plt.savefig(fname[:-4]+"_boxpd.png")




plt.draw()
plt.show()
