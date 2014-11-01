#!/usr/bin/python
# Produces the combined emission and spin-down plots, and plots of the individual observations

import argparse
import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt
import Vfunctions as Vf
import glob
import os

# Read command line arguments
parser = argparse.ArgumentParser(description='Produces plots for variability studies')
parser.add_argument('-p','--pulsar', help='Pulsar name', required=True)
parser.add_argument('-c','--combined', help='Produces combined emission and spindown plot', action='store_true')
args = parser.parse_args()

print 'Read arguments'
pulsar = args.pulsar
print "PULSAR IS",pulsar
interval = 1
outfiles_ext = np.sort(np.array(glob.glob('./{0}/zoomed*.txt' .format(pulsar))))
outfiles = []
datfiles = []
for i in range(outfiles_ext.shape[0]):
    outfiles.append(os.path.splitext(outfiles_ext[i])[0])
    datfiles.append(outfiles[i] + '.dat')

mjd = np.loadtxt('./{0}/mjd.txt'.format(pulsar))
mjdremoved = np.loadtxt('./{0}/mjdremoved.txt'.format(pulsar))
mjdinfer = np.arange(mjd[0],mjd[-1],interval)

for i in range(len(outfiles)):
    
    linferredarray = np.loadtxt('{0}_linferred_array.dat' .format(outfiles[i]))
    inferredarray = np.loadtxt('{0}_inferred_array.dat' .format(outfiles[i]))
    inferredvar = np.loadtxt('{0}_inferred_var.dat' .format(outfiles[i]))
    readbins = np.loadtxt('{0}'.format(datfiles[i]))

    leftbin = readbins[0]
    rightbin = readbins[1]
    allbins = readbins[2]
    bins = rightbin-leftbin
# yaxis = np.linspace(leftbin/allbins, rightbin/allbins, bins)
    if i ==0:
        yaxis = np.linspace(0, bins/allbins, bins)
    else:
        yaxis = np.linspace(leftbin/allbins, rightbin/allbins, bins)
    maxdifference = np.amax(inferredarray)
    mindifference = np.amin(inferredarray)
    limitdifference = np.max((maxdifference, np.abs(mindifference)))

    if i == 0:
        Vf.makemap(inferredarray, -limitdifference, limitdifference, mjdinfer, yaxis, mjd, 'MJD', 'Pulse Phase', pulsar, '{0}_inferreddata.png'.format(outfiles[i]), peakline=allbins/4-leftbin)
        Vf.makemap(linferredarray, -np.log10(limitdifference), np.log10(limitdifference), mjdinfer, yaxis, mjd, 'MJD', 'Pulse Phase', pulsar, '{0}_linferreddata.png'.format(outfiles[i]), peakline=allbins/4-leftbin)
        Vf.makemap(inferredvar, 0 , np.amax(inferredvar), mjdinfer, yaxis, mjd, 'MJD', 'Pulse Phase', pulsar, '{0}_inferredvariance.png'.format(outfiles[i]), peakline=allbins/4-leftbin)

    else:
        Vf.makemap(inferredarray, -limitdifference, limitdifference, mjdinfer, yaxis, mjd, 'MJD', 'Pulse Phase', pulsar, '{0}_inferreddata.png'.format(outfiles[i]))
        Vf.makemap(linferredarray, -np.log10(limitdifference), np.log10(limitdifference), mjdinfer, yaxis, mjd, 'MJD', 'Pulse Phase', pulsar, '{0}_linferreddata.png'.format(outfiles[i]))
        Vf.makemap(inferredvar, 0 , np.amax(inferredvar), mjdinfer, yaxis, mjd, 'MJD', 'Pulse Phase', pulsar, '{0}_inferredvariance.png'.format(outfiles[i]))
        


if (args.combined):

    nudot = np.loadtxt('{0}/{0}_nudot.dat' .format(pulsar))
    Llim2 = np.loadtxt('{0}/{0}_Llim2.dat' .format(pulsar))
    Ulim2 = np.loadtxt('{0}/{0}_Ulim2.dat' .format(pulsar))
    mjdinfer_spindown = np.loadtxt('./{0}/{0}_mjdinfer_spindown.dat'.format(pulsar))
    power = int((-1)*np.floor(math.log10(abs(np.mean(nudot))))) 
    nudot = nudot*10**power
    Llim2 = Llim2*10**power
    Ulim2 = Ulim2*10**power

    fig=plt.figure()
    fig.set_size_inches(16,10)
    ax=fig.add_subplot(2,1,1)
    ax.xaxis.set_visible(False)

    Vf.makemap(linferredarray, limitdifference, limitdifference, mjdinfer, yaxis, mjd, 'MJD', 'Pulse phase', pulsar, '.not_saved.png', combined=True)

    ax=fig.add_subplot(2,1,2)

    ax.grid()

    plt.plot(mjdinfer_spindown, nudot)
    plt.fill_between(mjdinfer_spindown, Llim2, Ulim2, color = 'b', alpha = 0.2)
    plt.xlim(mjdinfer[0],mjdinfer[-1])
    plt.xlabel('MJD')
    plt.ylabel(r'$\mathrm{{\dot{{\nu}}}}$ ($\mathrm{{10^{{-{0}}} s^{{-2}}}}$)'.format(power) ,fontsize=16)
    plt.subplots_adjust(hspace=0)
    cbaxes = fig.add_axes([0.65, 0.49, 0.25, 0.01])
    cb = plt.colorbar(cax = cbaxes,orientation="horizontal")
    cb.update_ticks()

    plt.savefig('./{0}/{0}_comparison.png'.format(pulsar))
