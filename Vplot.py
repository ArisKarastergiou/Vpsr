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
parser.add_argument('-f','--final', help='Produces final combined normalised and non-normalised plot', action='store_true')
args = parser.parse_args()

print 'Read arguments'
pulsar = args.pulsar
print "PULSAR IS",pulsar
interval = 1
outfiles_ext = np.sort(np.array(glob.glob('./{0}/zoomed*.txt' .format(pulsar))))
outfiles = []
datfiles = []
yaxis = []
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
    if i ==0 or i==1: # if looking at norm and not-normalised plots, i=0 and i=1 are MP and i=2 and i=3 are IP.
        yaxis.append(np.linspace(0, bins/allbins, bins))

    else:
        yaxis.append(np.linspace(leftbin/allbins, rightbin/allbins, bins))
    maxdifference = np.amax(inferredarray)
    mindifference = np.amin(inferredarray)
    limitdifference = np.max((maxdifference, np.abs(mindifference)))

    if i == 0 or i == 1:
        Vf.makemap(inferredarray, -limitdifference, limitdifference, mjdinfer, yaxis[i], mjd, 'MJD', 'Pulse Phase', pulsar, '{0}_inferreddata.png'.format(outfiles[i]), peakline=allbins/4-leftbin)
        Vf.makemap(linferredarray, -np.log10(limitdifference), np.log10(limitdifference), mjdinfer, yaxis[i], mjd, 'MJD', 'Pulse Phase', pulsar, '{0}_linferreddata.png'.format(outfiles[i]), peakline=allbins/4-leftbin)
        Vf.makemap(inferredvar, 0 , np.amax(inferredvar), mjdinfer, yaxis[i], mjd, 'MJD', 'Pulse Phase', pulsar, '{0}_inferredvariance.png'.format(outfiles[i]), peakline=allbins/4-leftbin)

    else:
        Vf.makemap(inferredarray, -limitdifference, limitdifference, mjdinfer, yaxis[i], mjd, 'MJD', 'Pulse Phase', pulsar, '{0}_inferreddata.png'.format(outfiles[i]))
        Vf.makemap(linferredarray, -np.log10(limitdifference), np.log10(limitdifference), mjdinfer, yaxis[i], mjd, 'MJD', 'Pulse Phase', pulsar, '{0}_linferreddata.png'.format(outfiles[i]))
        Vf.makemap(inferredvar, 0 , np.amax(inferredvar), mjdinfer, yaxis[i], mjd, 'MJD', 'Pulse Phase', pulsar, '{0}_inferredvariance.png'.format(outfiles[i]))
        
#     if (args.combined):

#         nudot = np.loadtxt('{0}/{0}_nudot.dat' .format(pulsar))
#         Llim2 = np.loadtxt('{0}/{0}_Llim2.dat' .format(pulsar))
#         Ulim2 = np.loadtxt('{0}/{0}_Ulim2.dat' .format(pulsar))
#         mjdinfer_spindown = np.loadtxt('./{0}/{0}_mjdinfer_spindown.dat'.format(pulsar))
#         power = int((-1)*np.floor(math.log10(abs(np.median(nudot)))))
#         nudot = nudot*10**power
#         Llim2 = Llim2*10**power
#         Ulim2 = Ulim2*10**power

#         fig=plt.figure()
#         fig.set_size_inches(16,10)
#         ax=fig.add_subplot(2,1,1)
#         ax.xaxis.set_visible(False)

#         if i == 0 or i == 1:
#             Vf.makemap(inferredarray, -limitdifference, limitdifference, mjdinfer, yaxis[i], mjd, 'MJD', 'Pulse phase', pulsar, '.not_saved.png',peakline=allbins/4-leftbin, combined=True)
#         else:
#             Vf.makemap(inferredarray, -limitdifference, limitdifference, mjdinfer, yaxis[i], mjd, 'MJD', 'Pulse phase', pulsar, '.not_saved.png',combined=True)


#         ax=fig.add_subplot(2,1,2)
        
#         ax.grid()

#         plt.plot(mjdinfer_spindown, nudot)
#         plt.fill_between(mjdinfer_spindown, Llim2, Ulim2, color = 'b', alpha = 0.2)
#         plt.xlim(mjdinfer[0],mjdinfer[-1])
#         start_mjd_diff = int(abs(mjdinfer[0]-mjdinfer_spindown[0]))
#         end_mjd_diff = int(abs(mjdinfer[-1]-mjdinfer_spindown[-1]))
#         if end_mjd_diff == 0:
#             end_mjd_diff = 1
# #        plt.ylim(np.median(nudot)-2*np.std(nudot),np.median(nudot)+2*np.std(nudot))
#         plt.ylim(np.min(nudot[start_mjd_diff:-end_mjd_diff]),np.max(nudot[start_mjd_diff:-end_mjd_diff]))

#         plt.xlabel('MJD')
#         plt.ylabel(r'$\mathrm{{\dot{{\nu}}}}$ ($\mathrm{{10^{{-{0}}} s^{{-2}}}}$)'.format(power) ,fontsize=16)
#         y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
#         ax.yaxis.set_major_formatter(y_formatter)
#         plt.subplots_adjust(hspace=0)
#         cbaxes = fig.add_axes([0.65, 0.49, 0.25, 0.01])
#         cb = plt.colorbar(cax = cbaxes,orientation="horizontal")
#         cb.update_ticks()

#         plt.savefig('./{0}/{0}_comparison_{1}.png'.format(pulsar,i))

        
    if i == 0:
        leftbin_initial = leftbin # to use later when plotting peakline


if (args.final):

    nudot = np.loadtxt('{0}/{0}_nudot.dat' .format(pulsar))
#   Llim2 = np.loadtxt('{0}/{0}_Llim2.dat' .format(pulsar))
#   Ulim2 = np.loadtxt('{0}/{0}_Ulim2.dat' .format(pulsar))
    errorbars = np.loadtxt('{0}/{0}_errorbars.dat' .format(pulsar))
    mjdinfer_spindown = np.loadtxt('./{0}/{0}_mjdinfer_spindown.dat'.format(pulsar))
    mjd_norm = np.loadtxt('./{0}/mjd_norm.txt'.format(pulsar))
        
    for i in range(outfiles_ext.shape[0]/2):

        inferredarray = np.loadtxt('{0}_inferred_array.dat' .format(outfiles[2*i]))
        inferredarray_norm = np.loadtxt('{0}_inferred_array.dat' .format(outfiles[2*i+1]))
        maxdifference_norm = np.amax(inferredarray_norm)
        mindifference_norm = np.amin(inferredarray_norm)
        limitdifference_norm = np.max((maxdifference_norm, np.abs(mindifference_norm)))

        maxdifference = np.amax(inferredarray)
        mindifference = np.amin(inferredarray)
        limitdifference = np.max((maxdifference, np.abs(mindifference)))

        if i == 0:
            Vf.combined_map(i,inferredarray_norm,inferredarray, -limitdifference_norm, limitdifference_norm, -limitdifference, limitdifference, mjdinfer, yaxis[i], mjd_norm,mjd,nudot,errorbars,mjdinfer_spindown, 'MJD', 'Phase fraction', pulsar, peakline=allbins/4-leftbin_initial)
        else:
            Vf.combined_map(i,inferredarray_norm,inferredarray, -limitdifference_norm, limitdifference_norm, -limitdifference, limitdifference, mjdinfer, yaxis[i+1], mjd_norm,mjd,nudot,errorbars,mjdinfer_spindown, 'MJD', 'Phase fraction', pulsar)
        

