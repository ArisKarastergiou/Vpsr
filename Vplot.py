#!/usr/bin/python
# Produces the combined emission and spin-down plots, and plots of the individual observations

import argparse
import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt
import Vfunctions as Vf

# Read command line arguments
parser = argparse.ArgumentParser(description='Produces plots for variability studies')
parser.add_argument('-p','--pulsar', help='Pulsar name', required=True)
parser.add_argument('-c','--combined', help='Produces combined emission and spindown plot', action='store_true')
# parser.add_argument('-d','--diagnosticplots', help='make image plots', action='store_true')
args = parser.parse_args()

print 'Read arguments'
pulsar = args.pulsar
interval = 1
datfile = 'zoomed_J1602-5100_bins_203-320.dat'

linferredarray = np.loadtxt('{0}/{0}_linferred_array.dat' .format(pulsar))
inferredarray = np.loadtxt('{0}/{0}_inferred_array.dat' .format(pulsar))
inferredvar = np.loadtxt('{0}/{0}_inferred_var.dat' .format(pulsar))
f = open('{0}/{0}_outfile.dat' .format(pulsar), "r")
outfile = f.read()
f.close()
# outfile = np.loadtxt('{0}/{0}_outfile.dat' .format(pulsar))

mjd = np.loadtxt('./{0}/mjd.txt'.format(pulsar))
mjdremoved = np.loadtxt('./{0}/mjdremoved.txt'.format(pulsar))
readbins = np.loadtxt('./{0}/{1}'.format(pulsar,datfile))
mjdinfer = np.arange(mjd[0],mjd[-1],interval)

leftbin = readbins[0]
rightbin = readbins[1]
allbins = readbins[2]
bins = rightbin-leftbin
yaxis = np.linspace(leftbin/allbins, rightbin/allbins, bins)
maxdifference = np.amax(inferredarray)
mindifference = np.amin(inferredarray)
limitdifference = np.max((maxdifference, np.abs(mindifference)))


yaxis = np.linspace(leftbin/allbins, rightbin/allbins, bins)
Vf.makemap(inferredarray, -limitdifference, limitdifference, mjdinfer, yaxis, mjd, mjdremoved, 'MJD', 'Pulse Phase', pulsar, './{0}/{1}_inferreddata.png'.format(pulsar,outfile))
Vf.makemap(linferredarray, -np.log10(limitdifference), np.log10(limitdifference), mjdinfer, yaxis, mjd, mjdremoved, 'MJD', 'Pulse Phase', pulsar, './{0}/{1}_linferreddata.png'.format(pulsar,outfile))
Vf.makemap(inferredvar, 0 , np.amax(inferredvar), mjdinfer, yaxis, mjd, mjdremoved, 'MJD', 'Pulse Phase', pulsar, './{0}/{1}_inferredvariance.png'.format(pulsar,outfile))



if (args.combined):

    nudot = np.loadtxt('{0}/{0}_nudot.dat' .format(pulsar))
    Llim2 = np.loadtxt('{0}/{0}_Llim2.dat' .format(pulsar))
    Ulim2 = np.loadtxt('{0}/{0}_Ulim2.dat' .format(pulsar))
    mjdinfer_spindown = np.loadtxt('./{0}/{0}_mjdinfer_spindown.dat'.format(pulsar))
    power = int((-1)*np.floor(math.log10((-1)*np.mean(nudot)))) 
    nudot = nudot*10**power
    Llim2 = Llim2*10**power
    Ulim2 = Ulim2*10**power

    fig=plt.figure()
    fig.set_size_inches(16,10)
    ax=fig.add_subplot(2,1,1)
    ax.xaxis.set_visible(False)

    Vf.makemap(linferredarray, -np.log10(limitdifference), np.log10(loglimitdifference), mjdinfer, yaxis, mjd, mjdremoved, 'MJD', 'Pulse phase', pulsar, '.not_saved.png', combined=True)

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
