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
args = parser.parse_args()

print 'Read arguments'
pulsar = args.pulsar
interval = 1
datfile = 'zoomed_J0738-4042_bins_0-378.dat'

linferredarray = np.loadtxt('{0}/{0}_linferred_array.dat' .format(pulsar))
nudot = np.loadtxt('{0}/{0}_nudot.dat' .format(pulsar))
Llim2 = np.loadtxt('{0}/{0}_Llim2.dat' .format(pulsar))
Ulim2 = np.loadtxt('{0}/{0}_Ulim2.dat' .format(pulsar))
mjd = np.loadtxt('./{0}/mjd.txt'.format(pulsar))
mjdremoved = np.loadtxt('./{0}/mjdremoved.txt'.format(pulsar))
readbins = np.loadtxt('./{0}/{1}'.format(pulsar,datfile))
mjdinfer_spindown = np.loadtxt('./{0}/{0}_mjdinfer_spindown.dat'.format(pulsar))
mjdinfer = np.arange(mjd[0],mjd[-1],interval)

leftbin = readbins[0]
rightbin = readbins[1]
allbins = readbins[2]
bins = rightbin-leftbin
yaxis = np.linspace(leftbin/allbins, rightbin/allbins, bins)
power = int((-1)*np.floor(math.log10((-1)*np.mean(nudot)))) 
nudot = nudot*10**power
Llim2 = Llim2*10**power
Ulim2 = Ulim2*10**power
logmaxdifference = np.amax(linferredarray)
logmindifference = np.amin(linferredarray)
loglimitdifference = np.max((logmaxdifference, np.abs(logmindifference)))

fig=plt.figure()
fig.set_size_inches(16,10)
ax=fig.add_subplot(2,1,1)
ax.xaxis.set_visible(False)

Vf.makemap(linferredarray, -loglimitdifference, loglimitdifference, mjdinfer, yaxis, mjd, mjdremoved, 'MJD', 'Pulse phase', pulsar, '.not_saved.png', combined=True)

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
