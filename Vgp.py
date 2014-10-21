#!/usr/bin/python
# Reads in 2d profile data(epoch, bin), aligns, normalizes, rejects outliers and writes out

import argparse
import pylab as pb
pb.ion()
import numpy as np
import GPy
import math
import matplotlib
import matplotlib.pyplot as plt
import scipy as sc
import sys
import Vfunctions as Vf
import os 
from os.path import basename

# Read command line arguments
parser = argparse.ArgumentParser(description='Pulsar profile variability studies using GPs')
parser.add_argument('-f','--filename', help='File 2d data set', required=True)
parser.add_argument('-p','--pulsar', help='Pulsar name', required=True)
parser.add_argument('-i','--interval', help='inference interval', type=int, required=True)
parser.add_argument('-d','--diagnosticplots', help='make image plots', action='store_true')

args = parser.parse_args()

print 'Read arguments'
filename = args.filename
filebase = basename(filename)
outfile = os.path.splitext(filebase)[0]
datfile = outfile + '.dat'
pulsar = args.pulsar
interval = args.interval
data = np.loadtxt(filename)


data = data/100
mjd = np.loadtxt('./{0}/mjd.txt'.format(pulsar))
readbins = np.loadtxt('./{0}/{1}'.format(pulsar,datfile))
leftbin = readbins[0]
rightbin = readbins[1]
allbins = readbins[2]

bins = data.shape[0]
profiles = data.shape[1]
print bins,profiles, np.std(data[0:20,:]), np.std(data[bins-21:-1,:])
# find properties of noise, assuming data have been prepared using Valign
noiserms = np.mean((np.std(data[0:20,:]), np.std(data[bins-21:-1,:])))
print 'RMS of the noise is:', noiserms 

# Normalize by rms, so everything is now in units of snr
#data = data/noiserms

# build template model to subtract
template = np.mean(data,1)

# create difference data
difference = np.zeros((bins,profiles))
for i in range(profiles):
    difference[:,i] = data[:,i] - template
    
maxdifference = np.amax(difference)
mindifference = np.amin(difference)
limitdifference = np.max((maxdifference, np.abs(mindifference)))

# mjds for inference
print 'MJD', mjd[0], mjd[-1]
mjdinfer = np.arange(mjd[0],mjd[-1],interval)
profilesinfer = mjdinfer.shape[0]
Llim = np.zeros(profilesinfer)
Ulim = np.zeros(profilesinfer)
xtraining = mjd
inferredarray = np.zeros((bins, profilesinfer))
inferredvar = np.zeros((bins, profilesinfer))

# Run through each bin, train a GP and infer at chosen interval
if (args.diagnosticplots):
    if not (os.path.exists('./{0}/bins/'.format(pulsar))):
        os.mkdir('./{0}/bins/'.format(pulsar))  

for i in range(bins):
    ytraining=difference[i,:]
    inferredarray[i,:], inferredvar[i,:] = Vf.gpinferred(xtraining, ytraining, mjdinfer, noiserms)
    Ulim[:] = inferredarray[i,:] + 2 * np.sqrt(inferredvar[i,:])
    Llim[:] = inferredarray[i,:] - 2 * np.sqrt(inferredvar[i,:])

    if (args.diagnosticplots):
        plt.plot(xtraining, difference[i,:],'r.')
        plt.plot(mjdinfer, inferredarray[i,:], 'b-')
        plt.fill_between(mjdinfer, Llim, Ulim, color = 'b', alpha = 0.2)
        x1,x2,y1,y2 = plt.axis()
        plt.axis((x1,x2,mindifference,maxdifference))
        plt.savefig('./{0}/bins/bin{1}.png'.format(pulsar,int(i+leftbin)))
        plt.clf()
        
inferredarray = inferredarray/noiserms
maxdifference = np.amax(inferredarray)
mindifference = np.amin(inferredarray)
limitdifference = np.max((maxdifference, np.abs(mindifference)))
linferredarray = np.zeros((inferredarray.shape[0],inferredarray.shape[1]))
for i in range(inferredarray.shape[0]):
    for j in range(inferredarray.shape[1]):
        if inferredarray[i,j] > -1.0 and inferredarray[i,j] < 1.0:
            linferredarray[i,j] = 0
        if inferredarray[i,j] > 1.0:
            linferredarray[i,j] = np.log10(inferredarray[i,j])
        if inferredarray[i,j] < -1.0:
            linferredarray[i,j] = -np.log10(-inferredarray[i,j])

print np.amax(linferredarray), np.log10(limitdifference)


if (args.diagnosticplots):
    yaxis = np.linspace(leftbin/allbins, rightbin/allbins, bins)
    Vf.makemap(inferredarray, -limitdifference, limitdifference, mjdinfer, yaxis, mjd, 'MJD dates', 'pulse phase', pulsar, './{0}/{1}_inferreddata.png'.format(pulsar,outfile))
    Vf.makemap(linferredarray, -np.log10(limitdifference), np.log10(limitdifference), mjdinfer, yaxis, mjd, 'MJD dates', 'pulse phase', pulsar, './{0}/{1}_linferreddata.png'.format(pulsar,outfile))
    Vf.makemap(inferredvar, 0 , np.amax(inferredvar), mjdinfer, yaxis, mjd, 'MJD dates', 'pulse phase', pulsar, './{0}/{1}_inferredvariance.png'.format(pulsar,outfile))

        
        
