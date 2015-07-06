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
import sys
import os.path
pb.ioff()

# Read command line arguments
parser = argparse.ArgumentParser(description='Pulsar profile variability studies using GPs')
parser.add_argument('-f','--filename', help='File 2d data set', required=True)
parser.add_argument('-p','--pulsar', help='Pulsar name', required=True)
parser.add_argument('-i','--interval', help='inference interval', type=int, required=True)

args = parser.parse_args()

print 'Read arguments'
filename = args.filename
filebase = basename(filename)
outfile = os.path.splitext(filebase)[0]
datfile = outfile + '.dat'
inferreddatafile = outfile + '_inferred_array.dat'
linferreddatafile = outfile + '_linferred_array.dat'
inferredvarfile = outfile + '_inferred_var.dat'
pulsar = args.pulsar
interval = args.interval
data = np.loadtxt(filename)

data = data/100
if os.path.exists('./{0}/mjd.txt'.format(pulsar)):
    mjd = np.loadtxt('./{0}/mjd.txt'.format(pulsar))
if os.path.exists('./{0}/mjd_norm.txt'.format(pulsar)):
    mjd = np.loadtxt('./{0}/mjd_norm.txt'.format(pulsar))
#mjdremoved = np.loadtxt('./{0}/mjdremoved.txt'.format(pulsar))
readbins = np.loadtxt('./{0}/{1}'.format(pulsar,datfile))
leftbin = readbins[0]
rightbin = readbins[1]
allbins = readbins[2]

bins = data.shape[0]
profiles = data.shape[1]

print "Vgp: Pulsar", pulsar
print "-----"
print "-----"

print bins,profiles, np.std(data[0:20,:]), np.std(data[bins-21:-1,:])
print 'size',data.shape
# find properties of noise, assuming data have been prepared using Valign
noiserms = np.mean((np.std(data[0:20,:]), np.std(data[bins-21:-1,:])))
print 'RMS of the noise is:', noiserms 

# Normalize by rms, so everything is now in units of snr
#data = data/noiserms

# build template model to subtract
template = np.median(data,1)

# create difference data
difference = np.zeros((bins,profiles))
for i in range(profiles):
    difference[:,i] = data[:,i] - template
    
maxdifference = np.amax(difference)
mindifference = np.amin(difference)
#limitdifference = np.max((maxdifference, np.abs(mindifference)))

# mjds for inference
#print 'MJD', mjd[0], mjd[-1]
mjdinfer = np.arange(mjd[0],mjd[-1],interval)
profilesinfer = mjdinfer.shape[0]
Llim = np.zeros(profilesinfer)
Ulim = np.zeros(profilesinfer)
xtraining = mjd
inferredarray = np.zeros((bins, profilesinfer))
inferredvar = np.zeros((bins, profilesinfer))

# Run through each bin, train a GP and infer at chosen interval
#if filename[-8:-4] == 'norm':

if not (os.path.exists('./{0}/{1}_bins/'.format(pulsar,outfile))):
    os.mkdir('./{0}/{1}_bins/'.format(pulsar,outfile))  

for i in range(bins):
    ytraining=difference[i,:]
    inferredarray[i,:], inferredvar[i,:], model_params = Vf.gpinferred(xtraining, ytraining, mjdinfer, noiserms)
    Ulim[:] = inferredarray[i,:] + 2 * np.sqrt(inferredvar[i,:])
    Llim[:] = inferredarray[i,:] - 2 * np.sqrt(inferredvar[i,:])
    print "********** GP operating on bin",i+1,"of",bins,"for pulsar",pulsar,"**********"
    plt.plot(xtraining, difference[i,:],'r.')
    plt.plot(mjdinfer, inferredarray[i,:], 'b-')
    plt.fill_between(mjdinfer, Llim, Ulim, color = 'b', alpha = 0.2)
    x1,x2,y1,y2 = plt.axis()
#    plt.text(55000,0.005,'{0}' .format(model_params),fontsize=8)
    plt.axis((x1,x2,mindifference,maxdifference))
    plt.savefig('./{0}/{1}_bins/bin{2}.png'.format(pulsar,outfile,int(i+leftbin)))
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

outputfile = './{0}/{1}' .format(pulsar, linferreddatafile)
np.savetxt(outputfile, linferredarray)
outputfile = './{0}/{1}' .format(pulsar, inferreddatafile)
np.savetxt(outputfile, inferredarray)
outputfile = './{0}/{1}' .format(pulsar, inferredvarfile)
np.savetxt(outputfile, inferredvar)

f = open('./{0}/{0}_outfile.dat' .format(pulsar), 'w')
f.write(outfile)
f.close()

firstmjd = [mjd[0]]
np.savetxt('./{0}/first_emission_mjd.txt' .format(pulsar), firstmjd)
