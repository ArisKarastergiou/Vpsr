#!/usr/bin/python
# Reads in 2d profile data(epoch, bin), aligns, normalizes, rejects outliers and writes out

import argparse
#import pylab as pb
#pb.ion()
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
#pb.ioff()

# Read command line arguments
parser = argparse.ArgumentParser(description='Pulsar profile variability studies using GPs')
parser.add_argument('-f','--filename', help='File 2d data set', required=True)
parser.add_argument('-p','--pulsar', help='Pulsar name', required=True)
parser.add_argument('-i','--interval', help='inference interval', type=int, required=True)
parser.add_argument('-m','--modelprofiles', help='recontruct the model profiles', action='store_true')
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

if outfile[-5:] == '_norm':
    suffix = '_norm'
else:
    suffix = ''

#data = data/100
if os.path.exists('./{0}/mjd.txt'.format(pulsar)):
    mjd = np.loadtxt('./{0}/mjd.txt'.format(pulsar))
if os.path.exists('./{0}/mjd_norm.txt'.format(pulsar)):
    mjd = np.loadtxt('./{0}/mjd_norm.txt'.format(pulsar))

noiserms = np.loadtxt('./{0}/{0}_off_pulse_std{1}.txt'.format(pulsar,suffix))

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
#noiserms = np.mean((np.std(data[0:20,:]), np.std(data[bins-21:-1,:])))
#print 'RMS of the noise is:', noiserms

# build template model to subtract and find mean rms of off-peak
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
    inferredarray[i,:], inferredvar[i,:], model_params = Vf.gpinferred(xtraining, ytraining, mjdinfer)
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

#Plots the simulation variability maps for the paper

# max_of_inferredarray = np.max(inferredarray)
# min_of_inferredarray = np.min(inferredarray)
# limitdifference = np.max((max_of_inferredarray, np.abs(min_of_inferredarray)))
# my_sim_max = limitdifference
# my_sim_min = -limitdifference

# firstmjd = [mjd[0]]
# np.savetxt('./{0}/first_emission_mjd.txt' .format(pulsar), firstmjd)
# fig = plt.figure()
# plt.imshow(inferredarray,aspect = 'auto',cmap = "RdBu_r", vmin = my_sim_min, vmax = my_sim_max)
# for i in range(profiles):
#     print 'mjd',i,mjd[i]
#     plt.vlines(mjd[i],0,203,linestyles='dotted')
# plt.ylabel(r'Phase Fraction',fontsize=14)
# plt.xlabel(r'Days',fontsize=14)
# ylocs = [0,50,100,150,200]
# ytlabs = [0/1024.0,50/1024.0,100/1024.0,150/1024.0,200/1024.0]
# ytlabs = np.round(ytlabs,2)
# plt.yticks(ylocs,ytlabs)
# #plt.colorbar()
# plt.hlines(104,0,1770)
# cbaxes_corr = fig.add_axes([0.45, 0.93, 0.45, 0.01])
# cb_corr = plt.colorbar(cax = cbaxes_corr,orientation="horizontal")
# cbaxes_corr.tick_params(labelsize=10) 
# plt.savefig('var_map.png' .format(pulsar))
# plt.clf()

# Plots the variability map in isolation

maxdifference = np.amax(inferredarray)
mindifference = np.amin(inferredarray)
limitdifference = np.max((maxdifference, np.abs(mindifference)))

readbins = np.loadtxt('./{0}/{1}'.format(pulsar,datfile))

leftbin = readbins[0]
rightbin = readbins[1]
allbins = readbins[2]
bins = rightbin-leftbin

yaxis=[]
yaxis.append(np.linspace(0, bins/allbins, bins))

print 'yaxis is',yaxis
print 'yaxis zero is',yaxis[0]

Vf.makemap(inferredarray, -limitdifference, limitdifference, mjdinfer, yaxis[0],mjd, 'Modified Julian Date', 'Fraction Of Pulse Period', pulsar, './{0}/{1}_inferreddata.png'.format(pulsar,outfile), peakline=allbins/4-leftbin)

if (args.modelprofiles):

    if not (os.path.exists('./{0}/model_profiles'.format(pulsar))):
        os.mkdir('./{0}/model_profiles'.format(pulsar))

    no_bins = inferredarray.shape[0]
    no_profiles = inferredarray.shape[1]

    model_interval = 10

    model_profiles = np.zeros((no_profiles,no_bins))

    pre_scaling = inferredarray*noiserms
    for i in range(no_profiles):
        for j in range(no_bins):
            model_profiles[i,j] = pre_scaling[j,i] + template[j]

    plot_no = int(np.floor(no_profiles/model_interval))

    for n in range (plot_no):
        print n+1,'of',plot_no
        fig=plt.figure()
        plt.plot(model_profiles[n*model_interval,:])
        plt.xlim(0,no_bins-1)
        plt.ylim(-0.1,1.3)
        plt.xlabel('Pulse Phase Bin')
        plt.ylabel('Normalised Flux Density')
        if n < 9:
            fig.text(0.68, 0.75, '00{0}/{1}' .format(n+1,plot_no), size=16)
            fig.text(0.68, 0.65, 'Freq. {0} Hz' .format(int(mjd[0])+n*model_interval), size=16)
            plt.savefig('./{0}/model_profiles/00{1}.png' .format(pulsar,n+1))
        if n > 8 and n < 99:
            fig.text(0.68, 0.75, '0{0}/{1}' .format(n+1,plot_no), size=16)
            fig.text(0.68, 0.65, 'MJD. {0}' .format(int(mjd[0])+n*model_interval), size=16)
            plt.savefig('./{0}/model_profiles/0{1}.png' .format(pulsar,n+1))
        if n > 98:
            fig.text(0.68, 0.75, '{0}/{1}' .format(n+1,plot_no), size=16)
            fig.text(0.68, 0.65, 'Freq. {0} Hz' .format(int(mjd[0])+n*model_interval), size=16)
            plt.savefig('./{0}/model_profiles/{1}.png' .format(pulsar,n+1))
        plt.close()
        plt.clf()

    os.system('convert -loop 0 ./{0}/model_profiles/*.png ./{0}/animation.gif' .format(pulsar))

    for i in range(100):
        print 'for profile',i,'sum is',np.sum(model_profiles[i,:])

    np.savetxt('{0}/{0}_model_profiles.txt' .format(pulsar),model_profiles)
