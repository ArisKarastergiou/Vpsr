#!/usr/bin/python

# Reads in a file of residuals and the nudot used to produce it and
# produces a smooth nudot vs time array and plot


import argparse
import pylab as pb
pb.ion()
import numpy as np
import GPy
import math
import matplotlib
import matplotlib.pyplot as plt
import scipy as sc
import scipy.signal as ss
import sys
import Vfunctions as Vf
import os 
from os.path import basename
from numpy.random import normal

# Read command line arguments
#------------------------------
parser = argparse.ArgumentParser(description='Pulsar nudot variability studies using GPs')
parser.add_argument('-f','--filename', help='File containing residuals', required=True)
parser.add_argument('-e','--parfile', help='ephemeris for nudot', required=True)
parser.add_argument('-p','--pulsar', help='Pulsar name', required=True)
parser.add_argument('-d','--diagnosticplots', help='make image plots', action='store_true')
parser.add_argument('-r1','--rbf1', nargs=2 ,default = (30,1000), help='lengthscale boundaries 1', type = float, required = False)
parser.add_argument('-toae','--toae', default = -1.0, help='toa error', type = float, required = False)
#parser.add_argument('-r2','--rbf2', nargs=2 ,default = (100, 1000),help='lengthscale boundaries 2', type = float, required = False)
#------------------------------

args = parser.parse_args()
filename = args.filename
parfile = args.parfile
filebase = basename(filename)
outfile = os.path.splitext(filebase)[0]
datfile = outfile + '.dat'
pulsar = args.pulsar
rbf1s = args.rbf1[0]
rbf1e = args.rbf1[1]
TOAerror = args.toae
#print rbf1s,rbf1e
#rbf2s = args.rbf2[0]
#rbf2e = args.rbf2[1]

if not (os.path.exists('./{0}/'.format(pulsar))):
    os.mkdir('./{0}/'.format(pulsar))

residtmp = np.loadtxt(filename)
# Read residuals with one day minimum distance
comparison = residtmp[0,0]
rowstodelete = []
for i in range(1, residtmp.shape[0]):
    if residtmp[i,0] < comparison + 10:
        rowstodelete.append(i)
    else:
        comparison = residtmp[i,0]

residuals = np.delete(residtmp, rowstodelete, axis=0)

# epoch and nudot
q = open(parfile)
for line in q:
    if line.startswith('F0'):
        f0_line = line.split()
        period = 1/float(f0_line[1])
    if line.startswith('F1'):
        f1_line = line.split()
        nudot0 = float(f1_line[1])
    if line.startswith('PEPOCH'):
        pepoch_line = line.split()
        epoch = float(pepoch_line[1])
q.close()

#nudot_infer=1
mjd = residuals[:,0] + epoch
mjdfirst = mjd[0]
mjdlast = mjd[-1]
xtraining = mjd

ytraining = residuals[:,1]
meanerror = np.std(residuals[:,1])
print "Median toa error:", np.median(residuals[:,2])
if TOAerror == -1.0 :
    TOAerror = np.median(residuals[:,2])
mjdinfer = np.arange(int(mjdfirst), int(mjdlast), 1)
#mjdinfer30 = np.arange(int(mjdfirst), int(mjdlast), nudot_infer)
# Train a Gaussian process on the residuals

kernel1 = GPy.kern.RBF(1)
kernel = kernel1

#kernel2 = GPy.kern.rbf(1)
#kernel = kernel1+kernel2

xtraining1 = xtraining.reshape(xtraining.shape[0],1)
ytraining1 = ytraining.reshape(ytraining.shape[0],1)
xnew = mjdinfer.reshape(mjdinfer.shape[0],1)
#xnew30 = mjdinfer30.reshape(mjdinfer30.shape[0],1)
model = GPy.models.GPRegression(xtraining1,ytraining1,kernel)
print model,rbf1s,rbf1e
model.rbf.lengthscale.constrain_bounded(rbf1s, rbf1e)
#model.constrain_bounded('rbf_lengthscale', rbf1s, rbf1e)
model.Gaussian_noise.variance.constrain_bounded(0,5*TOAerror*TOAerror)
#model.constrain_bounded('noise_variance', 0,5*TOAerror*TOAerror)
model.optimize()
model.optimize_restarts(num_restarts = 10)
print "MODEL FOR PULSAR",pulsar, model
ypredict, yvariance = model.predict(xnew)
ymodel, yvarmodel = model.predict(xtraining1)

f = open('model_params_1.txt','a')
f.write('{0}_1kern {1:.3e} {2:.2f} {3:.3e}\n' .format(pulsar,model[0],model[1],model[2]))
f.close()



# Compute residuals at training points
#resid_resid = np.array(ymodel - ytraining1)
#resid_resid_err = 2 * np.sqrt(yvarmodel)
#rr = np.squeeze(resid_resid.T)
#print "rr shape", rr.shape
## Compute autocorrelation function of residuals to test if it is white noise
#autocorr = np.correlate(rr,rr, mode='full')
#smoothac = ss.savgol_filter(autocorr[autocorr.shape[0]/2:],13, 4)
#previous = 1000000
#zerocrossing = 0
#i=0
#while zerocrossing < 3:
#    if (smoothac[i]*previous < 0):
#        print smoothac[i],i, zerocrossing
#        zerocrossing += 1
#        peak = i
#    previous = smoothac[i]
#    i+=1
#print 0.8*peak

#plt.plot(autocorr[autocorr.shape[0]/2:])
#plt.xlabel('Lag')
#plt.ylabel('Autocorrelation')
#plt.savefig('{0}_ac1kern.png' .format(pulsar))
#plt.clf()

Ulim = ypredict + 2*np.sqrt(yvariance)
Llim = ypredict - 2*np.sqrt(yvariance)

# Now use the Gaussian Process to obtain the derivatives, i.e. nudot
par = np.zeros(2)
K1 = kernel.K(xtraining1, xtraining1)
K1invOut = GPy.util.linalg.pdinv(np.matrix(K1))
K1inv = K1invOut[1]
X_TRAINING = np.matrix([xtraining]).T
#XPREDICT = np.matrix([mjdinfer]).T
XPREDICT = X_TRAINING
Y_TRAINING = np.matrix(np.array(ytraining.flatten())).T

# First lengthscale kernel
CovFunc = Vf.DKD
par[0] = model['rbf.variance'] # use the optimized amplitude
par[1] = model['rbf.lengthscale'] # use the optimized lengthscale
K_prime = CovFunc(XPREDICT, X_TRAINING, par) # training points
K_prime_p = 3*par[0]/par[1]**4 # These are the diagonal elements of the variance
# Second lengthscale kernel
# par[0] = model['rbf_2_variance'] # use the optimized amplitude
# par[1] = model['rbf_2_lengthscale'] # use the optimized lengthscale
# K_prime += CovFunc(XPREDICT, X_TRAINING, par) # training points
# K_prime_p += 3*par[0]/par[1]**4 # These are the diagonal elements of the variance


# Now work out nudot and errors
KiKx, _ = GPy.util.linalg.dpotrs(K1inv, np.asfortranarray(K_prime.T), lower = 1)
#-------
#mu = np.dot(KiKx.T, self.likelihood.Y)
nudot = np.array(nudot0  + np.dot(KiKx.T, Y_TRAINING)/period/(86400)**2)
#-------

#Kxx = self.kern.Kdiag(_Xnew, which_parts=which_parts)
#var = Kxx - np.sum(np.multiply(KiKx, Kx), 0)
#var = var[:, None]
nudot_err = np.array(np.sqrt(K_prime_p - np.sum(np.multiply(KiKx, K_prime.T),0).T)/(86400)**2)
print "Average nudot error is:", np.mean(nudot_err)

# Limits of nudot plot
#Ulim2 = np.array(nudot + 2*nudot_err)
#Llim2 = np.array(nudot - 2*nudot_err)
errorbars = np.array(2*nudot_err)

# Write outputs
outputfile = '{0}/{0}_nudot.dat' .format(pulsar)
np.savetxt(outputfile, nudot)
#outputfile = '{0}/{0}_Llim2.dat' .format(pulsar)
#np.savetxt(outputfile, Llim2)
#outputfile = '{0}/{0}_Ulim2.dat' .format(pulsar)
#np.savetxt(outputfile, Ulim2)
outputfile = '{0}/{0}_errorbars.dat' .format(pulsar)
np.savetxt(outputfile, errorbars)
outputfile = '{0}/{0}_mjdinfer_spindown.dat' .format(pulsar)
np.savetxt(outputfile, xtraining)

# Produce a plot showing difference between model and data
# resid_resid = []
# for i in range(xtraining.shape[0]):
#     idx = np.argmin(np.abs(mjdinfer - xtraining[i]))
#     resid_resid.append(ytraining[i]-ypredict[idx])

resid_resid = (ymodel -ytraining1)*1000
resid_resid_err = residuals[:,2]*1000
#2 * np.sqrt(yvarmodel)

# Make plots

if (args.diagnosticplots):
    fig=plt.figure()
    fig.set_size_inches(16,10)
    ax=fig.add_subplot(2,1,1)
    plt.plot(xtraining, ytraining,'r.')
    plt.plot(mjdinfer, ypredict, 'k-')
    plt.fill_between(xnew[:,0], Llim[:,0], Ulim[:,0], color = 'k', alpha = 0.2)
    plt.xlabel('Modified Julian Date', fontsize=16)
    #plt.ylabel('Timing Residuals (Sec)', fontsize=16)
    ax.xaxis.set_visible(False)
    ax.grid()
    plt.ylim(np.min(ypredict),np.max(ypredict))
    ax=fig.add_subplot(2,1,2)

    
    plt.plot(xtraining, resid_resid,'k-')
    plt.errorbar(xtraining, resid_resid, yerr=resid_resid_err, fmt='.',color = 'k') 
    ax.grid()
    plt.xlabel('Modified Julian Date', fontsize=16)
    #plt.ylabel('Data - Model (mS)', fontsize=16)
    #x1,x2,y1,y2 = plt.axis()
    #plt.axis((x1,x2,np.min(Llim), np.max(Ulim)))
    plt.ylim(np.min((resid_resid)-(2*resid_resid_err)),np.max((resid_resid)+(2*resid_resid_err)))

# makes histogram of data - model values

    # a = plt.axes([.65, .4, .2, .2])                                                      
    # n, bins, patches = plt.hist(resid_resid,50, color='k')                                          
    # plt.xlabel("Data - Model (mS)")                                                      
    # plt.ylabel("Frequency")

    # fig.text(0.7, 0.2, 'GP Noise Variance: {0:.3e}' .format(model[2]), ha='center', va='center', size=16)

    fig.text(0.06, 0.72, 'Timing Residuals (s)', ha='center', va='center', rotation='vertical', size=16)
    fig.text(0.06, 0.3, 'Data - Model (ms)', ha='center', va='center', rotation='vertical', size=16)
    ax.xaxis.labelpad = 20
    plt.subplots_adjust(hspace=0.1)
    plt.savefig('./{0}/residuals_1kern.png'.format(pulsar))
    plt.clf()
#    plt.plot(mjdinfer30, nudot, 'r+')
#    x1,x2,y1,y2 = plt.axis()
#    plt.axis((x1,x2,np.min(Llim2[10:-10]), np.max(Ulim2[10:-10])))
#    plt.fill_between(xnew30[30:-30,0], Llim2[30:-30,0], Ulim2[30:-30,0], color = 'b', alpha = 0.2)
#    plt.fill_between(xtraining1[3:-3,0], Llim2[3:-3,0], Ulim2[3:-3,0], color = 'b', alpha = 0.2)
#    x=np.squeeze(mjdinfer)
    x=np.squeeze(xtraining1)
    y=np.squeeze(nudot)
    ye=np.squeeze(2*nudot_err)
    plt.errorbar(x[3:-3], y[3:-3], yerr=ye[3:-3],linestyle='-')
    plt.xlabel('MJD of Simulation')			
    plt.ylabel(r'$\mathrm{{\dot{{\nu}}}}$ ($\mathrm{{s^{{-2}}}}$)')
    plt.subplots_adjust(hspace=0)
    plt.savefig('./{0}/nudot_1kern.png'.format(pulsar))
    plt.clf()
# model.constrain_bounded('rbf_1_lengthscale', rbf1s, rbf1e)
# model.constrain_bounded('rbf_2_lengthscale', rbf2s, rbf2e)
# model.constrain_bounded('noise_variance', 0,TOAerror*TOAerror)
# #model.constrain_fixed('rbf_2_variance',0)
# model.optimize()
# model.optimize_restarts(num_restarts = 10)
# print "MODEL FOR PULSAR",pulsar, model
    firstmjd = [mjdfirst]
    np.savetxt('./{0}/first_nudot_mjd.txt' .format(pulsar), firstmjd)
