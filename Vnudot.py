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
import sys
import Vfunctions as Vf
import os 
from os.path import basename

# Read command line arguments
#------------------------------
parser = argparse.ArgumentParser(description='Pulsar nudot variability studies using GPs')
parser.add_argument('-f','--filename', help='File containing residuals', required=True)
parser.add_argument('-e','--parfile', help='ephemeris for nudot', required=True)
parser.add_argument('-p','--pulsar', help='Pulsar name', required=True)
parser.add_argument('-d','--diagnosticplots', help='make image plots', action='store_true')
parser.add_argument('-r1','--rbf1', nargs=2 ,default = (30,90), help='lengthscale boundaries 1', type = float, required = False)
parser.add_argument('-r2','--rbf2', nargs=2 ,default = (100, 1000),help='lengthscale boundaries 2', type = float, required = False)
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
rbf2s = args.rbf2[0]
rbf2e = args.rbf2[1]

residuals = np.loadtxt(filename)
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

mjd = residuals[:,0] + epoch
mjdfirst = mjd[0]
mjdlast = mjd[-1]
xtraining = mjd
ytraining = residuals[:,1]
meanerror = np.std(residuals[:,1])
TOAerror = np.median(residuals[:,2])
mjdinfer = np.arange(int(mjdfirst), int(mjdlast), 1)

# Train a Gaussian process on the residuals

kernel1 = GPy.kern.rbf(1)
kernel2 = GPy.kern.rbf(1)
kernel = kernel1+kernel2


xtraining1 = xtraining.reshape(xtraining.shape[0],1)
ytraining1 = ytraining.reshape(ytraining.shape[0],1)
xnew = mjdinfer.reshape(mjdinfer.shape[0],1)
model = GPy.models.GPRegression(xtraining1,ytraining1,kernel, normalize_X=False)
model.constrain_bounded('rbf_1_lengthscale', rbf1s, rbf1e)
model.constrain_bounded('rbf_2_lengthscale', rbf2s, rbf2e)
#model.constrain_bounded('noise_variance',0, TOAerror)
model.optimize()
model.optimize_restarts(num_restarts = 10)
print "MODEL FOR PULSAR",pulsar, model
ypredict, yvariance, a, b = model.predict(xnew)

Ulim = ypredict + 2*np.sqrt(yvariance)
Llim = ypredict - 2*np.sqrt(yvariance)

# Now use the Gaussian Process to obtain the derivatives, i.e. nudot
par = np.zeros(2)
K1 = kernel.K(xtraining1, xtraining1)
K1invOut = GPy.util.linalg.pdinv(np.matrix(K1))
K1inv = K1invOut[1]
X_TRAINING = np.matrix([xtraining]).T
XPREDICT = np.matrix([mjdinfer]).T
Y_TRAINING = np.matrix(np.array(ytraining.flatten())).T

# First lengthscale kernel
CovFunc = Vf.DKD
par[0] = model['rbf_1_variance'] # use the optimized amplitude
par[1] = model['rbf_1_lengthscale'] # use the optimized lengthscale
K_prime = CovFunc(XPREDICT, X_TRAINING, par) # training points
K_prime_p = 3*par[0]/par[1]**4 # These are the diagonal elements of the variance
# Second lengthscale kernel
par[0] = model['rbf_2_variance'] # use the optimized amplitude
par[1] = model['rbf_2_lengthscale'] # use the optimized lengthscale
K_prime += CovFunc(XPREDICT, X_TRAINING, par) # training points
K_prime_p += 3*par[0]/par[1]**4 # These are the diagonal elements of the variance


# Now work out nudot and errors
KiKx, _ = GPy.util.linalg.dpotrs(K1inv, np.asfortranarray(K_prime.T), lower = 1)
#-------
#mu = np.dot(KiKx.T, self.likelihood.Y)
nudot = nudot0  + np.dot(KiKx.T, Y_TRAINING)/period/86400**2
#-------

#Kxx = self.kern.Kdiag(_Xnew, which_parts=which_parts)
#var = Kxx - np.sum(np.multiply(KiKx, Kx), 0)
#var = var[:, None]
nudot_err = np.sqrt(K_prime_p - np.sum(np.multiply(KiKx, K_prime.T),0).T)/86400**2
print "Average nudot error is:", np.mean(nudot_err)

# Limits of nudot plot
Ulim2 = np.array(nudot + 2*nudot_err)
Llim2 = np.array(nudot - 2*nudot_err)

# Write outputs
outputfile = '{0}/{0}_nudot.dat' .format(pulsar)
np.savetxt(outputfile, nudot)
outputfile = '{0}/{0}_Llim2.dat' .format(pulsar)
np.savetxt(outputfile, Llim2)
outputfile = '{0}/{0}_Ulim2.dat' .format(pulsar)
np.savetxt(outputfile, Ulim2)
outputfile = '{0}/{0}_mjdinfer_spindown.dat' .format(pulsar)
np.savetxt(outputfile, mjdinfer)

resid_resid = []
for i in range(xtraining.shape[0]):
    idx = np.argmin(np.abs(mjdinfer - xtraining[i]))
    resid_resid.append(ytraining[i]-ypredict[idx])

# Make plots
if (args.diagnosticplots):
    fig=plt.figure()
    fig.set_size_inches(16,10)
    ax=fig.add_subplot(2,1,1)
    plt.plot(xtraining, ytraining,'r.')
    plt.plot(mjdinfer, ypredict, 'b-')
    plt.fill_between(xnew[:,0], Llim[:,0], Ulim[:,0], color = 'b', alpha = 0.2)
    plt.xlabel('Modified Julian Days')
    plt.ylabel('Timing Residuals (Sec)')
    ax.xaxis.set_visible(False)
    ax=fig.add_subplot(2,1,2)
    plt.plot(xtraining, resid_resid,'k.')
    ax.grid()
    plt.xlabel('Modified Julian Days')
    plt.ylabel('Data - Model (Sec)')
    #x1,x2,y1,y2 = plt.axis()
    #plt.axis((x1,x2,np.min(Llim), np.max(Ulim)))
    plt.savefig('./{0}/residuals.png'.format(pulsar))
    plt.clf()
    plt.plot(mjdinfer, nudot)
    plt.fill_between(xnew[:,0], Llim2[:,0], Ulim2[:,0], color = 'b', alpha = 0.2)
    x1,x2,y1,y2 = plt.axis()
    plt.axis((x1,x2,np.min(Llim2[10:-10]), np.max(Ulim2[10:-10])))
    plt.subplots_adjust(hspace=0)
    plt.savefig('./{0}/nudot.png'.format(pulsar))
    plt.clf()
