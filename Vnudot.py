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
parser = argparse.ArgumentParser(description='Pulsar nudot variability studies using GPs')
parser.add_argument('-f','--filename', help='File containing residuals', required=True)
parser.add_argument('-e','--parfile', help='ephemeris for nudot', required=True)
parser.add_argument('-p','--pulsar', help='Pulsar name', required=True)
parser.add_argument('-d','--diagnosticplots', help='make image plots', action='store_true')

args = parser.parse_args()

print 'Read arguments'
filename = args.filename
parfile = args.parfile
filebase = basename(filename)
outfile = os.path.splitext(filebase)[0]
datfile = outfile + '.dat'
pulsar = args.pulsar

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
#assert mjdfirst == np.min(mjd)
mjdlast = mjd[-1]
#assert mjdlast == np.max(mjd)


xtraining = mjd
ytraining = residuals[:,1]
meanerror = np.std(residuals[:,1])
TOAerror = np.std(residuals[:,2])
mjdinfer = np.arange(int(mjdfirst), int(mjdlast), 1)

# Train a Gaussian process on the residuals

kernel1 = GPy.kern.rbf(1)
kernel2 = GPy.kern.white(1)

kernel = kernel1 + kernel2
xtraining1 = xtraining.reshape(xtraining.shape[0],1)
ytraining1 = ytraining.reshape(ytraining.shape[0],1)
xnew = mjdinfer.reshape(mjdinfer.shape[0],1)
model = GPy.models.GPRegression(xtraining1,ytraining1,kernel, normalize_X=False)
model.constrain_bounded('rbf_lengthscale', 15,1000)
#model.constrain_fixed('rbf_variance', 0.16)
model.constrain_bounded('noise_variance',0, TOAerror)
model.constrain_bounded('white_variance',0, TOAerror)
model.optimize()
model.optimize_restarts(num_restarts = 50)
ypredict, yvariance, a, b = model.predict(xnew)
Ulim = ypredict + 2*np.sqrt(yvariance)
Llim = ypredict - 2*np.sqrt(yvariance)
print model
K1 = kernel.K(xtraining1, xtraining1)
print np.diag(K1)
K1inv = np.linalg.inv( np.matrix(K1) )
X_TRAINING = np.matrix([xtraining]).T
XPREDICT = np.matrix([mjdinfer]).T
Y_TRAINING = np.matrix(np.array(ytraining.flatten())).T

CovFunc = Vf.DKD
par = np.zeros(2)
par[0] = model['rbf_variance'] # use the optimized amplitude
par[1] = model['rbf_lengthscale'] # use the optimized lengthscale

K_prime = CovFunc(XPREDICT, X_TRAINING, par) # training points
K_prime_p = 3*par[0]**2/par[1]**4 # These are the diagonal elements of the variance

nudot = np.array(nudot0 + (K_prime * K1inv * Y_TRAINING)/period/86400**2)
nudot_var_mat = np.abs(K_prime_p - K_prime * K1inv * K_prime.T)
nudot_err = np.array(np.sqrt(np.diag(nudot_var_mat)))/86400**2
nudot_err = nudot_err.reshape(nudot_err.shape[0],1)
print nudot.shape, nudot_err.shape
print np.diag(K_prime * K1inv * K_prime.T)
Ulim2 = nudot + 2*nudot_err
Llim2 = nudot - 2*nudot_err

outputfile = '{0}/{0}_nudot.dat' .format(pulsar)
np.savetxt(outputfile, nudot)
outputfile = '{0}/{0}_Llim2.dat' .format(pulsar)
np.savetxt(outputfile, Llim2)
outputfile = '{0}/{0}_Ulim2.dat' .format(pulsar)
np.savetxt(outputfile, Ulim2)
outputfile = '{0}/{0}_mjdinfer_spindown.dat' .format(pulsar)
np.savetxt(outputfile, mjdinfer)

if (args.diagnosticplots):
    plt.plot(xtraining, ytraining,'r.')
    plt.plot(mjdinfer, ypredict, 'b-')
    plt.fill_between(xnew[:,0], Llim[:,0], Ulim[:,0], color = 'b', alpha = 0.2)
    x1,x2,y1,y2 = plt.axis()
    plt.axis((x1,x2,np.min(Llim), np.max(Ulim)))
    plt.savefig('./{0}/residuals.png'.format(pulsar))
    plt.clf()
    plt.plot(mjdinfer, nudot)
    plt.fill_between(xnew[:,0], Llim2[:,0], Ulim2[:,0], color = 'b', alpha = 0.2)
    plt.savefig('./{0}/nudot.png'.format(pulsar))
