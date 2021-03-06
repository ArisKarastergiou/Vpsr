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
import pylab
import scipy as sc
import scipy.signal as ss
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
parser.add_argument('-pr','--probplot', help='make probability plots', action='store_true')
parser.add_argument('-r1','--rbf1', nargs=2 ,default = (30,1000), help='lengthscale boundaries 1', type = float, required = False)
parser.add_argument('-r2','--rbf2', nargs=2 ,default = (30,1000),help='lengthscale boundaries 2', type = float, required = False)
parser.add_argument('-toae','--toae' ,default = -1.0, help='toa error', type = float, required = False)
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
TOAerror = args.toae
residtmp = np.loadtxt(filename)

if not (os.path.exists('./{0}/'.format(pulsar))):
    os.mkdir('./{0}/'.format(pulsar))

# Read residuals with one day minimum distance
comparison = residtmp[0,0]
rowstodelete = []
for i in range(1, residtmp.shape[0]):
    if residtmp[i,0] < comparison + 10.0:
        rowstodelete.append(i)
    else:
        comparison = residtmp[i,0]

residuals = np.delete(residtmp, rowstodelete, axis=0)

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
if TOAerror == -1.0:
    TOAerror = np.median(residuals[:,2])
mjdinfer = np.arange(int(mjdfirst), int(mjdlast), 1)

# Train a Gaussian process on the residuals
kernel1 = GPy.kern.RBF(1)
kernel2 = GPy.kern.RBF(1)
kernel = kernel1+kernel2

xtraining1 = xtraining.reshape(xtraining.shape[0],1)
ytraining1 = ytraining.reshape(ytraining.shape[0],1)
xnew = mjdinfer.reshape(mjdinfer.shape[0],1)
model = GPy.models.GPRegression(xtraining1,ytraining1,kernel)

model.add.rbf_1.lengthscale.constrain_bounded(rbf1s, rbf1e)
model.add.rbf_2.lengthscale.constrain_bounded(rbf2s, rbf2e)
#model.Gaussian_noise.variance.constrain_bounded(0,5*TOAerror*TOAerror)
model.optimize()
model.optimize_restarts(num_restarts = 10)
print "MODEL FOR PULSAR",pulsar, model
ypredict, yvariance = model.predict(xnew)
ymodel, yvarmodel = model.predict(xtraining1)

#f = open('model_params_2.txt','a')
#f.write('{0}_2kern {1:.3e} {2:.2f} {3:.3e} {4:.2f} {5:.3e}\n' .format(pulsar,model[0],model[1],model[2],model[3],model[4]))
#f.close()

Ulim = ypredict + 2*np.sqrt(np.abs(yvariance))
Llim = ypredict - 2*np.sqrt(np.abs(yvariance))

# Now use the Gaussian Process to obtain the derivatives, i.e. nudot
parA = np.zeros(2)
parB = np.zeros(2)
K1 = kernel.K(xtraining1, xtraining1)
K1invOut = GPy.util.linalg.pdinv(np.matrix(K1))
K1inv = K1invOut[1]
X_TRAINING = np.matrix([xtraining]).T
#XPREDICT = np.matrix([mjdinfer]).T
XPREDICT = X_TRAINING
Y_TRAINING = np.matrix(np.array(ytraining.flatten())).T

# First lengthscale kernel
CovFunc = Vf.DKD
parA[0] = model['add.rbf_1.variance'] # use the optimized amplitude
parA[1] = model['add.rbf_1.lengthscale'] # use the optimized lengthscale
K_prime = CovFunc(XPREDICT, X_TRAINING, parA) # training points
K_prime_p = 3*parA[0]/parA[1]**4 # These are the diagonal elements of the variance
# Second lengthscale kernel
parB[0] = model['add.rbf_2.variance'] # use the optimized amplitude
parB[1] = model['add.rbf_2.lengthscale'] # use the optimized lengthscale
K_prime += CovFunc(XPREDICT, X_TRAINING, parB) # training points
K_prime_p += 3*parB[0]/parB[1]**4 # These are the diagonal elements of the variance


CovFunc2 = Vf.TKD # for nu double dot
# Just use the longest timescale for nu2dot
#if parA[1]<parB[1]:
#    parA[1] = parB[1]
#    parA[0] = parB[0]
K_prime2 = CovFunc2(XPREDICT, X_TRAINING, parA) # training points
K_prime_p2 = 15*parA[0]/parA[1]**6 # These are the diagonal elements of the variance
K_prime2 += CovFunc2(XPREDICT, X_TRAINING, parA) # training points
K_prime_p2 += 15*parB[0]/parB[1]**6 # These are the diagonal elements of the variance

# Now work out nudot and errors
KiKx, _ = GPy.util.linalg.dpotrs(K1inv, np.asfortranarray(K_prime.T), lower = 1)
KiKx2, _ = GPy.util.linalg.dpotrs(K1inv, np.asfortranarray(K_prime2.T), lower = 1)
#-------
#mu = np.dot(KiKx.T, self.likelihood.Y)
nudot = np.array(nudot0  + np.dot(KiKx.T, Y_TRAINING)/period/(86400)**2)
nu2dot = np.array(np.dot(KiKx2.T, Y_TRAINING)/period/(86400)**3)
#-------

#Kxx = self.kern.Kdiag(_Xnew, which_parts=which_parts)
#var = Kxx - np.sum(np.multiply(KiKx, Kx), 0)
#var = var[:, None]
nudot_err = np.array(np.sqrt(np.abs(K_prime_p - np.sum(np.multiply(KiKx, K_prime.T),0).T))/(86400.)**2)
nu2dot_err = np.array(np.sqrt(np.abs(K_prime_p2 - np.sum(np.multiply(KiKx2, K_prime2.T),0).T))/(86400.)**3)
print "Average nudot error is:", np.mean(nudot_err)
print "The model predicts sigma TOA of:", np.sqrt(model['Gaussian_noise.variance'])
print "The data claim a median sigma TOA of:", np.median(residuals[:,2])

br_index = nu2dot/nudot**2/period

# Limits of nudot plot
#Ulim2 = np.array(nudot + 2*nudot_err)
#Llim2 = np.array(nudot - 2*nudot_err)
errorbars = np.array(2*nudot_err)

# Limits of nu2dot plot
Ulim2n = np.array(nu2dot + 2*nu2dot_err)
Llim2n = np.array(nu2dot - 2*nu2dot_err)
#errorbars = np.array(5*nu2dot_err)

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

resid_resid = (ymodel -ytraining1)*1000
resid_resid_err = residuals[:,2]*1000

if (args.probplot):
    start=50
    end=1000
    step = 50
    dims = int((end-start)/step)+1
    loggrid = np.zeros((dims,dims))
    iaxis = 0
    jaxis = 0
    xtraining1 = xtraining.reshape(xtraining.shape[0],1)
    ytraining1 = ytraining.reshape(ytraining.shape[0],1)
    for i in np.arange(start,end,step):
        for j in np.arange(i,end,step):
            k1 = GPy.kern.RBF(1)
            k2 = GPy.kern.RBF(1)
            k3 = k1 + k2
            m3 = GPy.models.GPRegression(xtraining1,ytraining1,k3)
            m3.add.rbf_1.lengthscale.constrain_fixed(float(i))
            m3.add.rbf_2.lengthscale.constrain_fixed(float(j))
            m3.optimize()
            print iaxis,jaxis,m3
#            model.optimize_restarts(num_restarts = 5)
            loggrid[iaxis,jaxis] = float(m3.log_likelihood())
            del m3, k1, k2, k3
            jaxis += 1
        iaxis += 1
        jaxis = iaxis
        
    print "maximum at:", np.unravel_index(loggrid.argmax(), loggrid.shape), np.max(loggrid)
    plt.pcolormesh(loggrid,vmin=np.max(loggrid) - 100.)
#    plt.axis([start, end, start, end])
    plt.set_cmap('Blues')
    plt.colorbar()
    plt.savefig('./{0}/prob.png'.format(pulsar))
    plt.clf()


# Make plots
if (args.diagnosticplots):
    fig=plt.figure()
    fig.set_size_inches(16,10)
    ax=fig.add_subplot(2,1,1)
    plt.tick_params(axis='both', which='major', labelsize=20)
    plt.plot(xtraining, ytraining,'r.')
    plt.plot(mjdinfer, ypredict, 'k-')
    plt.fill_between(xnew[:,0], Llim[:,0], Ulim[:,0], color = 'k', alpha = 0.2)
    plt.xlabel('Modified Julian Date', fontsize=16)
    ax.xaxis.set_visible(False)
    ax.grid()
    plt.ylim(np.min(ypredict),np.max(ypredict))
    ax=fig.add_subplot(2,1,2)
    plt.tick_params(axis='both', which='major', labelsize=20)
    plt.plot(xtraining, resid_resid,'k.')
    plt.errorbar(xtraining, resid_resid, yerr=resid_resid_err, fmt='.',color = 'k')
    ax.grid()
    plt.xlabel('Modified Julian Date', fontsize=20)
    plt.ylim(np.min((resid_resid)-(2*resid_resid_err)),np.max((resid_resid)+(2*resid_resid_err)))

# makes histogram of data - model values

    # a = plt.axes([.65, .4, .2, .2])
    # n, bins, patches = plt.hist(resid_resid,50, color='k')
    # plt.xlabel("Data - Model (mS)")
    # plt.ylabel("Frequency")
    # fig.text(0.7, 0.2, 'GP Noise Variance: {0:.3e}' .format(model[4]), ha='center', va='center', size=16)
    fig.text(0.06, 0.72, 'Timing Residuals (s)', ha='center', va='center',rotation='vertical', size=20)
    fig.text(0.06, 0.3, 'Model - Data (ms)', ha='center', va='center', rotation='vertical', size=20)
    ax.xaxis.labelpad = 20

    plt.subplots_adjust(hspace=0.1)
    plt.savefig('./{0}/residuals_2kern.png'.format(pulsar))
    plt.close()
    x=np.squeeze(xtraining1)
#    x=np.squeeze(mjdinfer)
    y=np.squeeze(nudot)
    ye=np.squeeze(2*nudot_err)
    plt.errorbar(x[3:-3], y[3:-3], yerr=ye[3:-3],linestyle='-')
#    plt.errorbar(x[1:-1], y[1:-1], yerr=ye[1:-1],linestyle='-')
    plt.xlabel('MJD of Simulation')
    plt.ylabel(r'$\mathrm{{\dot{{\nu}}}}$ ($\mathrm{{s^{{-2}}}}$)')
    plt.subplots_adjust(hspace=0)
    plt.savefig('./{0}/nudot_2kern.png'.format(pulsar))
    plt.clf()

    y=np.squeeze(br_index)
    ye=np.squeeze(2*nu2dot_err)
    plt.errorbar(x[3:-3], y[3:-3], yerr=ye[3:-3],linestyle='-')
    plt.xlabel('MJD of Simulation')
    plt.ylabel(r'$\mathrm{{\dot{{\nu}}}}$ ($\mathrm{{s^{{-2}}}}$)')
    plt.subplots_adjust(hspace=0)
    plt.savefig('./{0}/brindex.png'.format(pulsar))
    plt.clf()

    y=np.squeeze(nu2dot)
    ye=np.squeeze(2*nu2dot_err)
    plt.errorbar(x[3:-3], y[3:-3], yerr=ye[3:-3],linestyle='-')
    plt.xlabel('MJD of Simulation')
    plt.ylabel(r'$\mathrm{{\dot{{\nu}}}}$ ($\mathrm{{s^{{-2}}}}$)')
    plt.subplots_adjust(hspace=0)
    plt.savefig('./{0}/nu2dot.png'.format(pulsar))
    plt.clf()

    # k1 = GPy.kern.RBF(1)
    # model1 = GPy.models.GPRegression(xtraining1,ytraining1,k1)
    # rbf1D = np.zeros(50)
    # for i in range(0,50):
    #     ii = 30 + 20*i
    #     model1.rbf.lengthscale.constrain_fixed(ii)
    #     model1.optimize()
    #     rbf1D[i]=model1.log_likelihood()
    firstmjd = [mjdfirst]
    np.savetxt('./{0}/first_nudot_mjd.txt' .format(pulsar), firstmjd)
