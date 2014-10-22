# Functions for V code
import numpy as np
import GPy
import matplotlib.pyplot as plt
import scipy.spatial as sp
import math
import os

def makeoriginalplots(pulsar, data, mjd):
    nbins = data.shape[0]
    nprofiles = data.shape[1]
    if not (os.path.exists('./{0}/original_profiles'.format(pulsar))):
        os.mkdir('./{0}/original_profiles'.format(pulsar))
    for i in range(nprofiles):
        plt.plot(data[:,i])
        plt.suptitle('{0}'.format(mjd[i]), fontsize=14, fontweight='bold')
        plt.savefig('./{0}/original_profiles/{1}_{2}.png' .format(pulsar,int(math.floor(mjd[i])),i))
        plt.clf()

def aligndata(baselineremoved, brightest):
    nbins = baselineremoved.shape[0]
    nprofiles = baselineremoved.shape[1]
    template = baselineremoved[:,brightest]
    # rotate template to put peak at 1/4
    peakbin = np.argmax(template)
    fixedlag = int(nbins/4)-peakbin
    aligned = np.zeros((nbins,nprofiles))
    newtemplate = np.roll(template, fixedlag)
    template = newtemplate
    plt.plot(newtemplate)
    plt.savefig('template.png')
    plt.clf()
    for i in range(nprofiles):
        xcorr = np.correlate(template,baselineremoved[:,i],"full")
        lag = np.argmax(xcorr)
        aligned[:,i] = np.roll(baselineremoved[:,i],lag)
    template = np.mean(aligned,1)
    # repeat with better template now and shift peak to 1/4 of the profile
    peakbin = np.argmax(template)
    fixedlag = int(nbins/4)-peakbin
    double = np.zeros(2*nbins)
    for i in range(nprofiles):
        double[0:nbins] = baselineremoved[:,i]
        double[nbins:2*nbins] = baselineremoved[:,i]
#        xcorr = np.correlate(template,baselineremoved[:,i],"full")
        xcorr = np.correlate(template,double,"full")
        lag = np.argmax(xcorr) + fixedlag
        aligned[:,i] = np.roll(baselineremoved[:,i],lag)
    newtemplate = np.roll(template, fixedlag)
    return np.array(aligned), np.array(newtemplate)

def removebaseline(data, outliers):
    # chop profile into 8 parts and check the part with the lowest rms.
    # Remove the mean of that from everything. Remove outliers based on rms.
    nbins = data.shape[0]
    nprofiles = data.shape[1]
    # initialize output array
    baselineremoved = data
    smallestrms = np.zeros(nprofiles)
    smallestmean = np.zeros(nprofiles)
    for i in range(nprofiles):
        rms = np.zeros(8)
        mean = np.zeros(8)
        section = nbins/8
        for j in range(8):
            rms[j] = np.std(data[j*section:(j+1)*section,i])
            mean[j] = np.mean(data[j*section:(j+1)*section,i])
        smallestrms[i] = np.min(rms)
        baseindex = np.argmin(rms)
        baseline = mean[baseindex]
        smallestmean[i] = baseline
        baselineremoved[:,i] = data[:,i] - baseline
    medianrms = np.median(smallestrms)

    outlierindex = []
    for i in range(nprofiles):
        if smallestrms[i] > outliers * medianrms:
            outlierindex.append(i)
        if smallestrms[i] * outliers < medianrms :
            outlierindex.append(i)
        if np.min(baselineremoved[:,i]) < - 5 * smallestrms[i]:
            outlierindex.append(i)

    ou = np.array(outlierindex)
    baselineoutlierremoved = np.delete(baselineremoved, ou, 1)
    rmsremoved = np.delete(smallestrms, ou)
    print 'Removed outliers: ', ou.shape[0]
    return baselineoutlierremoved, rmsremoved, ou


def findbrightestprofile(data,rmsdata):
    snr = np.zeros(rmsdata.shape[0])
    for i in range(data.shape[1]):
        snr[i] = np.max(data[:,i])/rmsdata[i]
    brightestindex = np.argmax(snr)
    return brightestindex

def binstartend(data,rms):
    peak = np.max(data)
    peakbin = np.argmax(data)
    bins = data.shape[0]
    peaksnr = peak/rms
    power = peaksnr
    peaks = 0
    lstart = 0
    lend = 0
    if peaksnr > 10:
        peaks = 1
        thisbin = peakbin - 1
        while power > 0:
            if peakbin != 0 :
                power = data[thisbin]
                data[thisbin] = 0
                thisbin = thisbin - 1
        lstart = thisbin
        thisbin = peakbin
        power = peaksnr
        while power > 0:
            if peakbin != bins - 1 :
                power = data[thisbin]
                data[thisbin] = 0
                thisbin = thisbin + 1
        lend = thisbin
    start = np.max((lstart - 20, 0))
    end = np.min((lend + 20, bins - 1))
    cuttemplate = np.array(data)
    return start, end, peaks, cuttemplate


def gpinferred(xtraining, ytraining, xnew, rmsnoise):
    # choose RBF (Gaussian) model
    kernel1 = GPy.kern.Matern32(1)
    kernel2 = GPy.kern.rbf(1)
    kernel3 = GPy.kern.white(1)
    kernel = kernel1 + kernel3
    # build model
    xtraining = xtraining.reshape(xtraining.shape[0],1)
    ytraining = ytraining.reshape(ytraining.shape[0],1)
    xnew = xnew.reshape(xnew.shape[0],1)
    model = GPy.models.GPRegression(xtraining,ytraining,kernel, normalize_X=False)
#    model.constrain_bounded('rbf_lengthscale', 100, 300)
    model.constrain_bounded('Mat32_lengthscale', 15, 300)
    model.optimize()
    model.optimize_restarts(num_restarts = 10)
    print model
    yp, yp_var, a, b = model.predict(xnew)  # GP at xtraining points for outlier detection
    return np.array(yp.T), np.array(yp_var.T)

def makemap(data, myvmin, myvmax, xaxis, yaxis, xlines, xlabel, ylabel, title, outfile):
    xbins = data.shape[1]
    ybins = data.shape[0]
    plt.imshow(data , aspect="auto",cmap = "RdBu_r", vmin = myvmin, vmax = myvmax)    
    for i in range(xlines.shape[0]):
        plt.vlines(xlines[i]-xaxis[0],0,ybins,linestyles='dotted')

    plt.ylabel(ylabel,fontsize=16)
    plt.xlabel(xlabel,fontsize=16)
    plt.suptitle(title,fontsize=16)
    xlocs = np.arange(xbins,step = 500)
    xticklabels = []
    for i in xlocs:
        xticklabels.append(np.int(xaxis[i]))
    plt.xticks(xlocs,xticklabels,rotation="horizontal")
    
    ylocs = np.arange(ybins,step = 10)
    yticklabels = []
    for i in ylocs:
        yticklabels.append(round(yaxis[i],3))
    plt.yticks(ylocs,yticklabels)
    plt.colorbar(orientation="vertical")
    plt.savefig(outfile)
    plt.clf()

def DKD(X1, X2, theta):

    X1, X2 = np.matrix(X1), np.matrix(X2) # ensure both sets of inputs are matrices
    
    D2 = sp.distance.cdist(X1, X2, 'sqeuclidean') # calculate squared Euclidean distance
    D1 = np.zeros((X1.shape[0],X2.shape[0]))
    for i in range(X1.shape[0]):
        for j in range(X2.shape[0]):
            D1[i,j] = X1[i] - X2[j]
# This is my second derivative of the SE kernel after the two differentiation calcs
    K = (theta[0]/(theta[1]**2)) * np.exp(- D2 / (2*(theta[1]**2))) * (1 - D2/(theta[1]**2)) 

    return np.matrix(K)
