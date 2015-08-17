# A growing collection of functions for V code
import numpy as np
import GPy
import matplotlib
import matplotlib.pyplot as plt
import scipy.spatial as sp
import math
import os
import scipy.signal as ss
import sys
def SE(X1, X2, amp, length, white_noise = False):
        '''
        Squared exponential covariance function
        '''
        X1, X2 = np.matrix(X1), np.matrix(X2) # ensure both sets of inputs are matrices
        D2 = sp.distance.cdist(X1, X2, 'sqeuclidean') # calculate squared Euclidean distance
        K = amp**2 * np.exp(- D2 / (2*(length**2))) # calculate covariance matrix
        return np.matrix(K)

def forpaper(pulsar,profile1, profile2, median, yllim=None, yulim=None):
    nbins = profile1.shape[0]
    xaxis = np.linspace(0,float(nbins)/1024,nbins)
    xlocs = np.linspace(0,nbins-1,10)
    xticklabels = []
    for i in xlocs:
	    xticklabels.append(np.round(xaxis[int(i)],3))
    ax = plt.figure().add_subplot(1,1,1)
    plt.plot(profile1,'b')
    plt.plot(profile2,'r')
    plt.plot(median,'k--')
    plt.vlines(np.argmax(median),yllim,yulim,linestyles='dotted')
    plt.ylim(yllim,yulim)
    plt.ylabel(r'Normalised Intensity',fontsize=14)
    plt.xlabel(r'Fraction of Pulse Period',fontsize=14)
    plt.xticks(xlocs,xticklabels)
    plt.savefig('./{0}/{0}_profiles.png' .format(pulsar))
    plt.clf()


def makeplots(pulsar, data, mjd, dir, allbins, template=None, yllim=None, yulim=None, peakindex=None,cal=None):
    nbins = data.shape[0]
    nprofiles = data.shape[1]
    xaxis = np.linspace(0,float(nbins)/allbins,nbins)
    xlocs = np.linspace(0,nbins-1,10)
    xticklabels = []
    for i in xlocs:
	    xticklabels.append(np.round(xaxis[int(i)],3))
    if not (os.path.exists('./{0}/{1}'.format(pulsar,dir))):
        os.mkdir('./{0}/{1}'.format(pulsar,dir))
    for i in range(nprofiles):
        plt.plot(data[:,i])
        if template!=None:
            plt.plot(template,'r')
	    plt.vlines(np.argmax(template),yllim,yulim,linestyles='dotted')
        if yllim !=None or yulim !=None:
            plt.ylim(yllim,yulim)
	if cal == None:
		plt.ylabel(r'Normalised Flux Density',fontsize=14)
#		plt.yticks([])
	else:
		plt.ylabel(r'Flux Density (mJy)',fontsize=14)
	plt.xlim(0,data.shape[0])
	plt.xlabel(r'Fraction of Pulse Period',fontsize=14)
	plt.xticks(xlocs,xticklabels)
#	plt.suptitle('{0}'.format(mjd[i]), fontsize=14, fontweight='bold')
        plt.savefig('./{0}/{1}/{2}_{3}.png' .format(pulsar,dir,int(math.floor(mjd[i])),i))
        plt.clf()




def goodplots_ip(pulsar, data_mp, data_ip, mjd, dir, allbins, startbin, template_mp, template_ip, yllim, yulim, peakindex, cal=None):
    if not (os.path.exists('./{0}/{1}'.format(pulsar,dir))):
        os.mkdir('./{0}/{1}'.format(pulsar,dir))
    nprofiles = data_mp.shape[1]

    nbins_mp = data_mp.shape[0]
    xaxis_mp = np.linspace(0,float(nbins_mp)/allbins,nbins_mp)
    xlocs_mp = np.linspace(0,nbins_mp-1,5)
    xticklabels_mp = []
    for i in xlocs_mp:
	    xticklabels_mp.append(np.round(xaxis_mp[int(i)],3))

    nbins_ip = data_ip.shape[0]
    xaxis_ip = np.linspace(float(startbin)/allbins,float(nbins_ip+startbin)/allbins,nbins_ip)
    xlocs_ip = np.linspace(0,nbins_ip-1,5)
    xticklabels_ip = []
    for i in xlocs_ip:
	    xticklabels_ip.append(np.round(xaxis_ip[int(i)],3))


    for i in range(nprofiles):
        fig=plt.figure()
        ax_mp=fig.add_subplot(1,2,1)
        plt.plot(data_mp[:,i])
        plt.plot(template_mp,'r')
	plt.xticks(xlocs_mp,xticklabels_mp)
        plt.ylim(yllim,yulim)
        plt.vlines(np.argmax(template_mp),yllim,yulim,linestyles='dotted')
	plt.xlabel(r'Fraction of Pulse Period',fontsize=14)
        if cal == None:
                plt.ylabel(r'Normalised Flux Density',fontsize=14)
                plt.yticks([])
        else:
                plt.ylabel(r'Flux Density (mJy)',fontsize=14)
        ax_ip=fig.add_subplot(1,2,2)
        plt.plot(data_ip[:,i])
        plt.plot(template_ip,'r')
	plt.xticks(xlocs_ip,xticklabels_ip)
        plt.ylim(yllim,yulim)
	plt.xlabel(r'Fraction of Pulse Period',fontsize=14)
	plt.yticks([])

        #plt.suptitle('{0}'.format(mjd[i]), fontsize=14, fontweight='bold')
        plt.savefig('./{0}/{1}/{2}_{3}.png' .format(pulsar,dir,int(math.floor(mjd[i])),i))
        plt.clf()

        plt.close(fig)

def aligndata(baselineremoved, brightest, pulsar):
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
    plt.savefig('./{0}/{0}_brightest.png' .format(pulsar))
    plt.clf()
    for i in range(nprofiles):
        xcorr = np.correlate(template,baselineremoved[:,i],"full")
        lag = np.argmax(xcorr)
        aligned[:,i] = np.roll(baselineremoved[:,i],lag)
    template = np.median(aligned,1)
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
        newtemplate = np.median(aligned,1)
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
    peak = np.zeros(nprofiles)
    for i in range(nprofiles):
        rms = np.zeros(8)
        mean = np.zeros(8)
        section = nbins/8
        for j in range(8):
            rms[j] = np.std(data[j*section:(j+1)*section,i])
            mean[j] = np.mean(data[j*section:(j+1)*section,i])
        smallestrms[i] = np.min(rms)
        peak[i] = np.max(data[:,i]) # remove low snr not just the std of noise
        baseindex = np.argmin(rms)
        baseline = mean[baseindex]
        smallestmean[i] = baseline
        baselineremoved[:,i] = data[:,i] - baseline
    medianrms = np.median(smallestrms)
    medianpeak = np.median(peak)
    outlierindex = []
    inlierindex = []
    for i in range(nprofiles):
	    if smallestrms[i]/np.max(data[:,i]) > outliers * medianrms/medianpeak:
#        if smallestrms[i] > outliers * medianrms or smallestrms[i] * outliers < medianrms or np.min(baselineremoved[:,i]) < - 5 * smallestrms[i]:
		    outlierindex.append(i)
	    else:
		    inlierindex.append(i)

    ou = np.array(outlierindex)
    inl = np.array(inlierindex)
    

    
    removedprofiles = np.delete(baselineremoved,inl,1)
    baselineoutlierremoved = np.delete(baselineremoved, ou, 1)
    rmsremoved = np.delete(smallestrms, ou)
    print 'Removed outliers: ', ou.shape[0]
    return baselineoutlierremoved, removedprofiles, rmsremoved, ou, inl


def findbrightestprofile(data,rmsdata):
    snr = np.zeros(rmsdata.shape[0])
    for i in range(data.shape[1]):
        snr[i] = np.max(data[:,i])/rmsdata[i]
    brightestindex = np.argmax(snr)
    return brightestindex

def binstartend(data,peakoriginal,rms):
    peak = np.max(data)
    peakbin = np.argmax(data)
    bins = data.shape[0]
    window = bins/20
 #   print peak, peakbin, rms
 #   print 'peak snr is',peak/rms
    peaksnr = peak/rms
    power = peaksnr
    peaks = 0
    lstart = 0
    lend = 0
    if peaksnr > 15:
        peaks = 1
        thisbin = peakbin - 1
        while power > 0.02*peakoriginal:
            if peakbin != 0 :
                power = data[thisbin]
                data[thisbin] = 0
                thisbin = np.max((thisbin - window,0))
        lstart = thisbin
        thisbin = peakbin
        power = peaksnr
        while power > 0.02*peakoriginal:
            if peakbin != bins - 1 :
                power = data[thisbin]
                data[thisbin] = 0
                thisbin = np.min((thisbin + window,bins-1))
        lend = thisbin
#    start = np.max((lstart - int(bins/50), 0))
#    end = np.min((lend + int(bins/50), bins - 1))
    start = lstart
    end = lend
    data[start:end] = 0
    cuttemplate = np.array(data)
#    print "start and end",start+70, end+50
    return start, end, peaks, cuttemplate


def gpinferred(xtraining, ytraining, xnew):
    # choose RBF (Gaussian) model
    kernel1 = GPy.kern.Matern32(1)
    kernel2 = GPy.kern.RBF(1)
    kernel3 = GPy.kern.White(1)
    kernel = kernel1 
#   + kernel3
    # build model
    xtraining = xtraining.reshape(xtraining.shape[0],1)
    ytraining = ytraining.reshape(ytraining.shape[0],1)
    xnew = xnew.reshape(xnew.shape[0],1)
    model = GPy.models.GPRegression(xtraining,ytraining,kernel)
    model.Mat32.lengthscale.constrain_bounded(30,300)
    model.optimize()
    model.optimize_restarts(num_restarts = 10)
    print model
    yp, yp_var = model.predict(xnew)  # GP at xtraining points for outlier detection
    return np.array(yp.T), np.array(yp_var.T), model

def makemap(data, myvmin, myvmax, xaxis, yaxis, xlines, xlabel, ylabel, title, outfile, peakline=None, combined=None):
    if combined == None:
	fig=plt.figure()
        fig.set_size_inches(16,10)
    xbins = data.shape[1]
    ybins = data.shape[0]
    plt.imshow(data , aspect="auto",cmap = "RdBu_r", vmin = myvmin, vmax = myvmax)    
    for i in range(xlines.shape[0]):
        plt.vlines(xlines[i]-xaxis[0],0,ybins,linestyles='dotted')
#    if xlinesremoved.shape[0]!=0:
#       for i in range(xlinesremoved.shape[0]):
#            if xlinesremoved[i]-xaxis[0] >= 0:
#                plt.vlines(xlinesremoved[i]-xaxis[0],0,ybins,linestyles='dotted', color = "r")

    plt.ylabel(ylabel,fontsize=16)
    plt.xlabel(xlabel,fontsize=16)
#    plt.suptitle(title,fontsize=16)
    xlocs = np.arange(xbins,step = 500)
    xticklabels = []
    for i in xlocs:
        xticklabels.append(np.int(xaxis[i]))
    plt.xticks(xlocs,xticklabels,rotation="horizontal")

    if peakline!=None:
        plt.hlines(peakline,0,xbins)
    
#    ylocs = np.arange(ybins,step = 30)
    ylocs = np.linspace(0,ybins,10)
    yticklabels = []
#    for i in ylocs:
#        print "KEK",yaxis.shape
#       yticklabels.append(round(yaxis[i],3))
#    plt.yticks(ylocs,yticklabels)
    if combined == None:
        plt.colorbar(orientation="vertical")
        plt.savefig(outfile)
        plt.clf()

def DKD(X1, X2, theta):

    X1, X2 = np.matrix(X1), np.matrix(X2) # ensure both sets of inputs are matrices
    
    D2 = sp.distance.cdist(X1, X2, 'sqeuclidean') # calculate squared Euclidean distance
    D1 = np.zeros((X1.shape[0],X2.shape[0]))
#    for i in range(X1.shape[0]):
#        for j in range(X2.shape[0]):
#            D1[i,j] = X1[i] - X2[j]
# This is my second derivative of the SE kernel after the two differentiation calcs
    K = theta[0] * np.exp(- D2 / (2*(theta[1]**2))) * ( theta[1]**2 - D2) / theta[1]**4
    return np.matrix(K)

def TKD(X1, X2, theta):

    X1, X2 = np.matrix(X1), np.matrix(X2) # ensure both sets of inputs are matrices
    
    D2 = sp.distance.cdist(X1, X2, 'sqeuclidean') # calculate squared Euclidean distance
    D1 = np.zeros((X1.shape[0],X2.shape[0]))
    for i in range(X1.shape[0]):
        for j in range(X2.shape[0]):
            D1[i,j] = -(X1[i] - X2[j])
# This is my second derivative of the SE kernel after the two differentiation calcs
    K = theta[0] * D1 * np.exp(- D2 / (2*(theta[1]**2))) * ( 3 * theta[1]**2 - D2) / theta[1]**6
    return np.matrix(K)

# Normalises the data to the peak

def norm_to_peak(data,onpulsedata):
	
	

    profiles = data.shape[1]

    for i in range(profiles):
#            data[:,i]/=np.mean(data[peakbin-5:peakbin+5,i])
	    data[:,i]/=np.mean(onpulsedata[:,i])
    return data
# Produces a plot with non-normalised variabiliy map, normalised variability map and spindown rate

def combined_map(zoom_region_no,data_norm,data_not, myvmin_norm, myvmax_norm, myvmin_not,myvmax_not, xaxis, yaxis, xlines_norm, xlines_not, nudot, errorbars, mjdinferspindown, xlabel, ylabel,pulsar, peakline=None):

    import matplotlib as mpl
    label_size = 14
    mpl.rcParams['xtick.labelsize'] = label_size
    mpl.rcParams['ytick.labelsize'] = label_size
	

    fig=plt.figure()
    fig.set_size_inches(16,10)

    ax=fig.add_subplot(3,1,1)
    ax.xaxis.set_visible(False)

    xbins = data_not.shape[1]
    ybins = data_not.shape[0]
    plt.imshow(data_not , aspect="auto",cmap = "RdBu_r", vmin = myvmin_not, vmax = myvmax_not)    

    for i in range(xlines_norm.shape[0]):
	    plt.vlines(xlines_norm[i]-xaxis[0],0,ybins,linestyles='dotted')

#    if xlinesremoved.shape[0]!=0:
#       for i in range(xlinesremoved.shape[0]):
#            if xlinesremoved[i]-xaxis[0] >= 0:
#                plt.vlines(xlinesremoved[i]-xaxis[0],0,ybins,linestyles='dotted', color = "r")

#    plt.ylabel(ylabel,fontsize=18)
    plt.xlabel(xlabel,fontsize=18)
#    plt.suptitle(title,fontsize=16)
    xlocs = np.arange(xbins,step = 500)
    xticklabels = []
    for i in xlocs:
        xticklabels.append(np.int(xaxis[i]))
    plt.xticks(xlocs,xticklabels,rotation="horizontal")

    if peakline!=None:
        plt.hlines(peakline,0,xbins)
    
#    ylocs = np.arange(ybins,step = int(data_norm.shape[0]/4.1))

    ylocs = np.linspace(0,ybins-1,num = 5)
    yticklabels = []
    print yaxis.shape, ylocs.shape
    for i in ylocs:
        yticklabels.append(round(yaxis[i],2))
    plt.yticks(ylocs,yticklabels)


    cbaxes1 = fig.add_axes([0.6, 0.645, 0.30, 0.01])
    cb1 = plt.colorbar(cax = cbaxes1,orientation="horizontal")
    cb1.update_ticks()
    cb1.ax.tick_params(labelsize=12)
    ax=fig.add_subplot(3,1,2)
    ax.xaxis.set_visible(False)

    xbins = data_norm.shape[1]
    ybins = data_norm.shape[0]
    plt.imshow(data_norm , aspect="auto",cmap = "RdBu_r", vmin = myvmin_norm, vmax = myvmax_norm, zorder = -5)    
    
    for i in range(xlines_norm.shape[0]):
	    plt.vlines(xlines_norm[i]-xaxis[0],0,ybins,linestyles='dotted')
#    if xlinesremoved.shape[0]!=0:
#       for i in range(xlinesremoved.shape[0]):
#            if xlinesremoved[i]-xaxis[0] >= 0:
#                plt.vlines(xlinesremoved[i]-xaxis[0],0,ybins,linestyles='dotted', color = "r")

#    plt.ylabel(ylabel,fontsize=18)
    plt.xlabel(xlabel,fontsize=18)
 #   plt.suptitle(title,fontsize=16)
    xlocs = np.arange(xbins,step = 500)
    xticklabels = []
    for i in xlocs:
        xticklabels.append(np.int(xaxis[i]))
    plt.xticks(xlocs,xticklabels,rotation="horizontal")

    if peakline!=None:
        plt.hlines(peakline,0,xbins)
    
#    ylocs = np.arange(ybins,step = int(data_norm.shape[0]/4.1))
    ylocs = np.linspace(0,ybins-1,num = 5) 
    yticklabels = []
    for i in ylocs:
        yticklabels.append(round(yaxis[i],2))
    plt.yticks(ylocs,yticklabels)


    ax=fig.add_subplot(3,1,3)

    power = int((-1)*np.floor(math.log10(abs(np.median(nudot)))))

    nudot = nudot*10**power
    #spinllim = spinllim*10**power
    #spinulim = spinulim*10**power

    errorbars = errorbars*10**power

    ax.grid()
#    plt.plot(mjdinferspindown, nudot)
#    plt.errorbar(mjdinferspindown, nudot, yerr=errorbars,linestyle='-')
    plt.errorbar(mjdinferspindown[3:-3], nudot[3:-3], yerr=errorbars[3:-3],linestyle='-')
#    plt.fill_between(mjdinferspindown, spinllim, spinulim, color = 'b', alpha = 0.2)
    plt.xlim(xaxis[0],xaxis[-1])
    #start_mjd_diff = int(abs(xaxis[0]-mjdinferspindown[0]))
    #end_mjd_diff = int(abs(xaxis[-1]-mjdinferspindown[-1]))
    #if end_mjd_diff == 0:
	   # end_mjd_diff = 1
    #plt.ylim(np.median(nudot)-2*np.std(nudot),np.median(nudot)+2*np.std(nudot))

# Make the y limits based on which nudot values are visible on plot, i.e. between mjdinferspindown[3:-3]

    nudot_visible = []

    for i in range(nudot.shape[0]):
	    if i > 3 and i < nudot.shape[0]-3 and mjdinferspindown[i] > xaxis[0] and mjdinferspindown[i] < xaxis[-1]:
		    nudot_visible.append(nudot[i])
        
    diff_ymin_ymax = np.max(nudot_visible)-np.min(nudot_visible)

    plt.ylim(np.min(nudot_visible)-0.1*diff_ymin_ymax,np.max(nudot_visible)+0.05*diff_ymin_ymax)


    fig.text(0.05, 0.777, 'Fraction of\nPulse Period', ha='center', va='center', rotation='vertical', size=18)
    fig.text(0.05, 0.503, 'Fraction of\nPulse Period', ha='center', va='center', rotation='vertical', size=18)
    fig.text(0.05, 0.23, r'$\mathrm{{\dot{{\nu}}}}$ ($\mathrm{{10^{{-{0}}} s^{{-2}}}}$)'.format(power), ha='center', va='center', rotation='vertical',size=18)


    plt.xlabel(xlabel,fontsize=18)
#    plt.ylabel(r'$\mathrm{{\dot{{\nu}}}}$ ($\mathrm{{10^{{-{0}}} s^{{-2}}}}$)'.format(power) ,fontsize=18)
    ax.xaxis.labelpad = 20
    y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    plt.locator_params(axis = 'y', nbins=5)
    ax.yaxis.set_major_formatter(y_formatter)
    plt.subplots_adjust(hspace=0.1)
    cbaxes2 = fig.add_axes([0.60, 0.37, 0.3, 0.01])
    cb2 = plt.colorbar(cax = cbaxes2,orientation="horizontal")
    cb2.ax.tick_params(labelsize=12) 
    cb2.update_ticks()

    plt.savefig('./{0}/{0}_final_{1}.png'.format(pulsar,zoom_region_no))

def autocorr(rr,plotname):
	# Compute autocorrelation function of residuals to test if it is white noise
	autocorr = np.correlate(rr,rr, mode='full')
	datapoints = rr.shape[0] 
	print datapoints
	smoothac = ss.savgol_filter(autocorr[datapoints+1:],11, 3)
	previous = 1000000
	zerocrossing = 0
	i=0
	point = np.zeros(4)
	while zerocrossing < 3 and i<datapoints:
		if (smoothac[i]*previous < 0):
			zerocrossing += 1
			point[zerocrossing] = i
			print "acc an", zerocrossing, i 
		previous = smoothac[i]
		i+=1
	plt.plot(smoothac[0:10*point[3]])
	plt.savefig(plotname)
	plt.clf()
	peak = 0.5*(point[3]-point[2])+point[2]
	noise = np.std(autocorr[-50:])
	snr = smoothac[int(peak)]/noise
	return peak, snr
