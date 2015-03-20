import math
import matplotlib.pyplot as plt
import pylab
import numpy as np
import scipy as sc
import decimal as dc
import random as rn
from scipy import stats
import operator as op

# FREQUENTLY CHANGED VARIABLES
####################################################
    
secs_bet_nucheck = dc.Decimal('3600.0') # seconds between nu value being updated. more often is more accurate
days_bet_test = dc.Decimal('1.0') # days between the mode having the opportunity to switch
secs_bet_test = dc.Decimal('3600.0') # days between the mode having the opportunity to switch
days_bet_measures = dc.Decimal('15.0') # days between measurement of nu dot and mode fraction
days_of_obs = dc.Decimal('800.0') # days of overall duration of experiment

test_fraction = 0.7
nudot_2_fact = dc.Decimal('0.1')

####################################################


# INFREQUENTLY CHANGED VARIABLES
####################################################

nu_init = (dc.Decimal('1'))
nudot_init = dc.Decimal('-1.05e-13')
nudot_test = dc.Decimal('-1.05e-13')
start_mjd = dc.Decimal('55000.0')
####################################################


# DEFINITIONS
####################################################

half = dc.Decimal('0.5')
two = dc.Decimal('2.0')
one = dc.Decimal('1.0')
four = dc.Decimal('4.0')
secs_in_day = dc.Decimal('86400.0')
#secs_bet_test = days_bet_test*secs_in_day
secs_bet_measures = days_bet_measures*secs_in_day
nudot_1 = nudot_init
nudot_2 = nudot_1 + nudot_1*nudot_2_fact
phi = []
phi_test = []
nu = []
nudot = []
toa_no_noise = []
toa = []
toa_no_noise_days = []
residuals = []
toa_days = []
phi_last = []
noise_sd = 1e-4
rate = []
rate.append(test_fraction)
####################################################

dc.setcontext(dc.Context(prec=100))

print "--------------------------------------------------"
print "--------------------------------------------------"

print "Initial frequency:",nu_init
print "nudot of mode 1:",nudot_1
print "nudot of mode 2:",nudot_2
print "Days of observation:",days_of_obs
print "Days between measurements:",days_bet_measures
print "Days between potential mode switch:",days_bet_test

print "--------------------------------------------------"
print "--------------------------------------------------"

# First TOA at 0
phi.append(dc.Decimal(0.0))
phi_test.append(dc.Decimal(0.0))
toa_no_noise.append(dc.Decimal(0.0))
toa_no_noise_days.append(toa_no_noise[-1]/secs_in_day)

# First value of nu and nudot
nu.append(dc.Decimal(nu_init))
nudot.append(dc.Decimal(nudot_init))
g = 0
length_nu = 0
# Calculate TOAs at the interval wanted

# j is the index of the output TOA
for j in range(days_of_obs/days_bet_measures):
    count_nudot_tests = 0
    count_nudot1 = 0
# i is the index of time intervals within a TOA interval
    for i in range (int(secs_bet_measures/secs_bet_nucheck)):
        
#        print i,'out of',int(secs_bet_measures/secs_bet_nucheck)

        length_nu+=1
        g+=1

        phi.append(phi[-1] + nu[-1]*(secs_bet_nucheck) + half*nudot[-1]*(secs_bet_nucheck)**two)
#        print 'phase 1:', phi[-1]
        phi_test.append(phi_test[-1] + nu_init*(secs_bet_nucheck) + half*nudot_test*(secs_bet_nucheck)**two)

# Update nu
        nu.append(nu[-1]+nudot[-1]*secs_bet_nucheck)
# Check if update to nudot is necessary
        test_rand = np.random.uniform(0,1,1)

        if (np.mod(secs_bet_nucheck*(i+1), secs_bet_test)==0):
            count_nudot_tests += 1.0
#            if test_rand[0] > test_fraction - 0.1 * np.sin(j/24.*2.*3.1415):
#            if test_rand[0] > test_fraction - 0.1 * j:
            if test_rand[0] > test_fraction - np.random.normal(0,0.1,1):
                count_nudot1 +=1.0
                nudot.append(nudot_2)
                print 'Changed nudot'
            else:
                nudot.append(nudot_1)

# end of i for loop
# now get x toas at this moment = time + (1-phase_decimal)/nu

    rate.append(count_nudot1/(count_nudot_tests))
    num_toas = 5
    for k in range(-num_toas, 0):
        time = (j+1) * secs_bet_measures + (k+1)*secs_bet_nucheck 
        phase_decimal = dc.Decimal(phi[k]-int(phi[k]))
        print 'phase: ', phase_decimal, phi[k],time,k
        phase_test_decimal = dc.Decimal(phi_test[k]-int(phi_test[k]))
#        print 'Time in seconds', time, phase_decimal
        toa_no_noise.append(time + (one-phase_decimal)/nu[k])
        model_toa = time + (one-phase_test_decimal)/nu_init
        toa_no_noise_days.append(toa_no_noise[-1]/secs_in_day)
        residuals.append(toa_no_noise_days[-1] - model_toa/secs_in_day)

    print j

#obs_code = [' @']*len(toa_no_noise)
#fakefile = [' mode_code.rf']*len(toa_no_noise)
obs_freq = [430.000]*len(toa_no_noise)
uncert = [100.00000]*len(toa_no_noise)
noise = np.random.normal(0,1e-4,len(toa_no_noise))
for i in range(len(toa_no_noise)):
    toa_no_noise[i] += dc.Decimal(noise[i])
width = [10.00000]*len(toa_no_noise)

for i in range(len(toa_no_noise)):
#    toa_days[i] = toa_days[i]+start_mjd
    toa_no_noise_days[i] = toa_no_noise_days[i]+start_mjd

np.savetxt('./sim_toas_no_noise.txt',toa_no_noise)
np.savetxt('./residuals.txt',residuals)
np.savetxt('./rates.txt',rate)
#np.savetxt('./sim_toas_no_noise_days.txt',toa_no_noise_days)
#np.savetxt('./sim_toas_no_noise_days.txt',np.vstack((obs_code,obs_freq,toa_no_noise_days,width,uncert)).T,fmt=['%s              ','%.3f ','%.14f ','%.4f              ','%.5f'])
with open('./sim_toas_no_noise_days.txt', 'wb') as f:
  f.write(b'FORMAT 1\n')
  np.savetxt(f,np.vstack((obs_freq,toa_no_noise_days,uncert)).T,fmt=[' mode_code.rf %.3f','%.14f ','%.4f @'])
 


    
