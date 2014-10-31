#!usr/bin/python
import os
import sys
import numpy as np

done=1

bins = 1024

# Open the folder list, go into each directory and return list of files at around 1400 MHz and their frequency

f=open('matthew.txt')

for line in f:
    pname = line.rstrip('\n')
    pname_no_J = pname[1:]
    os.chdir('./matthew_atnf_data/{0}'.format(pname))
    try:
        os.remove('mjd.txt')
    except os.error:
        pass
    
    os.system('vap -c freq *.dzTF | egrep "(1369.000|1433.000|1375.000|1368.875)" > 1400list.txt' )

# Create an array which separates the .SFTC filename and the frequency

    y = np.genfromtxt("1400list.txt", dtype=[('filename','S40'),('freq','S40')])
    numrows = y.shape[0]

    new_y = []
    new_mjd = []

# Sort into date order
    
    for i in range(numrows):
            new_y.append(y['filename'][i])
    print "new y before sort",new_y
    new_y.sort(key=lambda new_y: new_y[1:7])
    print "new y after sort",new_y
     
# Loop over each epoch (row) and perform pdv to find out how many are left after the cull of non-1024
    rem_count=0

    for i in range(numrows):
            epoch_name = new_y[i]
            os.system('pdv -Zt {0} > temp.txt '.format(epoch_name))
            stokes_line = np.genfromtxt('temp.txt', usecols=3, dtype=[('stokesI','float')], skip_header=1)
            if len(stokes_line)!=1024:
                rem_count+=1
    num_left=numrows-rem_count

# Loop over each epoch (row) and perform pdv and read into a text file           

    stokes_list = np.zeros(shape=(num_left,1024))
    b=0
    removed=0
    for i in range(numrows):

        epoch_name = new_y[i]
        os.system('pdv -Zt {0} > temp.txt '.format(epoch_name))
        stokes_line = np.genfromtxt('temp.txt', usecols=3, skip_header=1)
        if len(stokes_line) !=1024:
            removed+=1
        else:
            os.system('vap -nc "mjd, length" {0} >> mjd.txt'.format(new_y[i]))
            # os.system('vap -nc "length" {0} >> length.txt'.format(new_y[i]))
            stokes_list[b] = stokes_line
            b+=1
                
    mjdarray = np.genfromtxt("mjd.txt", dtype=[('filename2','S40'),('mjd','f8'),('len','f8')])
    for i in range(num_left):
        new_mjd.append(mjdarray['mjd'][i])
    new_mjd.sort()

    mjd_length = np.zeros((num_left,2))

    for i in range(num_left):
            mjd_length [i,0] = mjdarray['mjd'][i]
            mjd_length [i,1] = mjdarray['len'][i]
            
    mjd_length_sorted = mjd_length[mjd_length[:,0].argsort()]

    new_length = mjd_length_sorted[:,1]
    mjd_length = mjd_length_sorted[:,0]

    print "length",new_length
    print "mjd",mjd_length
                        
    print "Number of profiles removed from pulsar",pname,"=",removed,"out of",numrows,"original profiles   -   ",done,"of 180"  
    stokes_columns = np.transpose(stokes_list)
    stokes_columns2 = np.zeros((1026,num_left))
    stokes_columns2 = np.vstack([stokes_columns,new_length,new_mjd])

    np.savetxt('../../new_lists/{0}_1400list_{1}.txt'.format(pname_no_J,bins),stokes_columns2, delimiter='\t')

    t = np.genfromtxt("../../new_lists/{0}_1400list_{1}.txt".format(pname_no_J,bins))
    
    print "size",t.shape

    done+=1

    try:
        os.remove('temp.txt')
    except os.error:
        pass

    try:
        os.remove('length.txt')
    except os.error:
        pass

    try:
        os.remove('1400list.txt')
    except os.error:
        pass

    try:
        os.remove('mjd.txt')
    except os.error:
        pass
    
    os.chdir('../..')
    
f.close()


