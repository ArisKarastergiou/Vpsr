# Reads in aligned and normalized 2d profile data(epoch, bin), calculates median
# profile, then subtracts from all profiles. Finally, runs a Gaussian process regression
# model to infer the data on a different, probably even, grid.

#!/usr/bin/python
import argparse
import pylab as pb
pb.ion()
import numpy as np
import GPy

# Read command line arguments
parser = argparse.ArgumentParser(description='Pulsar profile variability map creator')
parser.add_argument('-f','--filename', help='File containing aligned profiles at multiple epochs', required=True)
args = vars(parser.parse_args())
