#!/bin/bash
# script used for running multiple files through Vpsr code
# use with care, i.e. fear!
for i in `ls -1 data_grids/*1024.txt`
do
    pulsar=${i:11:9}; echo $pulsar;
    rm -rf J$pulsar/
    ./Valign.py -f $i -p J$pulsar -o -r 2 -d -g -n & 
done
wait
for i in `ls -1 data_grids/*1024.txt`
do
    pulsar=${i:11:9}; echo $pulsar;
    for j in `ls -1 J$pulsar/zoomed*.txt`
    do
	./Vgp.py -f $j -p J$pulsar -i 1 &
    done
done

Vnudot.py_2kern -f ../matthew_residuals/J1602-5100_residuals_1400.dat -e ../matthew_residuals/J1602-5100.par -p J1602-5100 -d

Vplot.py -p J1602-5100 -c
