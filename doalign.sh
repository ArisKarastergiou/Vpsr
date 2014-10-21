#!/bin/bash
for i in `ls -1 data_grids/*1024.txt`
do
#i='data_grids/1825-0935_1400list_1024.txt' 
    pulsar=${i:11:9}; echo $pulsar;
    rm -rf J$pulsar/
    ./Valign.py -f $i -p J$pulsar -r 2 -d -a & 
done
wait
for i in `ls -1 data_grids/*1024.txt`
do
#i='data_grids/1825-0935_1400list_1024.txt' 
    pulsar=${i:11:9}; echo $pulsar;
    for j in `ls -1 J$pulsar/zoomed*.txt`
    do
	./Vgp.py -f $j -p J$pulsar -i 1 -d &
    done
done