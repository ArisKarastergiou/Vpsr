#!/bin/bash
# for i in `ls -1 ~/profile_code/data_grids/*1024.txt`
# i='/home/paul/profile_code/data_grids/0738-4042_1400list_1024.txt' 
# do 
    # pulsar=${i:35:9}; echo $pulsar;
    # rm -rf J$pulsar/
    # ./Valign.py -f $i -p J$pulsar -r 2 -d -a -o -g -b &
# done
# wait
# for i in `ls -1 ~/profile_code/data_grids/*1024.txt`
# i='/home/paul/profile_code/data_grids/0738-4042_1400list_1024.txt' 
# do
    # pulsar=${i:35:9}; echo $pulsar;
    # for j in `ls -1 J$pulsar/zoomed*.txt`
    # do
	# ./Vgp.py -f $j -p J$pulsar -i 1 &
    # done
# done
wait
for i in `ls -1 ~/profile_code/data_grids/*1024*txt`
# i='/home/paul/profile_code/data_grids/0738-4042_1400list_1024.txt' 
do 
    pulsar=${i:35:9}; echo $pulsar;
    ./Vplot.py -p J$pulsar &
done