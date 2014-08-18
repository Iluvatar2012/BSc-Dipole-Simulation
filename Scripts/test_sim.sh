#!/bin/bash

# name the simulation
#$ -N bin_shear_test

# get status email notifications
#$ -M aiko@thphy.uni-duesseldorf.de
#$ -m bea

# work in current directory with current bash
#$ -cwd
#$ -S /bin/bash

# throw stderr in a special file
#$ -e ./job.err
#$ -o ./job.out

# run parallel job, 8 threads on test.q
#$ -pe shm 8
#$ -q test.q

# start the actual simulation
source /export/apps/intel/bin/compilervars.sh intel64
./simulation Config_Files/test -1
