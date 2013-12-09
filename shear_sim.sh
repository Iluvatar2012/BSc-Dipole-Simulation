#!/bin/bash

# name the simulation
#$ -N bin_shear

# get status email notifications
#$ -M aiko@thphy.uni-duesseldorf.de
#$ -m bea

# work in current directory with current bash
#$ -cwd
#$ -S /bin/bash

# throw stderr in a special file
#$ -e ./job.err

# run parallel job, 8 threads on parallel.q or normal.q
#$ -pe shm 8
#$ -q parallel.q,normal.q

# iterate over the given number of config files, in this case 51
#$ -t 1-26:1

# start the actual simulation
source /export/apps/intel/bin/compilervars.sh intel64
./simulation Config_Files/$SGE_TASK_ID
