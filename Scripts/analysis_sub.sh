#!/bin/bash
 
# number of nodes, nodes per cpu and memory of each cpu
#PBS -l select=1:ncpus=8:mem=200gb

# set execution time
#PBS -l walltime=24:00:00

# set mail alert at bginning, end and abortion of jobs
#PBS -m bea
#PBS -M aiko@thphy.uni-duesseldorf.de

# set Job name and project name
#PBS -N analysis
#PBS -A TP2SHEAR

user=`whoami`
 
#make unique scratch directory on GPFS filesystem
SCRATCHDIR=/scratch_gs/$USER/$PBS_JOBID
mkdir -p "$SCRATCHDIR"

#some (useful?) output
LOGFILE=$PBS_O_WORKDIR/$PBS_JOBNAME"."$PBS_JOBID".log"
cd $PBS_O_WORKDIR
 
echo "$PBS_JOBID ($PBS_JOBNAME) @ `hostname` at `date` in "$PBS_O_WORKDIR" START" > $LOGFILE
echo "`date +"%d.%m.%Y-%T"`" >> $LOGFILE
 
echo >> $LOGFILE
echo "GLOBAL PARAMETERS">> $LOGFILE
echo "---------------------------" >> $LOGFILE
echo "Node       : "$HOSTNAME >> $LOGFILE
echo "Arch       : "$ARCH >> $LOGFILE
echo "---------------------------" >> $LOGFILE
echo "RunDir     : "$PBS_O_WORKDIR >> $LOGFILE

# switch to the (faster) scratch directory
cd $SCRATCHDIR

# copy and execute the analysis script
cp ~/analysis_script.sh . 
./analysis_script.sh

# change back to working directory
cd $PBS_O_WORKDIR
 
#print the last known statistics of the job (memory usage, cpu time, etc...)
echo >> $LOGFILE
qstat -f $PBS_JOBID >> $LOGFILE
 
echo "$PBS_JOBID ($PBS_JOBNAME) @ `hostname` at `date` in "$RUNDIR" END" >> $LOGFILE
echo "`date +"%d.%m.%Y-%T"`" >> $LOGFILE