#!/bin/bash
 
#PBS -l select=1:ncpus=24:mem=1gb
#PBS -l walltime=24:00:00
#PBS -r n
#PBS -N shear_of_binary_colloidal_suspension
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

#execute gaussian IN the (fast) scratch directory
cd $SCRATCHDIR
g09 < $PBS_O_WORKDIR/$GaussianInputFilename > $GaussianOutputFilename
 
#copy files back from scratch directory
cp -r $SCRATCHDIR/* $PBS_O_WORKDIR/.
cd $PBS_O_WORKDIR
 
#print the last known statistics of the job (memory usage, cpu time, etc...)
echo >> $LOGFILE
qstat -f $PBS_JOBID >> $LOGFILE
 
echo "$PBS_JOBID ($PBS_JOBNAME) @ `hostname` at `date` in "$RUNDIR" END" >> $LOGFILE
echo "`date +"%d.%m.%Y-%T"`" >> $LOGFILE