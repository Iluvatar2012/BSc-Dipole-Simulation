#!/bin/bash

# copy all done simulation data
cp /scratch_gs/aiber100/51583\[*\].hpc-batch/Results/* . &
cp /scratch_gs/aiber100/51584\[*\].hpc-batch/Results/* .

# analyse all sets
for F in .
do
	~/Work/analysis $F &
done