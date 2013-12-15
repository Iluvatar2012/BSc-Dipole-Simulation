#!/bin/bash

STEP=last

for F in Results/steps_3000000/N*

do
	echo "Processing $F"
	./picture $F $STEP
done
