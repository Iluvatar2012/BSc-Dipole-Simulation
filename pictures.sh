#!/bin/bash

STEP=last

for F in Results/varying_m/N*

do
	echo "Processing $F"
	./picture $F $STEP
done
