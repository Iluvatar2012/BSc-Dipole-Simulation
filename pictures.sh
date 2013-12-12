#!/bin/bash

STEP=last

for F in Results/N*

do
	echo "Processing $F"
	./picture $F $STEP
done