#!/bin/bash

start_eqq="Equalised/Eqq_m_0.50__G_"
end_eqq=".hdf5"

start_res="Equalised/m_0.50__GammaA_"
end_res=".0__shear_0.hdf5"

for ((i=1; i<=30; i=i+1))

do
	arg_eqq=$start_eqq$i$end_eqq
	arg_res=$start_res$i$end_res

	./equaliser $arg_res $arg_eqq

done
