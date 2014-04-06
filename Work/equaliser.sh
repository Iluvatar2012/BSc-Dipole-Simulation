#!/bin/bash

start_eqq="Equalised/Eqq_m_1.0__G_"
end_eqq=".hdf5"

start_res="Results/m_1.0/D_rat_1.70__GammaA_"
end_res=".0__shear_0.hdf5"

for ((i=5; i<=50; i=i+5))

do
	arg_eqq=$start_eqq$i$end_eqq
	arg_res=$start_res$i$end_res

	./equaliser $arg_res $arg_eqq

done
