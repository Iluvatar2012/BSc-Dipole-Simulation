# set standard compiler to Intel C Compiler and set Compiler Flags, export variables
CC		:= icc
CFLAGS	:= -c -O3 -std=c99 -g
export


# call make in subdirectory BSC Dipole Simulation and BSC Graphic Output
all:
	$(MAKE) -C Builder/ all

	@echo ' '
	@echo ' '
	@echo ' '
	@echo ' '
	@echo ' '

	$(MAKE) -C Equaliser/ all

	@echo ' '
	@echo ' '
	@echo ' '
	@echo ' '
	@echo ' '

	$(MAKE) -C Movies/ all

	@echo ' '
	@echo ' '
	@echo ' '
	@echo ' '
	@echo ' '


	$(MAKE) -C Simulation_Band/ all
