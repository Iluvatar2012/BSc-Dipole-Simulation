# set standard compiler to Intel C Compiler and set Compiler Flags, export variables
CC		:= icc
CFLAGS	:= -c -O3 -std=c99 -g
export


# call make in subdirectory BSC Dipole Simulation and BSC Graphic Output
all:
	$(MAKE) -C Simulation_Periodic/ all

	@echo ' '
	@echo ' '
	@echo ' '
	@echo ' '
	@echo ' '
	$(MAKE) -C Simulation_Walls/ all

	@echo ' '
	@echo ' '
	@echo ' '
	@echo ' '
	@echo ' '
	$(MAKE) -C Simulation_Band/ all

	@echo ' '
	@echo ' '
	@echo ' '
	@echo ' '
	@echo ' '
	$(MAKE) -C Graphic/ all

	@echo ' '
	@echo ' '
	@echo ' '
	@echo ' '
	@echo ' '
	$(MAKE) -C Builder/ all

	@echo ' '
	@echo ' '
	@echo ' '
	@echo ' '
	@echo ' '
	$(MAKE) -C Pictures/ all	

	@echo ' '
	@echo ' '
	@echo ' '
	@echo ' '
	@echo ' '
	$(MAKE) -C Analysis/ all

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
	$(MAKE) -C Statistics/ all

	@echo ' '
	@echo ' '
	@echo ' '
	@echo ' '
	@echo ' '
	$(MAKE) -C Write_to_File/ all