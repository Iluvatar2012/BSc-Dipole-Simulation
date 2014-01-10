# set standard compiler to Intel C Compiler and set Compiler Flags, export variables
CC		:= icc
CFLAGS	:= -c -O3 -std=c99
export


# call make in subdirectory BSC Dipole Simulation and BSC Graphic Output
all:
	$(MAKE) -C Simulation/ all

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