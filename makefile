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