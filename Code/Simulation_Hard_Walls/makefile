# name the resulting program
NAME    := ../../Work/simulation_with_walls

# include the necessary libraries and determine basic parameters
LIBS    := -lm -lpthread -lhdf5

# call the C compiler
all: main.o hdf5.o read.o math.o simulation.o
	@echo ' '
	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Invoking: C Linker'
	$(CC) -o $(NAME) main.o hdf5.o read.o math.o simulation.o $(LIBS)
	@echo 'Finished building: '$(NAME)
	@echo ' '

	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Cleaning up'
	rm -rf *.o
	@echo ' '

main.o:	main.c
	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Invoking: C Compiler'
	$(CC) $(CFLAGS) main.c

hdf5.o:	hdf5_output.c
	$(CC) $(CFLAGS) -o hdf5.o hdf5_output.c

read.o: read_config.c
	$(CC) $(CFLAGS) -o read.o read_config.c

math.o:	extendedmath.c
	$(CC) $(CFLAGS) -o math.o extendedmath.c

simulation.o: simulation.c
	$(CC) $(CFLAGS) -o simulation.o simulation.c
