# name the resulting program
NAME    := ../../Work/statistics

# include the necessary libraries and determine basic parameters
LIBS    := -lm -lhdf5

all: statistics.o utilities.o hdf5.o 
	@echo ' '
	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Invoking: C Linker'
	$(CC) -o $(NAME) statistics.o hdf5.o utilities.o $(LIBS)
	@echo 'Finished building: '$(NAME)
	@echo ' '

	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Cleaning up'
	rm -rf *.o
	@echo ' '

statistics.o: statistics.c
	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Invoking: C Compiler'
	$(CC) $(CFLAGS) -o statistics.o statistics.c

utilities.o: utilities.c
	$(CC) $(CFLAGS) -o utilities.o utilities.c

hdf5.o: hdf5.c
	$(CC) $(CFLAGS) -o hdf5.o hdf5.c