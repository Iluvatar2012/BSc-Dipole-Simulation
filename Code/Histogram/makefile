# name the resulting program
NAME    := ../../Work/histogram

# include the necessary libraries and determine basic parameters
LIBS    := -lhdf5 -lm

# call the C compiler
all: histogram.o hdf5.o
	@echo ' '
	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Invoking: C Linker'
	$(CC) -o $(NAME) histogram.o hdf5.o $(LIBS)
	@echo 'Finished building: '$(NAME)
	@echo ' '

	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Cleaning up'
	rm -rf *.o
	@echo ' '

histogram.o: histogram.c
	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Invoking: C Compiler'
	$(CC) $(CFLAGS) -o histogram.o histogram.c

hdf5.o:	hdf5.c
	$(CC) $(CFLAGS) -o hdf5.o hdf5.c