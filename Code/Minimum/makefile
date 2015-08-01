# name the resulting program
NAME    := ../../Work/minimum

# include the necessary libraries and determine basic parameters
LIBS    := -lhdf5 -lm

# call the C compiler
all: minimum.o hdf5.o
	@echo ' '
	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Invoking: C Linker'
	$(CC) -o $(NAME) minimum.o hdf5.o $(LIBS)
	@echo 'Finished building: '$(NAME)
	@echo ' '

	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Cleaning up'
	rm -rf *.o
	@echo ' '

minimum.o: minimum.c
	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Invoking: C Compiler'
	$(CC) $(CFLAGS) -o minimum.o minimum.c

hdf5.o:	hdf5.c
	$(CC) $(CFLAGS) -o hdf5.o hdf5.c