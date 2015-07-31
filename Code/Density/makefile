# name the resulting program
NAME    := ../../Work/density

# include the necessary libraries and determine basic parameters
LIBS    := -lhdf5 -lm

# call the C compiler
all: density.o hdf5.o
	@echo ' '
	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Invoking: C Linker'
	$(CC) -o $(NAME) density.o hdf5.o $(LIBS)
	@echo 'Finished building: '$(NAME)
	@echo ' '

	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Cleaning up'
	rm -rf *.o
	@echo ' '

density.o: density.c
	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Invoking: C Compiler'
	$(CC) $(CFLAGS) -o density.o density.c

hdf5.o:	hdf5.c
	$(CC) $(CFLAGS) -o hdf5.o hdf5.c