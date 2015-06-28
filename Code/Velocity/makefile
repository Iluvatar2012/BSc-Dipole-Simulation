# define output name
NAME 	:= ../../Work/velocity

# define necessary libraries
LIBS 	:= -lm -lhdf5

all: velocity.o hdf5.o
	@echo ' '
	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Invoking: GCC C Linker'
	$(CC) velocity.o hdf5.o -o $(NAME) $(LIBS)
	@echo 'Finished building: '$(NAME)
	@echo ' '

	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Cleaning up'
	rm -rf *.o
	@echo ' '

velocity.o: velocity.c
	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Invoking: GCC C Compiler'
	$(CC) $(CFLAGS) -o velocity.o velocity.c

hdf5.o: hdf5.c
	$(CC) $(CFLAGS) -o hdf5.o hdf5.c
