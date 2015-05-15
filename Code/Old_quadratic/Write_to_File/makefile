# name the resulting program
NAME    := ../../Work/write_to_file

# include the necessary libraries and determine basic parameters
LIBS    := -lhdf5

# call the C compiler
all: write.o hdf5.o
	@echo ' '
	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Invoking: C Linker'
	$(CC) -o $(NAME) write.o hdf5.o $(LIBS)
	@echo 'Finished building: '$(NAME)
	@echo ' '

	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Cleaning up'
	rm -rf *.o
	@echo ' '

write.o: write.c
	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Invoking: C Compiler'
	$(CC) $(CFLAGS) -o write.o write.c

hdf5.o:	hdf5.c
	$(CC) $(CFLAGS) -o hdf5.o hdf5.c