# name the resulting program
NAME    := ../../Work/analysis

# include the necessary libraries and determine basic parameters
LIBS    := -lm -lpthread -lhdf5

# call the C compiler
all: analysis.o hdf5.o functions.o
	@echo ' '
	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Invoking: C Linker'
	$(CC) -o $(NAME) analysis.o hdf5.o functions.o $(LIBS)
	@echo 'Finished building: '$(NAME)
	@echo ' '

	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Cleaning up'
	rm -rf *.o
	@echo ' '

analysis.o:	analysis.c
	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Invoking: C Compiler'
	$(CC) $(CFLAGS) -o analysis.o analysis.c

hdf5.o:	hdf5.c
	$(CC) $(CFLAGS) -o hdf5.o hdf5.c

functions.o: functions.c
	$(CC) $(CFLAGS) -o functions.o functions.c
