# name the resulting program
NAME    := ../analysis

# include the necessary libraries and determine basic parameters
LIBS    := -lm -lhdf5 -lSDL

# call the gnu C compiler
all: analysis.o hdf5.o picture.o assist.o
	@echo ' '
	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Invoking: $(CC) C Linker'
	$(CC) -o $(NAME) analysis.o hdf5.o picture.o assist.o $(LIBS)
	@echo 'Finished building: '$(NAME)
	@echo ' '

	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Cleaning up'
	rm -rf *.o
	@echo ' '

analysis.o:	analysis.c
	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Invoking: GCC C Compiler'
	$(CC) $(CFLAGS) -o analysis.o analysis.c

hdf5.o:	hdf5_read.c
	$(CC) $(CFLAGS) -o hdf5.o hdf5_read.c

picture.o: picture.c
	$(CC) $(CFLAGS) -o picture.o picture.c

assist.o: assist.c
	$(CC) $(CFLAGS) -o assist.o assist.c