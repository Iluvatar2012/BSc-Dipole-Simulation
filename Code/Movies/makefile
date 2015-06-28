# define output name
NAME 	:= ../../Work/movies

# define necessary libraries
LIBS 	:= -lm -lhdf5 -lSDL

all: movies.o hdf5.o picture.o
	@echo ' '
	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Invoking: GCC C Linker'
	$(CC) movies.o hdf5.o picture.o -o $(NAME) $(LIBS)
	@echo 'Finished building: '$(NAME)
	@echo ' '

	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Cleaning up'
	rm -rf *.o
	@echo ' '

movies.o: movies.c
	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Invoking: GCC C Compiler'
	$(CC) $(CFLAGS) -o movies.o movies.c

hdf5.o: hdf5.c
	$(CC) $(CFLAGS) -o hdf5.o hdf5.c

picture.o: picture.c
	$(CC) $(CFLAGS) -o picture.o picture.c
