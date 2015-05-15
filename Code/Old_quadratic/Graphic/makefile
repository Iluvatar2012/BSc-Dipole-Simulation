# define output name
NAME 	:= ../../Work/graphic

# define necessary libraries
LIBS 	:= -lm -lhdf5 -lSDL

all: graphic.o hdf5.o picture.o
	@echo ' '
	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Invoking: GCC C Linker'
	$(CC) graphic.o hdf5.o picture.o -o $(NAME) $(LIBS)
	@echo 'Finished building: '$(NAME)
	@echo ' '

	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Cleaning up'
	rm -rf *.o
	@echo ' '

graphic.o: graphic.c
	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Invoking: GCC C Compiler'
	$(CC) $(CFLAGS) -o graphic.o graphic.c

hdf5.o: hdf5.c
	$(CC) $(CFLAGS) -o hdf5.o hdf5.c

picture.o: picture.c
	$(CC) $(CFLAGS) -o picture.o picture.c
