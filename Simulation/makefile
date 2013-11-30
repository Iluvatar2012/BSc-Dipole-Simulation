# name the resulting program
NAME    := ../simulation

# include the necessary libraries and determine basic parameters
LIBS    := -lm -lpthread -lhdf5
CFLAGS  := -c -std=c99 -O3

# call the gnu C compiler
all: main.o hdf5.o read.o math.o
	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Invoking: GCC C Linker'
	gcc -o $(NAME) main.o hdf5.o read.o math.o $(LIBS)
	@echo 'Finished building: '$(NAME)
	@echo ' '

	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Cleaning up'
	rm -rf *.o
	@echo ' '

main.o:	main.c
	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Invoking: GCC C Compiler'
	gcc $(CFLAGS) main.c
	@echo ' '

hdf5.o:	hdf5_output.c
	gcc $(CFLAGS) -o hdf5.o hdf5_output.c
	@echo ' '

read.o: read_config.c
	gcc $(CFLAGS) -o read.o read_config.c
	@echo ' '

math.o:	extendedmath.c
	gcc $(CFLAGS) -o math.o extendedmath.c
	@echo ' '
