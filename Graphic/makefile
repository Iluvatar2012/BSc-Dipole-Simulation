# define output name
NAME 	:= ../graphic

# define necessary libraries and basic parameters
LIBS 	:= -lm -lhdf5 -lSDL -lpng
CFLAGS	:= -c -O3 -std=c99

# set current C Compiler, maybe another one shall be used in the future
CC		:= gcc

all: simulation.o save.o
	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Invoking: GCC C Linker'
	$(CC) simulation.o save.o -o $(NAME) $(LIBS)
	@echo 'Finished building: '$(NAME)
	@echo ' '

	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Cleaning up'
	rm -rf *.o
	@echo ' '

simulation.o: simulation_gui.c
	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Invoking: GCC C Compiler'
	$(CC) $(CFLAGS) -o simulation.o simulation_gui.c
	@echo ' '

save.o: savepng.c
	$(CC) $(CFLAGS) -o save.o savepng.c
