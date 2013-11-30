# define output name
NAME 	:= graphic

# define necessary libraries and basic parameters
LIBS 	:= -lm -lhdf5 -lSDL
CFLAGS	:= -c -O3 -Wall -std=c99

all: simulation.o
	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Invoking: GCC C Linker'
	gcc simulation.o -o $(NAME) $(LIBS)
	@echo 'Finished building: '$(NAME)
	@echo ' '

	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Cleaning up'
	rm -rf *.o
	@echo ' '

simulation.o: simulation_gui.c
	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Invoking: GCC C Compiler'
	gcc $(CFLAGS) -o simulation.o simulation_gui.c
	@echo ' '