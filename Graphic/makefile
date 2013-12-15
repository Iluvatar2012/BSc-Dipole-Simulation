# define output name
NAME 	:= ../graphic

# define necessary libraries
LIBS 	:= -lm -lhdf5 -lSDL

all: simulation.o
	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Invoking: GCC C Linker'
	$(CC) simulation.o -o $(NAME) $(LIBS)
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
