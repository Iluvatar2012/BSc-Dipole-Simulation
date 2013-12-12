# name the resulting program
NAME    := ../builder

# include the necessary libraries and determine basic parameters
CFLAGS  := -c -std=c99 -O3

# set the used compiler, maybe this will change in the future
CC		:= gcc

all:	main.o
	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Invoking: GCC C Linker'
	$(CC) -o $(NAME) main.o
	@echo 'Finished building: '$(NAME)
	@echo ' '

	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Cleaning up'
	rm -rf *.o
	@echo ' '

main.o: config_builder.c
	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Invoking: GCC C Compiler'
	$(CC) $(CFLAGS) -o main.o config_builder.c
	@echo ' '
