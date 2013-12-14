# define output name
NAME 	:= ../picture

# define necessary libraries and basic parameters
LIBS 	:= -lm -lhdf5 -lSDL
CFLAGS	:= -c -O3 -std=c99

# set current C Compiler, maybe another one shall be used in the future
CC		:= icc

all: picture.o
	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Invoking: GCC C Linker'
	$(CC) picture.o -o $(NAME) $(LIBS)
	@echo 'Finished building: '$(NAME)
	@echo ' '

	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Cleaning up'
	rm -rf *.o
	@echo ' '

picture.o: picturemaker.c
	@echo '--------------------------------------------------------------------------------------------------------------------------------------'
	@echo 'Invoking: GCC C Compiler'
	$(CC) $(CFLAGS) -o picture.o picturemaker.c
	@echo ' '
