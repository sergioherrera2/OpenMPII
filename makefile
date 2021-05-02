LIBS 	     = -lm -fopenmp
NAMEEXEC   = openmpi2
CC	        = gcc

all: $(NAMEEXEC)

$(NAMEEXEC): *.c 
	$(CC) -o $(NAMEEXEC) *.c $(LIBS) 
clean:
	rm -f *.o $(NAMEEXEC) *~ *.*~

