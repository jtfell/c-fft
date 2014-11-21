FLAGS = -std=c99

all: complex.c fft.c test.c
	gcc complex.c $(FLAGS) -c
	gcc fft.c $(FLAGS) -c 
	gcc test.c complex.o fft.o $(FLAGS) -o fftTest -lm

clean:
	rm complex.o fft.o fftTest
