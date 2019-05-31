obs = makePopulation.o G.o develop.o utilities.o
CFLAGS = -c -I/usr/local/include -g -Os -std=c99
LFLAGS = -lgsl -lm -lgslcblas -L/usr/local/lib

makePopulation : $(obs)
	$gcc $(LFLAGS) $(obs) -o makePopulation

makePopulation.o : makePopulation.c constants.h develop.h G.h
	$gcc $(CFLAGS) makePopulation.c

G.o : G.c G.h constants.h utilities.h develop.h 
	$gcc $(CFLAGS) G.c develop.c

develop.o : develop.c develop.h G.h constants.h
	$gcc $(CFLAGS) develop.c

utilities.o : utilities.c develop.h
	$gcc $(CFLAGS) utilities.c
