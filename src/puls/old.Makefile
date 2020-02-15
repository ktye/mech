INCS = -I/usr/X11R6/include
LIBS = -lm -lg2c -L/usr/X11R6/lib -lX11

puls: plot.o amoswrap.o main.o libamos.a complex.o
	cc $(INCS) $(LIBS) main.o plot.o amoswrap.o libamos.a complex.o timer.o -o puls
main.o: main.c
	cc -c $(INCS) main.c
plot.o: plot.c
	cc -c $(INCS) plot.c
complex.o: complex.c
	cc -c complex.c
ctest: ctest.c complex.o
	cc -c ctest.c
	cc -lm ctest.o complex.o -o ctest
inflexion: puls
	ln -s puls inflexion
flowfft2args: amoswrap.o libamos.a flowfft2args.o complex.o
	cc -lm -lg2c amoswrap.o libamos.a complex.o flowfft2args.o -o flowfft2args
flowfft2args.o:	flowfft2args.c
	cc -c flowfft2args.c
