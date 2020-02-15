PROG=	aderamos

#copy:
#	cp $(.CURDIR)/libamos.a $(.CURDIR)/obj/

#all: copy

#SRCS=	puls.c timer.c plot.c amoswrap.c complex.c

CFLAGS+=	-I/usr/X11R6/include -I/usr/local/include
LDADD+= -L/usr/local/lib -lg2c $(HOME)/lib/complex.o $(HOME)/lib/amoswrap.o $(HOME)/lib/libamos.a -lm

.include <$(HOME)/.elmar.prog.mk>
