PROG=	puls

#copy:
#	cp $(.CURDIR)/libamos.a $(.CURDIR)/obj/

#all: copy

SRCS=	puls.c timer.c plot.c amoswrap.c complex.c

CFLAGS+=	-I/usr/X11R6/include -I/usr/local/include
#LDADD+= -lg2c -L/usr/X11R6/lib -lX11 -lXpm -L/usr/local/lib -ljpeg -lpng -ltiff -lz -lgif -lImlib $(HOME)/lib/libamos.a
LDADD+= -lg2c -L/usr/X11R6/lib -lX11 -lXpm -L/usr/local/lib -ljpeg -lpng -ltiff -lz -lImlib $(HOME)/lib/libamos.a

.include <$(HOME)/.elmar.prog.mk>
