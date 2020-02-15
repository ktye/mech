PROG=	flowfft2args

SRCS=	flowfft2args.c amoswrap.c complex.c

#copy:
#	cp $(.CURDIR)/libamos.a $(.CURDIR)/obj/

#all: copy

#CFLAGS+=	-I/usr/X11R6/include -I/usr/local/include
LDADD+= -lm -lg2c $(HOME)/lib/libamos.a

.include <$(HOME)/.elmar.prog.mk>
