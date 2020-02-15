PROG=	pp

CFLAGS+=	`pkg-config --cflags gtk+-2.0`
LDADD+=		`pkg-config --libs gtk+-2.0` -lm

SRCS= pp.c readfile.c

.include <$(HOME)/.elmar.prog.mk>
