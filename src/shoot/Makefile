all:
	@echo use mk instead
	exit 1

PROG=	shoot
LDADD+=	-lm
#CFLAGS+=
#LDFLAGS+=
SRCS=	ode.c baseflow.c v.c uw.c uvw.c shoot.c

LINKS=	shoot shootinter shoot shootsingle shoot shootplot

.include <$(HOME)/.elmar.prog.mk>
