/* $Id: bin2dat.c,v 1.1 2004/04/21 13:53:05 elmar Exp $ */

/*
 * bin2dat - double binary stream to text
 */

static char id[]="$Id: bin2dat.c,v 1.1 2004/04/21 13:53:05 elmar Exp $";
#include <stdio.h>

void usage(char **argv) {
	fprintf(stderr,"usage: %s ", argv[0]);
	exit(1);
}

int main (int args, char **argv) {
	int i;
	double f;


	while (1) {
		if (fread(&f, sizeof(double), 1, stdin) != 1) {
			return 0;
		}
		printf("%g\n",f);
	}
	return 0;
}
