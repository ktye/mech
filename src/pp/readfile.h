/* $Id: readfile.h,v 1.2 2004/04/02 19:25:13 elmar Exp $ */

#include <time.h>
struct datafile {
	int endian;
	int rows;
	int cols;
	char xlabel[100], ylabel[100], zlabel[100];
	double xmin, xmax, ymin, ymax;
	char header[1024*10];	// continous string
	int hlen[10];		// lengths of lines
	char filename[512];
	double *data;
};

void readfile( struct datafile *, const char *);

