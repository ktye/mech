/* $Id: timer.c,v 1.1.1.1 2003/09/28 22:43:47 elmar Exp $ */

/*
 * timer - get elapsed_time between calls
 */

#include <stdio.h>
#include <sys/time.h>
#include "timer.h"
double elapsed_time(void) {
	double dt;
	static int init=1;
	static long tsec, tusec;
	double new_tsec, new_tusec;
	struct timeval tp;

	if (init) {
		gettimeofday(&tp, NULL);
		tsec = tp.tv_sec;
		tusec = tp.tv_usec;
		init = 0;
		return 0.0;
	}
	gettimeofday(&tp,NULL);
	dt = (double)tp.tv_sec - (double)tsec  + (double)(tp.tv_usec - tusec)/1000000;
//	printf(">>> %f\n",(double)(tp.tv_usec - tusec)/1000000);

	tsec = tp.tv_sec;
	tusec = tp.tv_usec;
	return dt;
}
