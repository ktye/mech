/* $Id: readfile.c,v 1.4 2004/04/06 13:54:40 elmar Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

struct datafile {
	int endian;
	int rows;
	int cols;
	char xlabel[100], ylabel[100], zlabel[100];
	double xmin, xmax, ymin, ymax;
	char header[1024*10];
	int hlen[10];
	char filename[512];
	double *data;
};

int mendian(void) {
	int x = 1;
	if(*(char *)&x == 1) return 0; // little-endian
	else return 1; // big-endian
}

void dswap(double *x) {
        union {
                double d;
                unsigned char b[8];
        } u,v;
        u.d = *x ;
        v.b[0]=u.b[7] ;
        v.b[1]=u.b[6] ;
        v.b[2]=u.b[5] ;
        v.b[3]=u.b[4] ;
        v.b[4]=u.b[3] ;
        v.b[5]=u.b[2] ;
        v.b[6]=u.b[1] ;
        v.b[7]=u.b[0] ;
        *x=v.d ;
}

void readfile( struct datafile *f, const char *filename ) {
	int i,j,space;
	FILE *fp;
	size_t len;
	char *lp;
	char *s;
	int c;
	int linelength;
	int start, end;
	char endian[15];
	double x;
	fp = fopen(filename, "r");
	if (fp == NULL) {
		fprintf(stderr,"error: cannot open file %s\n", filename);
		exit(1);
	}

	snprintf(f->filename,512,"%s",filename);

	j=0;
	i=0;
	linelength=0;
	while(i<10*1024-1) {
		c = fgetc(fp);
		f->header[i] = c;
		linelength++;
		if (c == '\n') {
			f->hlen[j] = linelength;
			linelength = 0;
			if (j==9) {
				f->header[i] = 0;
				break;
			}
			j++;
		}
		i++;
	}

	sscanf(f->header,"%15s",&endian);
	if ( !strcmp(endian,"little-endian") ) f->endian = 0;
	else if ( !strcmp(endian,"big-endian") ) f->endian = 1;
	else {
		fprintf(stderr,"error: cannot read endian in %s\n",filename);
		exit(1);
	}

	j=0;
	i=0;
	space=0;
	while (1) {
		if (j == 3) break;
		c = f->header[i];
		if ((c!=' ') && (c!='\t') && (c!='\n')) {
			if (space) {
				j++;
				space = 0;
			}
		}
		else {
			space = 1;
		};
		i++;
	}
	i--;

	if ( sscanf(f->header+i,"%d %d",&f->rows,&f->cols) != 2) {
		fprintf(stderr,"error in header reading rows, cols in %s\n",filename);
		exit(1);
	}

	sscanf(f->header+f->hlen[0],"%10s %lf %lf %10s %lf %lf %10s",&f->xlabel,&f->xmin,&f->xmax,&f->ylabel,&f->ymin,&f->ymax,&f->zlabel);
	j=0; 
	for (i=0; i<10; i++)  
		j += f->hlen[i];
	if(fseek(fp, (long)j, SEEK_SET)!=0) {
		fprintf(stderr,"error seeking in file %s to pos %d\n",filename,j);
		exit(1);
	}

	printf("malloc f->data: %d",f->rows*f->cols*(size_t)sizeof(double));
	f->data = (double *)malloc(f->rows*f->cols*(size_t)sizeof(double));
	if (f->data == NULL) { perror("f->data"); exit(1); }
	printf(".\n");

	// read data, fill lines bottom up
	if (f->endian == mendian() ) {
		for (i=0; i<f->rows; i++) for (j=0; j<f->cols; j++) {
			if (fread(&f->data[(f->rows-i-1)*f->cols+j], sizeof(double), 1, fp) != 1) {
				fprintf(stderr,"error reading binary stream at pos %i from %s\n",i,filename);
				exit(1);
			}
		}
	} else {
		for (i=0; i<f->rows; i++) for (j=0; j<f->cols; j++) {
			if (fread(&x, sizeof(double), 1, fp) != 1) {
				fprintf(stderr,"error reading binary stream at pos %i from %s\n",i,filename);
				exit(1);
			}
			dswap(&x);
			f->data[(f->rows-i-1)*f->cols+j] = x;
		}
	}

	fclose(fp);
}
