/* $Id: endian.c,v 1.1 2004/09/08 14:48:23 elmar Exp $ */

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

int mendian(void) {
        int x = 1;
        if(*(char *)&x == 1) return 0; // little-endian
        else return 1; // big-endian
}       
