/*
 * cpoly - roots of complex polynoms
 */

#include <iostream>
#include <complex>
#include "cpoly.h"
using namespace std;
#define MAX 8
int main(int args, char **argv) {
	//int cpoly( const double *opr, const double *opi, int degree, double *zeror, double *zeroi )	
	double coeff_real[MAX];	
	double coeff_imag[MAX];	
	double zeros_real[MAX-1];	
	double zeros_imag[MAX-1];	
	double real, imag;
	double k,R;
	int i,n,ret;

	if (args!=3) return 1;
	R=atof(argv[1]);
	k=atof(argv[2]);

	/* r=3 -> linear */
	cout << ":0:17:brown" << endl;
	cout << k << "\t" << -(k*k+8.0)/R << endl;
	cout << endl;

	/* r=5 -> omega quadratisch */
	/* coeff_[0]: coeff vor hoechster Potenz! */
	coeff_real[0] = 1.0;
	coeff_imag[0] = 0.0;
	coeff_real[1] = -2.0*k;
	coeff_imag[1] = (2.0*k*k+24.0)/R;
	coeff_real[2] = k*k - (pow(k,4)+24.0*k*k+192.0)/(R*R);
	coeff_imag[2] = -(2.0*k*k*k + 16.0*k)/R;
	ret = cpoly(coeff_real, coeff_imag, 2, zeros_real, zeros_imag);
	cout << ":0:7:red" << endl;
	cout << zeros_real[0] << "\t" << zeros_imag[0] << endl;
	cout << zeros_real[1] << "\t" << zeros_imag[1] << endl;
	cout << endl;

	/* r=7 -> omega kubisch */
	coeff_real[0] = 1.0;
	coeff_imag[0] = 0.0;
	coeff_real[1] = -3.0*k;
	coeff_imag[1] = (3.0*k*k + 48.0)/R;
	coeff_real[2] = 3.0*k*k-(3.0*pow(k,4)+96.0*k*k+1152.0)/(R*R);
	coeff_imag[2] = -(6.0*pow(k,3)+64.0*k)/R;
	coeff_real[3] = -pow(k,3) + (3.0*pow(k,5)+64.0*pow(k,3)+768*k)/(R*R);
	coeff_imag[3] = (3.0*pow(k,4)+16.0*k*k)/R-(pow(k,6)+48.0*pow(k,4)+1152*k*k+9216.0)/pow(R,3);
	ret = cpoly(coeff_real, coeff_imag, 3, zeros_real, zeros_imag);
	cout << ":0:7:green" << endl;
	cout << zeros_real[0] << "\t" << zeros_imag[0] << endl;
	cout << zeros_real[1] << "\t" << zeros_imag[1] << endl;
	cout << zeros_real[2] << "\t" << zeros_imag[2] << endl;
	cout << endl;

	/* r=9 -> omega quartisch */
	coeff_real[0] = 1.0;
	coeff_imag[0] = 0.0;
	coeff_real[1] = -4.0*k;
	coeff_imag[1] = (4.0*k*k+80.0)/R;
	coeff_real[2] = 6.0*k*k-(6.0*pow(k,4)+240.0*k*k+3840.0)/(R*R);
	coeff_imag[2] = -(12.0*k*k*k+160.0*k)/R;
	coeff_real[3] = -4.0*k*k*k+(12.0*pow(k,5)+320.0*pow(k,3)+5120.0*k)/(R*R);
	coeff_imag[3] = (12.0*pow(k,4)+80.0*k*k)/R-(4.0*pow(k,6)+240.0*pow(k,4)+7680.0*k*k+92160.0)/(R*R*R);
	coeff_real[4] = pow(k,4)+(pow(k,8)+80.0*pow(k,6)+3840.0*pow(k,4)+92160.0*k*k+737280.0)/pow(R,4)-(6.0*pow(k,6)+80.0*pow(k,4)+1664.0*k*k)/(R*R);
	coeff_imag[4] = -4.0*pow(k,5)/R+(4.0*pow(k,7)+160.0*pow(k,5)+5120.0*pow(k,3)+61440.0*k)/pow(R,3);
	ret = cpoly(coeff_real, coeff_imag, 4, zeros_real, zeros_imag);
	cout << ":0:7:yellow" << endl;
	cout << zeros_real[0] << "\t" << zeros_imag[0] << endl;
	cout << zeros_real[1] << "\t" << zeros_imag[1] << endl;
	cout << zeros_real[2] << "\t" << zeros_imag[2] << endl;
	cout << zeros_real[3] << "\t" << zeros_imag[3] << endl;
	cout << endl;

	/* r=11 -> omega quintisch */
	coeff_real[0] = 1.0;
	coeff_imag[0] = 0.0;
	coeff_real[1] = -5.0*k;
	coeff_imag[1] = (5.0*k*k+120.0)/R;
	coeff_real[2] = 10.0*k*k-(10.0*pow(k,4)+480*k*k+9600.0)/(R*R);
	coeff_imag[2] = -(20.0*k*k*k+320.0*k)/R;
	coeff_real[3] = -10.0*k*k*k + (30.0*pow(k,5)+960.0*k*k*k+19200.0*k)/(R*R);
	coeff_imag[3] = (30.0*pow(k,4)+240.0*k*k)/R-(10.0*pow(k,6)+720.0*pow(k,4)+28800.0*k*k+460800.0)/(R*R*R);
	coeff_real[4] = 5.0*pow(k,4)-(30.0*pow(k,6)+480.0*pow(k,4)+12544.0*k*k)/(R*R) + (5.0*pow(k,8)+480.0*pow(k,6)+28800.0*pow(k,4)+921600.0*k*k+11059200.0)/pow(R,4);
	coeff_imag[4] = -20.0*pow(k,5)/R+(20.0*pow(k,7)+960.0*pow(k,5)+38400.0*pow(k,3)+614400.0*k)/pow(R,3);
	coeff_real[5] = -pow(k,5)+(10.0*pow(k,7)+2944.0*pow(k,3))/(R*R)-(5.0*pow(k,9)+320.0*pow(k,7)+19200.0*pow(k,5)+614400.0*pow(k,3)+73772800.0*k)/pow(R,4);
	coeff_imag[5] = (5.0*pow(k,6)-40.0*pow(k,4))/R -(10.0*pow(k,8)+240.0*pow(k,6)+12544.0*pow(k,4)+199680.0*k*k)/pow(R,3) + (pow(k,10)+120.0*pow(k,8)+9600.0*pow(k,6)+460800.0*pow(k,4)+11059200.0*k*k+88473600.0)/pow(R,5);
	ret = cpoly(coeff_real, coeff_imag, 5, zeros_real, zeros_imag);
	cout << ":0:7:pink" << endl;
	cout << zeros_real[0] << "\t" << zeros_imag[0] << endl;
	cout << zeros_real[1] << "\t" << zeros_imag[1] << endl;
	cout << zeros_real[2] << "\t" << zeros_imag[2] << endl;
	cout << zeros_real[3] << "\t" << zeros_imag[3] << endl;
	cout << zeros_real[4] << "\t" << zeros_imag[4] << endl;
	cout << endl;

	/* r=13 -> omega sixtisch */
	coeff_real[0] = 1.0;
	coeff_imag[0] = 0.0;
	coeff_real[1] = -6.0*k;
	coeff_imag[1] = (6.0*k*k+168.0)/R;
	coeff_real[2] = 15.0*k*k-(15.0*pow(k,4)+840.0*k*k+20160.0)/(R*R);
	coeff_imag[2] = -(30.0*k*k*k-560.0*k)/R;
	coeff_real[3] = -20.0*k*k*k+(60.0*pow(k,5)+2240.0*k*k*k+53760.0*k)/(R*R);
	coeff_imag[3] = (60.0*pow(k,4)+560.0*k*k)/R - (20.0*pow(k,6)+1680.0*pow(k,4)+80640.0*k*k+1612800.0)/(R*R*R);
	coeff_real[4] = 15.0*pow(k,4)-(90.0*pow(k,6)+1680.0*pow(k,4)+52864.0*k*k)/(R*R)+(15.0*pow(k,8)+1680.0*pow(k,6)+120960.0*pow(k,4)+4838400.0*k*k+77414400.0)/pow(R,4);
	coeff_imag[4] = -60.0*pow(k,5)+(60.0*pow(k,7)+3360.0*pow(k,5)+161280.0*k*k*k+3225600.0*k)/pow(R,3);
	coeff_real[5] = -6.0*pow(k,5)+(60.0*pow(k,7)+25088.0*pow(k,3))/(R*R)-(30.0*pow(k,9)+2240.0*pow(k,7)+161280.0*pow(k,5)+6451200.0*k*k*k+103219200*k)/pow(R,4);
	coeff_imag[5] = (30.0*pow(k,6)-280.0*pow(k,4))/R-(60.0*pow(k,8)+1680*pow(k,6)+105728.0*pow(k,4)+2107392.0*pow(k,2))/(R*R*R)+(6.0*pow(k,10)+840.0*pow(k,8)+80640.0*pow(k,6)+4838400.0*pow(k,4)+154828800.0*k*k+1857945600.0)/pow(R,5);
	coeff_real[6] = pow(k,6)-(15.0*pow(k,8)+280.0*pow(k,6)-5824.0*pow(k,4))/(R*R)+(15.0*pow(k,10)+560.0*pow(k,8)+52864.0*pow(k,6)+210707392.0*pow(k,4)+33546240.0*k*k)/pow(R,4)-(pow(k,12)+168*pow(k,10)+20160.0*pow(k,8)+1612800.0*pow(k,6)+77414400.0*pow(k,4)+1857945600.0*pow(k,2)+14863564800.0)/pow(R,6);
	coeff_imag[6] = (112.0*pow(k,5)-6.0*pow(k,7))/R+(10.0*pow(k,9)+25088.0*pow(k,5)+448512.0*pow(k,3))/pow(R,3)-(6.0*pow(k,11)+560.0*pow(k,9)+53760.0*pow(k,7)+3225600.0*pow(k,5)+103219200.0*pow(k,3)+1238630400.0*k)/pow(R,5);
	ret = cpoly(coeff_real, coeff_imag, 6, zeros_real, zeros_imag);
	cout << ":0:7:white" << endl;
	cout << zeros_real[0] << "\t" << zeros_imag[0] << endl;
	cout << zeros_real[1] << "\t" << zeros_imag[1] << endl;
	cout << zeros_real[2] << "\t" << zeros_imag[2] << endl;
	cout << zeros_real[3] << "\t" << zeros_imag[3] << endl;
	cout << zeros_real[4] << "\t" << zeros_imag[4] << endl;
	cout << zeros_real[5] << "\t" << zeros_imag[5] << endl;
	cout << endl;


	return 0;
}
