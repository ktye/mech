#include <iostream>
#include <complex>
#include <vector>
#include <gmpxx.h>
using namespace std;
typedef complex<mpq_class> crat;

#define NMAX	21
int main(int args, char **argv) {
	mpq_class R,k,omega_real,omega_imag,abs;
	crat a,b,omega;
	vector<crat> C(NMAX);
	int N;

	if (args!=4) return 1;
	R = argv[1];
	k = argv[2];
	N = atoi(argv[3]);
	if (N>NMAX)	N=NMAX;
	if (N<3)	N=3;

	C[1] = (1,0);
	for (omega_imag = -3; omega_imag <= 3; omega_imag += 0.05 ) {
		for (omega_real = 0; omega_real <= 8; omega_real += 0.05 ) {
			a = crat(-R*omega_imag-k*k,R*omega_real-R*k);
			b = crat(0,k*R);
			C[3] = -a/(mpq_class)(8);
			for (int i=3; i<N-2; i++) {
				C[i+2] = -(a*C[i] + b*C[i-2])/(mpq_class)(i*i+4*i+3);
			}
			C[0] = (0,0);
			for (int i=1; i<N; i+=2) {
				C[0] += C[i];		
			}	
			abs = real(C[0])*real(C[0]) + imag(C[0])*imag(C[0]);
//			cout << abs.get_d() << "\t" << a << "\t" << b << endl;
			cout << abs.get_d() << endl;

		}
	}
	
	return 0;
}
