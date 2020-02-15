/* http://research.microsoft.com/~hollasch/cgindex/math/cubic.c */
/* This code snippet is used to solve for the roots of a cubic polynomial. */
#define pi 3.1415926535897
#include <math.h>
solve3(coeff, root)
	double          coeff[4];
	double          root[3];
{
	double          a = coeff[3], b = coeff[2], c = coeff[1], d = coeff[0];
	int             n_roots;
	double          p, q, disc, b_over_3a, c_over_a, d_over_a;
	if (a == 0)
		return solve2(coeff, root);
	b_over_3a = b / (3 * a);
	c_over_a = c / a;
	d_over_a = d / a;
	p = b_over_3a * b_over_3a;
	q = 2 * b_over_3a * p - b_over_3a * c_over_a + d_over_a;
	p = c_over_a / 3 - p;
	disc = q * q + 4 * p * p * p;
	if (disc < 0) {
		double          r, theta, temp;
		r = .5 * sqrt(-disc + q * q);
		theta = atan2(sqrt(-disc), -q);
		temp = 2 * cbrt(r);
		root[0] = temp * cos(theta / 3);
		root[1] = temp * cos((theta + pi + pi) / 3);
		root[2] = temp * cos((theta - pi - pi) / 3);
		n_roots = 3;
	} else {
		double          alpha, beta;
		alpha = .5 * (sqrt(disc) - q);
		beta = -q - alpha;
		root[0] = cbrt(alpha) + cbrt(beta);
		if (disc > 0)
			n_roots = 1;
		else {
			root[1] = root[2] = -.5 * root[0];
			n_roots = 3;
		}
	}
	{
		int             i;
		for (i = 0; i < n_roots; i++)
			root[i] -= b_over_3a;
	}
	return n_roots;
}
solve2(coeff, root)
	double          coeff[3];
	double          root[2];
{
	double          a = coeff[2], b = coeff[1], c = coeff[0];
	double          disc, b_over_2a, c_over_a;
	if (a == 0)
		return solve1(coeff, root);
	b_over_2a = b / (2 * a);
	c_over_a = c / a;
	disc = b_over_2a * b_over_2a - c_over_a;
	if (disc < 0)
		return 0;
	else if (disc == 0) {
		root[0] = -b_over_2a;
		return 1;
	} else {
		root[0] = -b_over_2a + sqrt(disc);
		root[1] = -2 * b_over_2a - root[0];
		return 2;
	}
}
solve1(coeff, root)
	double          coeff[2];
	double          root[1];
{
	double          a = coeff[1], b = coeff[0];
	if (a == 0)
		return 0;
	root[0] = -b / a;
	return 1;
}
