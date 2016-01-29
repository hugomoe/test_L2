#ifndef AUX_FUN
#define AUX_FUN

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "parameters.h"

//a relaxed equality on doubles
bool eq(double a,double b){if(a<b+pow(2,-PREC_EQ) && a>b-pow(2,-PREC_EQ)){return true;}{return false;}}

//transport a point by an homography
void apply_homography_1pt(double y[2], double H[3][3], double x[2])
{
	double X = H[0][0] * x[0] + H[0][1] * x[1] + H[0][2];
	double Y = H[1][0] * x[0] + H[1][1] * x[1] + H[1][2];
	double Z = H[2][0] * x[0] + H[2][1] * x[1] + H[2][2];
	y[0] = X / Z;
	y[1] = Y / Z;
}

//a good definition of n%p, even for huge and negative numbers
int good_modulus(int nn, int p)
{
	if (!p) return 0;
	if (p < 1) return good_modulus(nn, -p);

	unsigned int r;
	if (nn >= 0)
		r = nn % p;
	else {
		unsigned int n = nn;
		r = p - (-n) % p;
		if (r == p)
			r = 0;
	}
	return r;
}

#endif //AUX_FUN
