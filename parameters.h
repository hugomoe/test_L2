#ifndef PARAMETERS
#define PARAMETERS

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

//size of the output image
#define WOUT 512
#define HOUT 512

//parameters for raised cosine-weighted sinc
#define SUPP 4 //the support of the filter is included in [-SUPP,SUPP]
#define PERIOD 1. //period of the raised cosine-weighted sinc
#define BETA 0.36 //roll-off factor of the raised cosine-weighted sinc
#define PREC 10 //the precomputation of the values of the filter will be at precision 2^(-PREC)

//parameters for the gaussian kernel in ripmap.h
#define TAPSR 5 //number of non zero values of the kernel (must be odd)
#define SIG 0.6 //variance

//parameters for the distance D in ripmap.h
#define D_BIAS 0. //added to the distance
#define D_COEFF 1.

//zoom for the ground truth
#define ZOOM 5

//precision of the computation
#define PREC_EQ 20 //an equality will be at precision 2^(-PREC_EQ)
#define PI 3.14159265358979323

#endif //PARAMETERS
