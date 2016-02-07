/**
  * compare a sub-image of img_dec.png, img_grd.png and img_rip.png (created by viho_alt)
  * viho_alt must have been runned before running compare_sub
  * 
  * compilation
  * 	gcc-5 -O3 compare_sub.c -ltiff -ljpeg -lpng -o compare_sub
  * usage
  * 	./compare_sub
  *		or
  * 	./compare_sub i1 j1 i2 j2
  */



#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "iio.c"
#include "parameters.h"



int main(int argc,char *argv[]){
	if (argc != 1 && argc != 5) {
		printf("usage :\n\tno args\nor\n\ti1 j1 i2 j2"); 
		return 1;
	}



	int i1,j1,i2,j2; //coordinates of two opposite vertices of a sub-image
	switch(argc){
		case 5:
		i1 = strtol(argv[1],NULL,10);
		j1 = strtol(argv[2],NULL,10);
		i2 = strtol(argv[3],NULL,10);
		j2 = strtol(argv[4],NULL,10);
		break;
		case 1:
		default:
		i1 = 0;
		j1 = 0;
		i2 = WOUT-1;
		j2 = HOUT-1;
		break;
	}



	int w, h, pd,
		w_temp, h_temp, pd_temp;

	float *img_dec;
	img_dec = iio_read_image_float_vec("img_dec.png", &w, &h, &pd);

	float *img_grd;
	img_grd = iio_read_image_float_vec("img_grd.png", &w_temp, &h_temp, &pd_temp);
	if(w != w_temp || h != h_temp || pd != pd_temp){
		printf("error : img_dec.png and img_grd.png do not have the same size");
		return 1;
	}

	float *img_rip;
	img_rip = iio_read_image_float_vec("img_rip.png", &w_temp, &h_temp, &pd_temp);
	if(w != w_temp || h != h_temp || pd != pd_temp){
		printf("error : img_rip.png does not have the same size than img_dec and img_grd");
		return 1;
	}



	int i_tl, j_tl, //coordinates of the top left pixel of the sub-image
		i_br, j_br; //coordinates of the bottom right pixel of the sub-image
	if(i1<i2){
		i_tl = i1;
		i_br = i2;
	}else{
		i_tl = i2;
		i_br = i1;
	}
	if(j1<j2){
		j_tl = j1;
		j_br = j2;
	}else{
		j_tl = j2;
		j_br = j1;
	}



	//the sub-image must be a subset of the output image
	if(i_tl<0){i_tl=0;}
	if(j_tl<0){j_tl=0;}
	if(i_br>w-1){i_br=w-1;}
	if(j_br>h-1){j_br=h-1;}



	float l1_dg = 0., l2_dg = 0., //compare decomposition and ground truth
		l1_dr = 0., l2_dr = 0., //compare decomposition and ripmap
		l1_rg = 0., l2_rg = 0.; //compare ripmap and ground truth
	float diff;
	int idx;

	for(int l=0 ; l<3 ; l++)
		for(int i=i_tl ; i<=i_br ; i++)
			for(int j=j_tl ; j<=j_br ; j++){
				idx = 3*(i+w*j) + l;

				diff = fabs(img_dec[idx]-img_grd[idx]);
				l1_dg += diff;
				l2_dg += pow(diff,2);

				diff = fabs(img_dec[idx]-img_rip[idx]);
				l1_dr += diff;
				l2_dr += pow(diff,2);

				diff = fabs(img_rip[idx]-img_grd[idx]);
				l1_rg += diff;
				l2_rg += pow(diff,2);
	}

//normalize
	int N = (i_br-i_tl+1)*(j_br-j_tl+1)*3; //size of the sub-image
	l1_dg = l1_dg/(float) N;
	l2_dg = sqrt(l2_dg/(float) N);
	l1_dr = l1_dr/(float) N;
	l2_dr = sqrt(l2_dr/(float) N);
	l1_rg = l1_rg/(float) N;
	l2_rg = sqrt(l2_rg/(float) N);

//print
	printf("l1-error :\n\t\tGround Truth\tRipMap\nDecomposition\t%f\t%f\nRipMap\t\t%f\n\n",l1_dg,l1_dr,l1_rg);
	printf("l2-error :\n\t\tGround Truth\tRipMap\nDecomposition\t%f\t%f\nRipMap\t\t%f\n\n",l2_dg,l2_dr,l2_rg);
	return 0;
}
