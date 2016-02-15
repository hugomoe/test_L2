#include "fft_zoom.h"
#include "parameters.h"
#include "aux_fun.h"



//A RENDRE PROPRE


//extrapolate by a constant value
float getsample_cons(float *x, int w, int h, int pd, int i, int j, int l)
{
	static float value = 0;
	if (w == 0 && h == 0)
		value = *x;
	if (i < 0 || i >= w || j < 0 || j >= h)
		return 0;
	if (l >= pd)
		l = pd - 1;
	return x[(i+j*w)*pd + l];
}

// auxiliary function for bilinear interpolation
static float evaluate_bilinear_cell(float a, float b, float c, float d,
							float x, float y)
{
	return a * (1-x) * (1-y)
	     + b * ( x ) * (1-y)
	     + c * (1-x) * ( y )
	     + d * ( x ) * ( y );
}

// bilinear interpolation
float bilinear_interpolation_at(float *x, int w, int h, int pd,
		float p, float q, int l)
{
	int ip = floor(p);
	int iq = floor(q);
	float a = getsample_cons(x, w, h, pd, ip  , iq  , l);
	float b = getsample_cons(x, w, h, pd, ip+1, iq  , l);
	float c = getsample_cons(x, w, h, pd, ip  , iq+1, l);
	float d = getsample_cons(x, w, h, pd, ip+1, iq+1, l);
	return evaluate_bilinear_cell(a, b, c, d, p-ip, q-iq);
}






//"VRAI" code

int apply_homo_ground_truth(float *img,float *img_f,int w,int h,int w_f,int h_f,double H[3][3]){	
	
	double fzoom=ZOOM;
	double HH[3][3];
	HH[0][0]=H[0][0];
	HH[0][1]=H[0][1];
	HH[0][2]=H[0][2]*fzoom;
	HH[1][0]=H[1][0];
	HH[1][1]=H[1][1];
	HH[1][2]=H[1][2]*fzoom;
	HH[2][0]=H[2][0]/fzoom;
	HH[2][1]=H[2][1]/fzoom;
	HH[2][2]=H[2][2];
	
	

	float *img_zoom = malloc(3*sizeof(float)*w*h*ZOOM*ZOOM);
	float *img_f_zoom = malloc(3*sizeof(float)*w_f*h_f*ZOOM*ZOOM);
	
	//bilinear zoom-in
	/*
	for(int l=0;l<3;l++){
	for(int i=0;i<w;i++){
		for(int j=0;j<h;j++){
			for(int u=0;u<ZOOM;u++){
				for(int v=0;v<ZOOM;v++){
					float x=u/fzoom;
					float y=v/fzoom;
					int id,jd;
					if(i==w-1){id=0;}else{id=i+1;}
					if(j==h-1){jd=0;}else{jd=j+1;}
					img_zoom[3*(i*ZOOM+u+(j*ZOOM+v)*w*ZOOM)+l]=(1-x)*(1-y)*(img[(i+j*w)*3+l])+(1-x)*y*(img[(i+jd*w)*3+l])
					+x*(1-y)*(img[(id+j*w)*3+l])+x*y*(img[(id+jd*w)*3+l]);
				}
			}
		}
	}
	}
	*/
	
	//symmetrize the input
	float *img_sym = malloc(2*w*2*h*3*sizeof(float));
	int i_sym, j_sym;
	for(int l=0;l<3;l++){
		for(int i=0;i<2*w;i++){
			i_sym = i-w/2;
			if(i_sym<0){i_sym = -1-i_sym;}
			else if(i_sym>w-1){i_sym = 2*w-1-i_sym;}
			for(int j=0;j<2*h;j++){
				j_sym = j-h/2;
				if(j_sym<0){j_sym = -1-j_sym;}
				else if(j_sym>w-1){j_sym = 2*h-1-j_sym;}
				img_sym[3*(i+2*w*j)+l] = img[3*(i_sym+w*j_sym)+l];
			}
		}
	}
	//zoom-in by zero-padding
	float *img_sym_zoom = malloc(2*w*ZOOM*2*h*ZOOM*3*sizeof(float));
	zoom(img_sym,2*w,2*h,3,2*w*ZOOM,2*h*ZOOM,img_sym_zoom);
	//erase the symmetrization
	for(int l=0;l<3;l++){
		for(int i=0;i<w*ZOOM;i++){
			for(int j=0;j<h*ZOOM;j++){
				img_zoom[3*(i + w*ZOOM * j) + l] = img_sym_zoom[3*(i+w*ZOOM/2 + 2*w*ZOOM * (j+h*ZOOM/2)) + l];
			}
		}
	}

	
	
	//naive warping
	for(int l=0;l<3;l++)
	for (int j = 0; j < ZOOM*h_f; j++)
	for (int i = 0; i < ZOOM*w_f; i++)
	{
		double p[2] ={i,j};
		
		apply_homography_1pt(p, HH, p);
		int idx = 3*(ZOOM * w_f * j + i)+l;
		img_f_zoom[idx] = bilinear_interpolation_at(img_zoom, ZOOM*w, ZOOM*h, 3, p[0], p[1], l);//POURQUOI ?
	}
	


	//zoom-out by convolution with a gaussian kernel
	double sigma = 0.5 * ZOOM; //double sigma = 0.6 * ZOOM; //results look better with 0.5
	int taps = ceil(6*sigma);
	double *gauss = malloc((2*taps+1)*sizeof(double));
	double tot = 0;
	for(int k=-taps;k<=taps;k++){
			tot += (gauss[k+taps] = exp(-(pow(k,2))/(2*pow(sigma,2))));
	}
	for(int k=-taps;k<=taps;k++){gauss[k+taps] /= tot;}


	float *img_aux = malloc(3*sizeof(float)*w_f*h_f*ZOOM); //size w_f x h_f*ZOOM
	if(img_aux==NULL){
        printf("img_aux is too large");
        exit(1);
    	}
	int idx, i0, j0;
	for(int l=0;l<3;l++){
		//horizontal convolution
		for (int i = 0; i < w_f; i++){
			for (int j = 0; j < h_f*ZOOM; j++){
				float v = 0;
				for(int k=-taps;k<=taps;k++){
					i0 = good_modulus(i*ZOOM+k,w_f*ZOOM);
					v += gauss[k+taps]*img_f_zoom[(i0+j*ZOOM*w_f)*3+l];
				}
				int idx = l+3*(w_f*j+i);
				img_aux[idx] = v;
			}
		}
		//vertical convolution
		for (int i = 0; i < w_f; i++){
			for (int j = 0; j < h_f; j++){
				float v = 0;
				for(int k=-taps;k<=taps;k++){
					j0 = good_modulus(j*ZOOM+k,h_f*ZOOM);
					v += gauss[k+taps]*img_aux[(i+j0*w_f)*3+l];
				}
				int idx = l+3*(w_f*j+i);
				img_f[idx] = v;
			}
		}
	}

	free(img_aux);
	free(img_sym);
	free(img_sym_zoom);
	free(img_zoom);
	free(img_f_zoom);
	free(gauss);
	return 0;
}
