#define ZOOM 5
#include "fft_zoom.h"
#include "parameters.h"
#include "aux_fun.h"

//la fft fait beaucoup de ringing



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
	
	
	float *img_aux = malloc(3*sizeof(float)*w*h*ZOOM*ZOOM);
	float *img_aux2 = malloc(3*sizeof(float)*w_f*h_f*ZOOM*ZOOM);
	for(int i=0;i<w*h*ZOOM*ZOOM*3;i++){img_aux[i]=0;}
	
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
					img_aux[3*(i*ZOOM+u+(j*ZOOM+v)*w*ZOOM)+l]=(1-x)*(1-y)*(img[(i+j*w)*3+l])+(1-x)*y*(img[(i+jd*w)*3+l])
					+x*(1-y)*(img[(id+j*w)*3+l])+x*y*(img[(id+jd*w)*3+l]);
				}
			}
		}
	}
	}
	*/
	
	//zoom-in by zero-padding
	int w_zoom = ZOOM*w; int h_zoom = ZOOM*h;
	zoom(img,w,h,3,w_zoom,h_zoom,img_aux);
	
	
	
	//naive warping
	for(int l=0;l<3;l++)
	for (int j = 0; j < ZOOM*h_f; j++)
	for (int i = 0; i < ZOOM*w_f; i++)
	{
		double p[2] ={i,j};
		
		apply_homography_1pt(p, HH, p);
		int idx = 3*(ZOOM * w_f * j + i)+l;
		img_aux2[idx] = bilinear_interpolation_at(img_aux, ZOOM*w, ZOOM*h, 3, p[0], p[1], l);//POURQUOI ?
	}
	


	//zoom-out by convolution with a gaussian kernel
	int taps = 3*ZOOM;
	double sigma = 0.6 * ZOOM;
	double *gauss = malloc(pow((2*taps+1),2)*sizeof(float));
	double tot = 0;
	for(int i=-taps;i<=taps;i++){
		for(int j=-taps;j<=taps;j++){
			tot += (gauss[i+taps+(j+taps)*(2*taps+1)] = exp(-(pow(i,2)+pow(j,2))/(2*pow(sigma,2))));
		}
	}	
	for(int u=0;u<pow(2*taps+1,2);u++){gauss[u]=gauss[u]/tot;}

	
	for(int l=0;l<3;l++){
		for (int j = 0; j < h_f; j++){
			for (int i = 0; i < w_f; i++){
				float v=0;
				for(int i2=-taps;i2<=taps;i2++){
					for(int j2=-taps;j2<=taps;j2++){
						int i1 = good_modulus(i*ZOOM+i2,w_f*ZOOM);
						int j1 = good_modulus(j*ZOOM+j2,h_f*ZOOM);
						v += gauss[i2+taps+(j2+taps)*(2*taps+1)]*img_aux2[(i1+j1*ZOOM*w_f)*3+l];
						int idx = l + 3 * (w_f * j + i);
						img_f[idx] = v;
					}
				}
			}
		}
	}
	
	
	return 0;
}
