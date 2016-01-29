// On fait un zoom in par zero padding.


#include "iio.h"
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <fftw3.h>


 int zoom (float * img, int w,int h,int pd,int kw,int kh,float *img_final)
 { void fourierForward(float* in,float* reOut,float* imOut,unsigned int largeur,unsigned int hauteur) ;
   void fourierBackward(float* reIn,float* imIn,float* out,unsigned int largeur,unsigned int hauteur) ;

 /* Soit img le pointeur de l'image de départ */
 /* soit w la longueur et h la largeur de l'image de départ */
 // kw et kh seront respectivement la longueur et la largeur de l'image finale.


	float *refftimg= malloc(w*h*sizeof(float));
	float *imfftimg= malloc(w*h*sizeof(float));
	float *img_aux= malloc(w*h*sizeof(float));
	float *refftimg2=malloc(kw*kh*sizeof(float));
	float *imfftimg2=malloc(kw*kh*sizeof(float));
	float *img_final_aux=malloc(sizeof(float)*kw*kh);
	
for(int l=0;l<pd;l++){ 

	//met a zero
	for(int i=0;i<w*h;i++){img_aux[i]=img[pd*i+l];imfftimg[i]=0;refftimg[i]=0;}
	for(int i=0;i<kw*kh;i++){refftimg2[i]=0;imfftimg2[i]=0;img_final_aux[i]=0;}
   
	fourierForward(img_aux,refftimg,imfftimg,w,h) ; // transformé de Fourier directe

    // début zéro-padding
	for (int i=0;i<kw;i++){
		for (int j=0;j<kh;j++){  
			refftimg2[i+j*kw]=0 ; 
            imfftimg2[i+j*kw]=0 ;     
        }
    }

    for (int i=0;i<w;i++){
       for (int j=0;j<h;j++){  
       		refftimg2[kw/2-w/2+i+kw*(kh/2-h/2+j)]=refftimg[w*j+i] ;     
        	imfftimg2[kw/2-w/2+i+kw*(kh/2-h/2+j)]=imfftimg[w*j+i] ; 
    	}
    }
    // fin zéro-padding
	fourierBackward(refftimg2, imfftimg2, img_final_aux, kw,kh); // transformée de Fourier inverse
	

	for(int i=0;i<kw;i++){
		for(int j=0;j<kh;j++){
			img_final[(i+j*kw)*pd+l]= (kw/(float)w)*(kh/(float)h)*img_final_aux[i+j*kw];
		}
	}

}

	free(refftimg2);
	free(imfftimg2);
	free(refftimg);
	free(imfftimg);
	free(img_aux);
	free(img_final_aux);
	
return 0 ;
}




// transformées de Fourier :

void fourierForward(float* in,
                    float* reOut,
                    float* imOut,
                    unsigned int largeur,
                    unsigned int hauteur)
{
   fftw_complex* spatial_repr;
   fftw_complex* frequency_repr;
   unsigned int i;
   unsigned int j;
   fftw_plan plan;
   int x,y;
 
   spatial_repr= fftw_malloc(sizeof(fftw_complex)*largeur*hauteur);
   frequency_repr= fftw_malloc(sizeof(fftw_complex)*largeur*hauteur);
 
 
   for(i=0;i<largeur*hauteur;i++)
   {
      spatial_repr[i][0] = in[i];
      spatial_repr[i][1] =  0.0f;
   }
 
   /*on calcule le plan d'exécution*/
   plan=fftw_plan_dft_2d(hauteur, largeur, spatial_repr, frequency_repr, FFTW_FORWARD, FFTW_ESTIMATE);
 
   /*on calcule la transformée*/
   fftw_execute(plan);
 
  for(j=0;j<hauteur;j++)
      for(i=0;i<largeur;i++)
      {
	        /*on recentre l'image*/
	      x=good_modulus(i+largeur/2,largeur);
	      y=good_modulus(j+hauteur/2,hauteur);
          reOut[y*largeur+x]=frequency_repr[j*largeur+i][0];
          imOut[y*largeur+x]=frequency_repr[j*largeur+i][1];
      }
 
   fftw_destroy_plan(plan);
   fftw_free(spatial_repr);
   fftw_free(frequency_repr);
 
}
		
 




/* reIn : partie réel de l'image dans l'espace de Fourier
 * imIn : partie imaginaire de l'image dans l'espace de Fourier
 * out : image de sortie
 * largeur : largeur des images d'entrée et de sortie
 */
void fourierBackward(float* reIn,
                     float* imIn,
                     float* out,
                     unsigned int largeur,
                     unsigned int hauteur)
{
   fftw_complex* spatial_repr;
   fftw_complex* frequency_repr;
   unsigned int i;
   unsigned int j;
   int x,y;
   fftw_plan plan;
 
   spatial_repr= fftw_malloc(sizeof(fftw_complex)*largeur*hauteur);
   frequency_repr= fftw_malloc(sizeof(fftw_complex)*largeur*hauteur);
 
   for(j=0;j<hauteur;j++)
      for(i=0;i<largeur;i++)
      {
          /*on décentre*/
	      x=good_modulus(i+largeur/2,largeur);
	      y=good_modulus(j+hauteur/2,hauteur);

	      frequency_repr[j*largeur+i][0]=reIn[y*largeur+x];
	      frequency_repr[j*largeur+i][1]=imIn[y*largeur+x];
      }
 
  plan=fftw_plan_dft_2d(hauteur, largeur, frequency_repr, spatial_repr, FFTW_BACKWARD, FFTW_ESTIMATE);
 
  fftw_execute(plan);
 
   /*on retranscrit l'image complexe en image réelle, sans oublier de diviser par largeur*hauteur*/
   for(i=0;i<largeur*hauteur;i++)
   {
      out[i]=spatial_repr[i][0]/(largeur*hauteur);
   }
 
   fftw_destroy_plan(plan);
   fftw_free(spatial_repr);
   fftw_free(frequency_repr);
}
