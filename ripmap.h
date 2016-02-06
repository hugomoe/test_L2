//SECTION 5.5 Ripmap
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "parameters.h"
#include "aux_fun.h"

//construit le ripmap (qui est une image 4 fois plus grande, toujours de profondeur 3)
//On a un ripmap de type float*
//Les deux premiers entiers correspondent a la taille du rectangle
//on a une fonction qui donne les coordonné dans le ripmap en fonction de la taille des rectangles





int coord(int i,int j,int u,int v,int w,int l){
	int x = good_modulus(u,w/pow(2,i));
	int y = good_modulus(v,w/pow(2,j));
	return (2*(1-1/pow(2,i))*w + 4*(1-1/pow(2,j))*w*w + x + y*2*w)*3+l;
}


//On construit le ripmap. On a ici fait un flitre gaussien dans la direction de la compression

//deux fonctions de flitrage verticale et horizontale

float gaussian_filter_h(float *r,int i,int j,int u,int v,int w,int l,float *g){
	float a,b;
	a = b = 0;
	for(int f=0;f<TAPSR;f++){a += r[coord(i-1,j,2*u+f-(TAPSR-1)/2,v,w,l)]*g[f];}
	for(int f=0;f<TAPSR;f++){b += r[coord(i-1,j,2*u+1+f-(TAPSR-1)/2,v,w,l)]*g[f];}
	return (a+b)/2;
}

float gaussian_filter_v(float *r,int i,int j,int u,int v,int w,int l,float *g){
	float a,b;
	a = b = 0;
	for(int f=0;f<TAPSR;f++){a += r[coord(i,j-1,u,2*v+f-(TAPSR-1)/2,w,l)]*g[f];}
	for(int f=0;f<TAPSR;f++){b += r[coord(i,j-1,u,2*v+1+f-(TAPSR-1)/2,w,l)]*g[f];}
	return (a+b)/2;
}


void build_ripmap(float *img,float *r,int w,int h,int pd,int logw){
	int l,ll,i,j,u,v,w1,w2;
	
	//on déclare un filtre gaussien 1D;
	float gauss1D[TAPSR];
	for(int f=0;f<TAPSR;f++){
		gauss1D[f]=exp(-pow((TAPSR-1)/2-f,2)/(2*pow(SIG,2)))/(sqrt(2*PI)*SIG);
	}
	float total = 0;
	for(int f=0;f<TAPSR;f++){total = total + gauss1D[f];}
	for(int f=0;f<TAPSR;f++){gauss1D[f] = gauss1D[f]/total;}
	
	for(u=0;u<w;u++){
		for(v=0;v<h;v++){
			for(l=0;l<3;l++){
				if(l>pd-1){ll=pd-1;}else{ll=l;}
				r[coord(0,0,u,v,w,l)]=img[(u+v*w)*pd+ll];
			}
		}
	}
//on construit l'image d'origine en (0,0)	

	for(i=1;pow(2,i)<=w;i++){
		w1=w/pow(2,i);
		for(u=0;u<w1;u++){
			for(v=0;v<w;v++){
				for(l=0;l<3;l++){    
					r[coord(i,0,u,v,w,l)] = gaussian_filter_h(r,i,0,u,v,w,l,gauss1D);
				}
			}
		}
	}
//ici on moyenne simplement pour avoir tous les rectangles aplatis en largeur	

	for(i=0;pow(2,i)<=w;i++){
		w1=w/pow(2,i);
		for(j=1;pow(2,j)<=w;j++){
			w2=w/pow(2,j);
			for(u=0;u<w1;u++){
				for(v=0;v<w2;v++){
					for(l=0;l<3;l++){
						r[coord(i,j,u,v,w,l)] = gaussian_filter_v(r,i,j,u,v,w,l,gauss1D);
					}
				}
			}
		}
	}
//ici on obtient les rectangles aplatis en hauteur a partir de ceux calculés avant 
	//return r;
}


//Pour la fonction de distance, on prend dh = |du/dx|+|du/dy| et dl = |dv/dx| + |dv/dy|
//On a ainsi le plus petit rectangle qui contient l'affinité localement


// on a utilisé cela comme ref pour écrire les fonctions
// c = D[9]x+D[10]y+D[11]
// du/dx = D[1]y+D[2]/c
// du/dy = D[5]x+D[6]/c
// dv/dx = D[3]y+D[4]/c
// dv/dy = D[7]x+D[8]/c


void precal_D(double H[3][3],double *D){
	// 0 1 2
	// 3 4 5
	// 6 7 8
	double *a = H[0];
	D[1]=a[0]*a[7]-a[6]*a[1];
	D[2]=a[0]*a[8]-a[6]*a[2];
	D[3]=a[3]*a[7]-a[6]*a[4];
	D[4]=a[3]*a[8]-a[6]*a[5];
	D[5]=-D[1];
	D[6]=a[1]*a[8]-a[7]*a[2];
	D[7]=-D[3];
	D[8]=a[4]*a[8]-a[7]*a[5];
	D[9]=a[6];
	D[10]=a[7];
	D[11]=a[8];
}


/*double fabs(double p){if(p>0){return p;}{return -p;}} *///une valeur absolu pour les double


//un variable globale que l'on additionne a la distance
#define D_BIAS 0.
#define D_COEFF 1.

//la fonction de calcul de D. On a pris le plus petit rectangle qui encadre le parallélogramme
void cal_D(int x,int y,double *D,double *d,double *coo){
	double p[2]={x,y};
	double a;
	a = pow(D[9]*p[0]+D[10]*p[1]+D[11],2);
	
//On declare les dérivé partielles	
	double dudx,dudy,dvdx,dvdy;
	dudx = (D[1]*p[1]+D[2])/a;
	dudy = (D[5]*p[0]+D[6])/a;
	dvdx = (D[3]*p[1]+D[4])/a;
	dvdy = (D[7]*p[0]+D[8])/a;
	
//On decale l'origine du rectangle pour que le parallélogramme soit à l'intérieur	
	if(dudx<0){coo[0] += dudx*D_COEFF;}
	if(dudy<0){coo[0] += dudy*D_COEFF;}
	if(dvdx<0){coo[1] += dvdx*D_COEFF;}
	if(dvdy<0){coo[1] += dvdy*D_COEFF;}
	
	d[0] = ( fabs(dudx) + fabs(dudy) + D_BIAS )*D_COEFF;
	d[1] = ( fabs(dvdx) + fabs(dvdy) + D_BIAS )*D_COEFF;
}


//l'interpolation bilinéaire
static float bilinear_ripmap(float *x, int w,
		float p, float q, int l, int d1, int d2){	
	float pp = p - 1/2*(pow(2,d1)-1)/pow(2,d1);
	float qq = q - 1/2*(pow(2,d2)-1)/pow(2,d2);   //on decale car les rectangles représente plusieurs pixels
	int ip = floor(p);
	int iq = floor(q);
	float a = x[coord(d1,d2,ip,iq,w,l)];
	float b = x[coord(d1,d2,ip+1,iq,w,l)];
	float c = x[coord(d1,d2,ip,iq+1,w,l)];
	float dd = x[coord(d1,d2,ip+1,iq+1,w,l)];
	return evaluate_bilinear_cell(a, b, c, dd, pp-ip, qq-iq);
}


//la fonction principale. On a distingué beaucoup de cas suivant dh<=1, dh>=logw, etc...

static float ripmap_interpolation_at(float *r, int w, int h,
		float x, float y, int l,double d[2]){
		int logw = (int) log2(w);
		int dh = floor(log2(d[0])) + 1;
		int dl = floor(log2(d[1])) + 1;
		float a,b,c,dd,u,v;
		if((pow(2,dh))>w && (pow(2,dl)>w)){return r[coord(logw,logw,0,0,w,l)];}
		if(dh<=1 && dl<=1){return bilinear_ripmap(r,w,x,y,l,0,0);}
		if(pow(2,dh)>w){
			if(dl<=1){return bilinear_ripmap(r,w,0,y,l,logw,0);}
			{
				a=bilinear_ripmap(r,w,0,y/pow(2,dl-1),l,logw,dl-1);
				b=bilinear_ripmap(r,w,0,y/pow(2,dl-2),l,logw,dl-2);
				return (dl-log2(d[1]))*b + (log2(d[1])-dl+1)*a;
			}
		}
		if(pow(2,dl)>w){
			if(dh<=1){return bilinear_ripmap(r,w,x,0,l,0,logw);}
			{
				a=bilinear_ripmap(r,w,x/pow(2,dh-1),0,l,dh-1,logw);
				b=bilinear_ripmap(r,w,x/pow(2,dh-2),0,l,dh-2,logw);
				return (dh-log2(d[0]))*b + (log2(d[0])-dh+1)*a;
			}
		}
		if(dh<=1){
			a=bilinear_ripmap(r,w,x,y/pow(2,dl-1),l,0,dl-1);
			b=bilinear_ripmap(r,w,x,y/pow(2,dl-2),l,0,dl-2);
			return (dl-log2(d[1]))*b + (log2(d[1])-dl+1)*a;}
		if(dl<=1){
			a=bilinear_ripmap(r,w,x/pow(2,dh-1),y,l,dh-1,0);
			b=bilinear_ripmap(r,w,x/pow(2,dh-2),y,l,dh-2,0);
			return (dh-log2(d[0]))*b + (log2(d[0])-dh+1)*a;
		
		}
		a=bilinear_ripmap(r,w,x/pow(2,dh-1),y/pow(2,dl-1),l,dh-1,dl-1);
		b=bilinear_ripmap(r,w,x/pow(2,dh-2),y/pow(2,dl-1),l,dh-2,dl-1);
		c=bilinear_ripmap(r,w,x/pow(2,dh-1),y/pow(2,dl-2),l,dh-1,dl-2);
		dd=bilinear_ripmap(r,w,x/pow(2,dh-2),y/pow(2,dl-2),l,dh-2,dl-2);
		u=(pow(2,dh)-d[0])/pow(2,dh);
		v=(pow(2,dl)-d[1])/pow(2,dl);
		return u*v*dd + (1-u)*v*b + (1-v)*u*c + (1-u)*(1-v)*a;
	}
	
	
	
	
	
int apply_homo_ripmap(float *img,float *img_f,int w,int h,int w_out,int h_out,double H[3][3]){
	
	double D[12];
	precal_D(H,D);

//on va construire le ripmap de img.
	float *ripmap=malloc(12*w*w*sizeof(float));
	int logw = (int) log2(w);
	build_ripmap(img,ripmap,w,h,3,logw);


	for (int j = 0; j < h_out; j++)
	for (int i = 0; i < w_out; i++)
	{
		double p[2] = {i, j};
		apply_homography_1pt(p, H, p);
		//p[0] = (p[0] - 0.5) * w / (w - 1.0);
		//p[1] = (p[1] - 0.5) * h / (h - 1.0);
		
		double z[2]={i,j};
		double a;
		a = pow(D[9]*z[0]+D[10]*z[1]+D[11],2);	
		double dudx,dudy,dvdx,dvdy;
		dudx = (D[1]*z[1]+D[2])/a;
		dudy = (D[5]*z[0]+D[6])/a;
		dvdx = (D[3]*z[1]+D[4])/a;
		dvdy = (D[7]*z[0]+D[8])/a;
	
		double det = fabs( dudy*dvdx - dvdy*dudx );
		double d[2];
		d[0] = fabs(dudx) + fabs(dudy);
		d[1] = fabs(dvdx) + fabs(dvdy);

		//on doit faire une copie de p, sinon il ne veut pas décaler...
		double coo[2];
		coo[0]=p[0]; coo[1]=p[1];
	
		//On decale l'origine du rectangle pour que le parallélogramme soit à l'intérieur	
		if(dudx<0){coo[0] += dudx;}
		if(dudy<0){coo[0] += dudy;}
		if(dvdx<0){coo[1] += dvdx;}
		if(dvdy<0){coo[1] += dvdy;}
		
		for (int l = 0; l < 3; l++)
		{
			int idx = l + 3 * (w_out * j + i);
			float v = ripmap_interpolation_at(ripmap, w, h, coo[0], coo[1], l, d);
			img_f[idx] = v;
		}
	}
	
	
	//truncate in the end to have non-periodic result
	double p[2];
	for(int i=0;i<w_out;i++){
		for(int j=0;j<h_out;j++){
			p[0]=i; p[1]=j;
			apply_homography_1pt(p,H,p);
			//p[0] = (p[0] - 0.5) * w / (w - 1.0);
			//p[1] = (p[1] - 0.5) * h / (h - 1.0);
			if(p[0]<0 || p[0]>w || p[1]<0 || p[1]>h ){
				for(int l=0;l<3;l++){img_f[(j*w_out+i)*3+l]=0;}
			}
		}
	}
	
	return 0;
		
}
