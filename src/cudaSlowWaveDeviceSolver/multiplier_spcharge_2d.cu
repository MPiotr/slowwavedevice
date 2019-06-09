#include <stdio.h>
#include <complex>
#include "cu_mult.h"
#include "multiplier_spcharge_2d.h"
//#include "print.h"
  

#define PAR ParamsM
#ifndef dm_Pi
__device__ const double dm_Pi = 3.141592653589793;
#endif

__device__ inline double fR_center(unsigned int q_init, unsigned int Nq_max,  double wall)
{
	return 0.5*wall + wall*((double)q_init-0.5*(double)Nq_max)/(double)(Nq_max);
}
__device__ inline double fInitial_angle(double R_center)
{
	return (0.0810194 - 2.05972*R_center +  28.0433*R_center*R_center);
}
__device__ inline double fR_cyclotron(double initial_angle, unsigned int v_init, unsigned int Nv_max)
{
	return (0.568*initial_angle + 0*0.035156*((double)v_init)/double(Nv_max));
}

__global__ void
motionstep_spacecharge_2d(PAR *par, int i, int k, double K, double A, double rB, double iB, int Nharm, double coupling = 1)
{
	
	unsigned int p0 = threadIdx.x;	    unsigned int Np = blockDim.x;
	unsigned int q0 = threadIdx.y;		unsigned int Nq = blockDim.y;
	unsigned int s0 = threadIdx.z;		unsigned int Ns = blockDim.z;

	unsigned int X = blockIdx.x;		unsigned int NX = gridDim.x;
	unsigned int Y = blockIdx.y;		unsigned int NY = gridDim.y;
	unsigned int Z = blockIdx.z;	//	unsigned int NZ = gridDim.z;

	unsigned int q_init =		 Nq*X + q0;		unsigned int Nq_max =		 Nq*gridDim.x;
	unsigned int s_init =		 Ns*Y + s0;		unsigned int Ns_max =		 Ns*gridDim.y;
	unsigned int v_init =			Z;			unsigned int Nv_max =		    gridDim.z;

	unsigned int warpsize = Ns*Nq*Np;
	unsigned int log2warpsize = round(log2((double)warpsize));
	int xi = (p0 + Np*q0 + s0*Np*Nq);



	double h, k1, voltage, g1, g3;
	h  = par->h; k1 = par->k1; 
	g1 = par->g1; voltage = par->voltage;
	g3 = par->g3;
	

	int N = par->Nz; 
	double dz =  par->L/(double)N;

	double z;
	
	double Q, Qk;
	double W, Wk;

	double rA = A, iA = 0; 
	double rQ1, iQ1, rQ3, iQ3, r;
	
	double R_cyclotron, R_center, kappa_cyclotron, phase_cyclotron, initial_angle;
	double wall = par->wall;
	
	R_center =  fR_center(q_init, Nq_max, wall);
	initial_angle = fInitial_angle(R_center);
	R_cyclotron =	fR_cyclotron(initial_angle, v_init, Nv_max);
	kappa_cyclotron = 1.758;
	phase_cyclotron = 2.*dm_Pi*(double)s_init/(double)Ns_max;

	double en0 = 1. + voltage/511.;	
	double angle_spread_factor = 1./sqrt(1. + initial_angle*initial_angle);




	__shared__ double    sh_sinQ[NS][NQ][NP];
	__shared__ double    sh_cosQ[NS][NQ][NP];
	__shared__ double sh_rAk[NR];
	__shared__ double sh_iAk[NR];

	double PH, EN, cosPH, sinPH, sinPS, cosPS;



	double H = h;//+dh(delta);
//	double NU = 1./(en0*(en0*en0 - 1.));
	double DELTA = 1-k1*en0/sqrt(en0*en0 - 1.);

	double *d_Qk = par->Qk;	double *Q0 = par->Q0;
	double *d_Wk = par->Wk;	double *W0 = par->W0;
	double *rAq3k = par->rAk;	double *iAq3k = par->iAk;
	double *rAq1k = par->rAq1k;	double *iAq1k = par->iAq1k;
	double *radii = par->radii;

	int gstride = (NX*NY*Z + NX*Y+X)*warpsize;

	int *destroyed = par->ifnotdestroyed;
	int ifnotdestroyed = destroyed[gstride + xi];

	
	

	if(i != 0)
	{
		Q = Q0[gstride + xi];
		W = W0[gstride + xi];
	}
	else
	{
		Q = 2.*dm_Pi/double(Np)*double(p0);
		W = 0; 
		Q0[gstride + xi] = Q;
		W0[gstride + xi] = W;

	}

	if( (i == 0)&&(xi ==  0))
	{
		par->avEN[NX*NY*Z + NX*Y+X] = 0; 
		par->int_rJ3_1[NX*NY*Z + NX*Y+X] = 0;
		par->int_iJ3_1[NX*NY*Z + NX*Y+X] = 0;
		par->int_rQ1[NX*NY*Z + NX*Y+X] = 0;
		par->int_iQ1[NX*NY*Z + NX*Y+X] = 0;
	}
	__syncthreads();


	Wk = (k > 0) ? d_Wk[4*gstride + warpsize*(k-1) + xi] : 0;
	Qk = (k > 0) ? d_Qk[4*gstride + warpsize*(k-1) + xi] : 0;

	z = ((double)i+ K)*dz; 
	r = (R_center + R_cyclotron*cos(kappa_cyclotron*z + phase_cyclotron));
	//ifnotdestroyed *= (r > -wall)? 1. : 0.;
	//ifnotdestroyed = 1;
	
	PH = Q + K*Qk;
	EN = W + K*Wk + en0;

	int ind_R = floor((r + wall)/(3.*wall)*NR);
	float dR = 3.*wall/(float)NR;
	double anchR = ind_R*dR - wall;

	double rDivDr = (r - anchR)/dR;
	
	if(ind_R < 0 ) ind_R = 0;

	double rak, iak, iq0, rq0;


	 
	if(k == 0) {	rak = 0; iak = 0;	} else												
	{																							//Загрузка (и интерполяция) амплитуд пространственного заряда
		if(xi < NR)    sh_rAk[xi] = rAq3k[xi+NR*(k-1)]; else									//Амплитуда пространственного заряда описывается таблицей						
		{if(xi < 2*NR) sh_iAk[xi-NR] = iAq3k[xi-NR+NR*(k-1)];}										//длиной NR и охватывает радиусы от -wall до 2*wall (в местных координатах)
		__syncthreads();																			//Отдельно берём значения k-коэффициента, отдельно само значение поля (для 1-й и 
		if(ind_R < NR - 1)																			//3-й гармоник (чтобы сэкономить shared memory)
		{																						
			rak = sh_rAk[ind_R] + (sh_rAk[ind_R+1]-sh_rAk[ind_R])*rDivDr;
			iak = sh_iAk[ind_R] + (sh_iAk[ind_R+1]-sh_iAk[ind_R])*rDivDr;
		}
		else
		{
			rak = sh_rAk[NR*(k-1) + ind_R];
			iak = sh_iAk[NR*(k-1) + ind_R];
		}
		__syncthreads();

	}

	
	

	if((xi < NR))    sh_rAk[xi] = par->ar0[xi];
	if((xi < 2*NR)&&(xi >=NR)) sh_iAk[xi-NR] = par->ai0[xi-NR];

	__syncthreads();												

	if(ind_R < NR - 1)
	{
		rq0 = sh_rAk[ind_R] + (sh_rAk[ind_R+1]-sh_rAk[ind_R])*rDivDr;
		iq0 = sh_iAk[ind_R] + (sh_iAk[ind_R+1]-sh_iAk[ind_R])*rDivDr;
	}
	else
	{
		rq0 = sh_rAk[ind_R];
		iq0 = sh_iAk[ind_R];
	}

	
	rQ3 = (rq0 + K*rak);
	iQ3 = (iq0 + K*iak);

//	if((i < 2)&&(p0 == 0)) printf("%i,%i, %i\t%g\t%g\t%g\t%g\n",i, k, q0, rQ3, iQ3,  iq0, iak);

	

	if(k == 0) {	rak = 0; iak = 0;	} else												//Амплитуда пространственного заряда описывается таблицей
	{																						//длиной NR и охватывает радиусы от -wall до 2*wall (в местных координатах)
		if(xi < NR)    sh_rAk[xi] = rAq1k[xi+NR*(k-1)]; else								//Отдельно берём значения k-коэффициента, отдельно само значение поля (для 1-й и 
		{if(xi < 2*NR) sh_iAk[xi - NR] = iAq1k[xi - NR +NR*(k-1)];}									//3-й гармоник (чтобы сэкономить shared memory)
		__syncthreads();																	
		
		if(ind_R < NR - 1)
		{
			rak = sh_rAk[ind_R] + (sh_rAk[ind_R+1]-sh_rAk[ind_R])*rDivDr;
			iak = sh_iAk[ind_R] + (sh_iAk[ind_R+1]-sh_iAk[ind_R])*rDivDr;
		}
		else
		{
			rak = sh_rAk[NR*(k-1) + ind_R];
			iak = sh_iAk[NR*(k-1) + ind_R];
		}
	}

	if(xi < NR)    sh_rAk[xi] = par->rAq1[xi]; else
	{if(xi < 2*NR) sh_iAk[xi - NR] = par->iAq1[xi - NR];}		
	__syncthreads();												

	if(ind_R < NR - 1)
	{
		rq0 = sh_rAk[ind_R] + (sh_rAk[ind_R+1]-sh_rAk[ind_R])*rDivDr;
		iq0 = sh_iAk[ind_R] + (sh_iAk[ind_R+1]-sh_iAk[ind_R])*rDivDr;
	//	if((i < 2)&&(xi == 0)) printf("%g\t%g\t%g\t%g  %g %g\n", sh_iAk[ind_R],sh_iAk[ind_R+1], iq0,  r, anchR, dR);
	//	printf("%i,%i\t%i,%i\t%g,%g\n", i,k, s0, ind_R, sh_rAk[ind_R] , (sh_rAk[ind_R+1]-sh_rAk[ind_R])/rDivDr);
	}
	else
	{
		rq0 = sh_rAk[ind_R];
		iq0 = sh_iAk[ind_R];
		
	}

	rQ1 = (rq0 + K*rak);
	iQ1 = (iq0 + K*iak);

//	if((i < 2)&&(p0 == 0)&&(q0 < 2)) printf("%i,%i,%i\t%g\t%g\t%g\t%g\n",i, k, q0, rQ1, iQ1,  iq0, iak);

	
	sincos(PH, &sinPH, &cosPH);
	sincos(3.*PH, &sinPS, &cosPS);

	double DQ = (H - k1*angle_spread_factor*EN/sqrt(EN*EN-1.));
	rA *= exp(-g1*r);
	double tmp = exp(-g3*r);
	rB *= tmp;
	iB *= tmp;

	
	
	Qk = dz*DQ;
//	Qk = dz*(H*NU*W + DELTA);
	Wk = -dz*((rA - iQ1)*cosPH - (iA + rQ1)*sinPH + (rB - iQ3)*cosPS -(iB +  rQ3)*sinPS)*double(ifnotdestroyed);

		
	d_Qk[4*gstride + warpsize*k + xi] = Qk;
	d_Wk[4*gstride + warpsize*k + xi] = Wk;

///	printf("%i,%i\t%i,%i\t%g,%g\n", i,k, s0, q0, d_Qk[4*gstride + warpsize*k + xi], d_Wk[4*gstride + warpsize*k + xi]);

	

	

	

	unsigned int stride;
	unsigned int log2Np = round(log2((float)Np));
	int ind_r0 = (q0  + Nq*s0) + Nq*Ns*(X + NX*Y+NX*NY*Z);//(q0 + Nq*X) + (Nq*X)*(s0 + Ns*Y) + (Nq*NX* Ns*NY)*Z;
//	int tot_r0 = Nq*NX * Ns*NY *NZ;

    

	sh_sinQ[s0][q0][p0] = DQ*sinPS*double(ifnotdestroyed);
	sh_cosQ[s0][q0][p0] = DQ*cosPS*double(ifnotdestroyed);

	__syncthreads();
	

	
	
	stride = Np;
	for(int q = 1; q <= log2Np; q++)
	{
		stride = stride >> 1;
		if(p0 < stride)
		{
			sh_sinQ[s0][q0][p0] += sh_sinQ[s0][q0][p0 + stride];
			
		}
		else
		{
			if(p0 < 2*stride)
			{
				
				sh_cosQ[s0][q0][p0 - stride] += sh_cosQ[s0][q0][p0];
				
			}
		}
		__syncthreads();	
	}


	if(p0 == 0)
	{
		
		par->int_rJ3[ind_r0] =  -3.*sh_sinQ[s0][q0][0]/Np; 
		par->int_iJ3[ind_r0] =  -3.*sh_cosQ[s0][q0][0]/Np; 
		
	
	}
	__syncthreads();

	sh_sinQ[s0][q0][p0] = DQ*sinPH*double(ifnotdestroyed);
	sh_cosQ[s0][q0][p0] = DQ*cosPH*double(ifnotdestroyed);
	
	
	stride = Np;
	for(int q = 1; q <= log2Np; q++)
	{
		stride = stride >> 1;
		if(p0 < stride)
		{
			sh_sinQ[s0][q0][p0] += sh_sinQ[s0][q0][p0 + stride];
			
		}
		else
		{
			if(p0 < 2*stride)
			{
				
				sh_cosQ[s0][q0][p0 - stride] += sh_cosQ[s0][q0][p0];
				
			}
		}
		__syncthreads();	
	}


	if(p0 == 0)
	{
		par->int_rQ1[ind_r0] = -sh_sinQ[s0][q0][0]/Np;
		par->int_iQ1[ind_r0] = -sh_cosQ[s0][q0][0]/Np;
		radii[ind_r0] = r;
//		if(i < 5) printf("motion out: %i,%i\t%i, %i\t %g, %g\n", i,k, s0, q0, par->int_rQ1[ind_r0], par->int_iQ1[ind_r0]);
	
	}

	if(X == 0)
	{
		if((xi < NR))					rAq3k[xi + NR*k] = 0; 
		if((xi < 2*NR)&&(xi >= NR))		iAq3k[xi - NR+NR*k] = 0;
		if((xi < 3*NR)&&(xi >= 2*NR))	rAq1k[xi-2*NR+NR*k] = 0; 
		if((xi < 4*NR)&&(xi >= 3*NR))	iAq1k[xi-3*NR+NR*k] = 0;
	
	}


}
__global__ void
amplitudestep_spacecharge_2d(PAR *par, int i,  int k, double K, double spchQ1, double spchQ3)
{
	unsigned int xi = threadIdx.x;
	unsigned int warpsize = blockDim.x;

//	unsigned int ind = xi + warpsize*blockIdx.x;	
//	unsigned int tot = warpsize*gridDim.x; 
	
	unsigned int log2warpsize = round(log2((float)warpsize));

	double wall = par->wall;

	double *rAq3k = par->rAk;
	double *iAq3k = par->iAk;
	double *rAq1k = par->rAq1k;
	double *iAq1k = par->iAq1k;

	double *rQ3 = par->int_rJ3;
	double *iQ3 = par->int_iJ3;

	double *rQ1 = par->int_rQ1;
	double *iQ1 = par->int_iQ1;

	int N = par->Nz;
	double dz = par->L/(double)N;
	double H = par->h;
	double r0, dR, FQ;
	double tmp_rQ, tmp_iQ;

	

	__shared__ double   sh_R[1024]; 
	__shared__ double   sh_I[1024]; 
	__shared__ float	sh_radii[1024]; // попытка впихнуть всё в shared память. Пока реально надо 1024 токи. Память позволяет увеличить их в два раза

	for(int m = 0; m < gridDim.x; m++) sh_radii[xi] = par->radii[xi + warpsize*m]; //точек 1024 (или может больше), а потоков только 512. 

	tmp_rQ = rQ3[xi]/warpsize;
	tmp_iQ = iQ3[xi]/warpsize;

	
	__syncthreads();

	dR = 3.*wall/(double) NR;

									//Здесь происходит вычисление амплитуды простраственного заряда
	for(int nr = 0; nr < NR; nr++)
	{
		r0 = dR*nr - wall;
		FQ = (sh_radii[xi] > r0) ? (sinh(3.*H*(r0 + wall))*exp(-3.*H*(sh_radii[xi]+wall))/(sh_radii[xi]+0.7-wall)) : (sinh(3.*H*(sh_radii[xi] + wall))*exp(-3.*H*(r0+wall))/(sh_radii[xi]+0.7-wall)); 

		sh_R[xi] = tmp_rQ*FQ;	//надо сначала скопировать все значения в shared память
		sh_I[xi] = tmp_iQ*FQ;

				

		int stride = warpsize;
		__syncthreads();
	
		for(int q = 1; q <= log2warpsize; q++)
		{
			stride = stride >> 1;
			if(xi < stride)
			{
				
				sh_R[xi] += sh_R[xi + stride];
			}
			else
			{
				if(xi < 2*stride)
				{
					sh_I[xi - stride] += sh_I[xi];
				}
			}
			__syncthreads();
		}

		
	

		if(xi == 0) rAq3k[NR*k + nr] += dz*spchQ3*sh_R[0]/gridDim.x; // каждый блок вносит свой вклад. Вклады складываются.
		if(xi == 0) iAq3k[NR*k + nr] += dz*spchQ3*sh_I[0]/gridDim.x;
				
		
		__syncthreads(); 
	}

	tmp_rQ = rQ1[xi]/warpsize;
	tmp_iQ = iQ1[xi]/warpsize;
	__syncthreads();

//	if(i < 5) printf("amplitude in: %i\t %g, %g\n", xi, rQ1[xi], iQ1[xi]);

	for(int nr = 0; nr < NR; nr++)
	{
		r0 = dR*nr - wall;
		FQ = (sh_radii[xi] > r0) ? (sinh(H*(r0 + wall))*exp(-H*(sh_radii[xi]+wall))/(sh_radii[xi]+0.7-wall)) : (sinh(H*(sh_radii[xi] + wall))*exp(-H*(r0+wall))/(sh_radii[xi]+0.7-wall)); 

	//	if((i == 0)) printf("r = %g\t%i\t%g\t%g\t%g\n", r0, xi, sh_radii[xi], FQ*spchQ1, FQ);

		sh_R[xi] = tmp_rQ*FQ;
		sh_I[xi] = tmp_iQ*FQ;

	
        __syncthreads();

		int stride = warpsize;

		
		for(int q = 1; q <= log2warpsize; q++)
		{
			stride = stride >> 1;
			if(xi < stride)
			{
				sh_R[xi] += sh_R[xi + stride];
			}
			else
			{
				if(xi < 2*stride)
				{
					sh_I[xi - stride] += sh_I[xi];
				}
			}	
			__syncthreads();
		}
	
		if(xi == 0)	rAq1k[NR*k + nr] += dz*spchQ1*sh_R[0]/gridDim.x;
		if(xi == 0) iAq1k[NR*k + nr] += dz*spchQ1*sh_I[0]/gridDim.x; 
	//	if(((xi == 0)||(xi == NP))&&(i < 2)) printf("%i\t%i\t%i\t%i\t%g\t%g\n", i,  nr, k,	xi, rAq1k[NR*k + nr], iAq1k[NR*k + nr]);
		
		
	
		
	}



}
__global__ void
spacecharge_2d_endstep(PAR *par, int i, int Nharm = 3, double fB = 0, int N_en_averaging = -1)
{
	unsigned int p0 = threadIdx.x;	    unsigned int Np = blockDim.x;
	unsigned int q0 = threadIdx.y;		unsigned int Nq = blockDim.y;
	unsigned int s0 = threadIdx.z;		unsigned int Ns = blockDim.z;

	unsigned int X_ = blockIdx.x;		unsigned int NX = gridDim.x;
	unsigned int Y_ = blockIdx.y;		unsigned int NY = gridDim.y;
	unsigned int Z_ = blockIdx.z;	//	unsigned int NZ = gridDim.z;

	unsigned int q_init =		 Nq*X_ + q0;		unsigned int Nq_max =		 Nq*gridDim.x;
	unsigned int s_init =		 Ns*Y_ + s0;		unsigned int Ns_max =		 Ns*gridDim.y;
	unsigned int v_init =			Z_;			unsigned int Nv_max =		    gridDim.z;

	double wall = par->wall; 

	double R_center =  fR_center(q_init, Nq_max, wall);
	double initial_angle = fInitial_angle(R_center);
	double R_cyclotron =	fR_cyclotron(initial_angle, v_init, Nv_max);
	double kappa_cyclotron = 1.758;
	double phase_cyclotron = 2.*dm_Pi*(double)s_init/(double)Ns_max;
	
	double z = par->L/double(par->Nz)*double(i);
	double r = (R_center + R_cyclotron*cos(kappa_cyclotron*z + phase_cyclotron));
	double g = (Nharm == 1)? (par->g1) : (par->g3);

	unsigned int warpsize = Ns*Nq*Np;
	unsigned int log2warpsize = round(log2((double)warpsize));

	int xi = (p0 + Np*q0 + s0*Np*Nq);

	
//	unsigned int xi = threadIdx.x;	    unsigned int Np = blockDim.x;
	unsigned int X = X_+NX*Y_+NX*NY*Z_;	

//	unsigned int warpsize = Np;
//	unsigned int gridsize = NX*NY*NZ;
	int Nz = par->Nz;
	
	double *d_Qk = par->Qk; double *Q0 = par->Q0;
	double *d_Wk = par->Wk; double *W0 = par->W0;
	double *rAq3k = par->rAk;   double *iAq3k = par->iAk;
	double *rAq1k = par->rAq1k; double *iAq1k = par->iAq1k;
	double *rJ3 = par->rJ3; double *iJ3 = par->iJ3;

	unsigned int stride = 4*X*warpsize;
	unsigned int gstride = X*warpsize;

	if(r < -wall) par->ifnotdestroyed[gstride + xi] *= 0;
	int ifnotdestroyed = par->ifnotdestroyed[gstride + xi]; 
	
	Q0[gstride + xi] += (d_Qk[stride + xi] + 2.*d_Qk[stride+ xi + warpsize] + 2.*d_Qk[stride+ xi + 2*warpsize] + d_Qk[stride+ xi+ 3*warpsize])/6.;
	W0[gstride + xi] += (d_Wk[stride + xi] + 2.*d_Wk[stride +xi + warpsize] + 2.*d_Wk[stride +xi + 2*warpsize] + d_Wk[stride +xi+ 3*warpsize])/6.;

	
	if(X == 0)
	{
		if((xi < NR))					par->ar0[xi]        += (rAq3k[xi]		 + 2.*rAq3k[xi + NR]		+ 2.*rAq3k[xi  + NR*2]		 + rAq3k[xi + NR*3])/6.;
		if((xi < 2*NR)&&(xi >= NR))		par->ai0[xi - NR]   += (iAq3k[xi - NR]	 + 2.*iAq3k[xi - NR + NR]	+ 2.*iAq3k[xi  - NR + NR*2]	 + iAq3k[xi - NR + NR*3])/6.;
		if((xi < 3*NR)&&(xi >= 2*NR))   par->rAq1[xi - 2*NR]+= (rAq1k[xi - 2*NR] + 2.*rAq1k[xi - 2*NR + NR] + 2.*rAq1k[xi - 2*NR + NR*2] + rAq1k[xi - 2*NR + NR*3])/6.; 
		if((xi < 4*NR)&&(xi >= 3*NR))   par->iAq1[xi - 3*NR]+= (iAq1k[xi - 3*NR] + 2.*iAq1k[xi - 3*NR + NR]	+ 2.*iAq1k[xi - 3*NR + NR*2] + iAq1k[xi - 3*NR + NR*3])/6.;
//		if(xi == 0) printf("%g ... ", par->rAq1[xi]);
	}
	
	__shared__ double sh_sinQ[NS*NQ*NP];
	__shared__ double sh_cosQ[NS*NQ*NP];

	double sinPS, cosPS;

	sincos((double)Nharm*Q0[gstride + xi], &sinPS, &cosPS);			//Это общий ток. Усреднение по всем частицам.

	
	sh_sinQ[xi] = sinPS*exp(-g*r)*double(ifnotdestroyed);
	sh_cosQ[xi] = cosPS*exp(-g*r)*double(ifnotdestroyed);

	__syncthreads();	

	stride = warpsize;
	for(int q = 1; q <= log2warpsize; q++)
	{
		stride = stride >> 1; 
		if(xi < stride)
		{
			sh_sinQ[xi] += sh_sinQ[xi + stride];
			
		}
		else
		{
			if(xi < 2*stride)
			{
				
				sh_cosQ[xi - stride] += sh_cosQ[xi];
				
			}
		}
		__syncthreads();	
	}
	
	__syncthreads();	
	if(xi == 0)
	{
		rJ3[X*Nz + i] =  sh_cosQ[0]/warpsize;
		iJ3[X*Nz + i] = -sh_sinQ[0]/warpsize;
		

		par->int_rJ3_1[X] +=     sh_cosQ[0]/warpsize*fB;
		par->int_iJ3_1[X] += -1.*sh_sinQ[0]/warpsize*fB; 

	}
	__syncthreads();

	if(i == N_en_averaging)
	{
		

		sh_sinQ[xi] = W0[gstride + xi]*double(ifnotdestroyed);

		__syncthreads();

		stride = warpsize;
		for(int q = 1; q <=  log2warpsize; q++)
		{
			stride =  stride >> 1;
			if(xi < stride)
			{
				sh_sinQ[xi] += sh_sinQ[xi + stride];
			}

			__syncthreads();	
		}

		if(xi == 0)
		{
			par->avEN[X] = sh_sinQ[0]/(warpsize);
		}
		
	}

	
	
	

}


double  Multiplier_SpCharge_2D::DeltaEnergySpaceCharge2D(double A)
{
	int GQ = Nq/NQ; int GS = Ns/NS; int GV = Nv;

//	int warpsize = NQ*NS*NP;
	int gridsize = GQ*GS*GV;

	ParamsM par = setPar();

	double K[4] = {0, 0.5, 0.5, 1.};

	dim3 motionwarp = dim3(NP, NQ, NS);
	dim3 motiongrid = dim3(GQ, GS, GV);

	dim3 amplitudewarp = dim3(NQ, NS, GV);
	dim3 amplitudegrip = dim3(GQ, GS, 1);

	cudaMemset(d_Qk, 0, 4*sizeof(double)*Np*Nq*Ns*Nv);	cudaMemset(d_W0, 0, sizeof(double)*Np*Nq*Ns*Nv);
	cudaMemset(d_Wk, 0, 4*sizeof(double)*Np*Nq*Ns*Nv);	cudaMemset(d_rAk, 0, sizeof(double)*4*NR);
	cudaMemset(d_Q0, 0, sizeof(double)*Np*Nq*Ns*Nv);	cudaMemset(d_iAk, 0, sizeof(double)*4*NR);
	cudaMemset(d_rAq1k, 0, sizeof(double)*4*NR);
	cudaMemset(d_iAq1k, 0, sizeof(double)*4*NR);
	
	
	
	cudaMemset(d_ai0, 0, sizeof(double)*NR);
	cudaMemset(d_ar0, 0, sizeof(double)*NR);
	cudaMemset(d_rAq1, 0, sizeof(double)*NR);	
	cudaMemset(d_iAq1, 0, sizeof(double)*NR);	



	cudaMemset(d_avEN, 0, sizeof(double));


	double la = Nperiods*period;
	int Nz = Nmax;
	double dz = Lmax/double(Nz);

	int Nstop = ceil(la/dz);
	double z,  rB = 0, iB = 0;
	double res = 0; 
	double avEN[512];	

	cudaMemcpy(d_par, &par, sizeof(PAR), cudaMemcpyHostToDevice);

//	printf("la = %g, dz = %g, Nstop = %i, Nz = %i, GQ*GS*GV = %i\n", la, dz, Nstop, Nz, GQ*GS*GV);
//	printf("warpsize = %i, gridsize = %i\n", warpsize, gridsize);

	double *reQ1 = new double [NR*Nstop];
	double *imQ1 = new double [NR*Nstop];
//	FILE *file = fopen("F:\\Piotr\\Multiplier_260GHz_Data\\q2d.txt","w");
		
	for(int i = 0; i < Nstop; i++)
	{
		z = double(i)*dz;

		for(int k = 0; k < 4; k++)
		{
			motionstep_spacecharge_2d <<<motiongrid,motionwarp>>>(d_par, i, k, K[k],  A*sin(Pi/la*(z+K[k]*dz)), rB, iB, 1);
			amplitudestep_spacecharge_2d  <<</*Nq*Ns*Nv/512 + */1, Nq*Ns*Nv>>>(d_par, i, k, K[k], spchQ1, spchQ3);
		}

		spacecharge_2d_endstep<<<motiongrid, motionwarp>>>(d_par, i, 3, 0, Nstop-1);
		
	/*	cudaMemcpy(reQ1 + i*NR, d_rAq1, sizeof(double)*NR, cudaMemcpyDeviceToHost);
		cudaMemcpy(imQ1 + i*NR, d_iAq1, sizeof(double)*NR, cudaMemcpyDeviceToHost);
		for(int a = 0; a < NR; a++) {double r = 3.*wall/double(NR)*(a) - wall; fprintf(file, "%g,%g,%g, %g\n", z, r, reQ1[a+i*NR], imQ1[a+i*NR]);}*/

		
	}
	res = 0;
	cudaMemcpy(&avEN, (void*) d_avEN, sizeof(double)*gridsize, cudaMemcpyDeviceToHost);
	
	for(int q = 0; q < gridsize; q++) 
	{
		res += avEN[q];
	//	printf("avEn = %g\n", avEN[q]);
	}

	res /= gridsize;
/*
	fclose(file);
	delete[] reQ1;
	delete[] imQ1;*/

	return res;
	


}
std::complex<double>  Multiplier_SpCharge_2D::CurrentASpaceCharge2D(double reA, double imA)
{
		int GQ = Nq/NQ; int GS = Ns/NS; int GV = Nv;
	

//	int warpsize = NP*NQ*NS;
	int gridsize = GQ*GS*GV;

	ParamsM par = setPar();

	double K[4] = {0, 0.5, 0.5, 1.};

	dim3 motionwarp = dim3(NP, NQ, NS);
	dim3 motiongrid = dim3(GQ, GS, GV);

	cudaMemset(d_Qk, 0, 4*sizeof(double)*Np*Nq*Ns*Nv);	cudaMemset(d_W0, 0, sizeof(double)*Np*Nq*Ns*Nv);
	cudaMemset(d_Wk, 0, 4*sizeof(double)*Np*Nq*Ns*Nv);	cudaMemset(d_rAk, 0, sizeof(double)*4*NR);
	cudaMemset(d_Q0, 0, sizeof(double)*Np*Nq*Ns*Nv);	cudaMemset(d_iAk, 0, sizeof(double)*4*NR);
	cudaMemset(d_rAq1k, 0, sizeof(double)*4*NR);
	cudaMemset(d_iAq1k, 0, sizeof(double)*4*NR);
	

	cudaMemset(d_ai0, 0, sizeof(double)*NR);
	cudaMemset(d_ar0, 0, sizeof(double)*NR);
	cudaMemset(d_rAq1, 0, sizeof(double)*NR);	
	cudaMemset(d_iAq1, 0, sizeof(double)*NR);	

	cudaMemset(d_avEN, 0, sizeof(double));

	double La = Nperiods*period;
	double dz = Lmax/double(Nmax);

//	double *energy = new double [Np*Nq*Ns*Nv];
//	double *phase = new double [Np*Nq*Ns*Nv];
//	FILE *res_enr = fopen("F:\\Piotr\\bwo_Data\\enPhase.csv", "w");
	
	
	int Nstop = ceil(La/dz);
	int nd = round(period/dz);

	double z, fA, fB;
//	double debug[3000];

	cudaMemcpy(d_par, &par, sizeof(PAR), cudaMemcpyHostToDevice);
/*
	double *reQ1 = new double [NR*Nstop];
	double *imQ1 = new double [NR*Nstop];	
	FILE *file = fopen("F:\\Piotr\\Multiplier_260GHz_Data\\j1d.txt","w");*/

		
	double absA = sqrt(reA*reA + imA*imA);
	double argA = atan(imA/reA);
	for(int i = 0; i < Nstop; i++)
	{
		z = double(i)*dz;
		//if(z < La*0.5) fA = 0.6*sin(2.*Pi/La*z) + 0.08; else sin(2.*Pi/La*(z - La*0.5))+0.08;
		fA = sin(Pi/La*z);
		fA *= absA;
		if(z > La) fA = 0;
		if((z > La+Ld)&&(z < La+ Ld+ Lb)) fB = sin(Pi*(z - La - Ld)/Lb); else fB = 0;


		for(int k = 0; k < 4; k++)
		{
			motionstep_spacecharge_2d <<<motiongrid,motionwarp>>>(d_par, i, k, K[k], fA, 0, 0, 1);
			amplitudestep_spacecharge_2d  <<<1, Nq*Ns*Nv>>>(d_par, i, k, K[k], spchQ1, spchQ3);
		}
//		spacecharge_endstep<<<gridsize, warpsize>>>(d_par, i, fB);
		spacecharge_2d_endstep<<<motiongrid, motionwarp>>>(d_par, i, 1, sin(Pi/La*z));

/*		cudaMemcpy(reQ1 + i*NR, d_ar0, sizeof(double)*NR, cudaMemcpyDeviceToHost);
		cudaMemcpy(imQ1 + i*NR, d_ai0, sizeof(double)*NR, cudaMemcpyDeviceToHost);
		for(int a = 0; a < NR; a++) {double r = 3.*wall/double(NR)*(a) - wall; fprintf(file, "%g,%g,%g, %g\n", z, r, reQ1[a+i*NR], imQ1[a+i*NR]);}*/

/*		if((enPrint != 0)&&(i%nd == 0))
		{
			cudaMemcpy(energy, d_W0, sizeof(double)*Np*Nq*Ns*Nv, cudaMemcpyDeviceToHost);
			cudaMemcpy(phase, d_Q0, sizeof(double)*Np*Nq*Ns*Nv, cudaMemcpyDeviceToHost);
			for(int en = 0; en < Np*Nq*Ns*Nv; en++)
			{
				fprintf(res_enr, "%g,%g\n", fmod(phase[en]+ 0.5*Pi, 2.*Pi) , energy[en]);

			}

		}*/ 

	}
		
	double rJ = 0, iJ = 0;	
	double rJ_1[6000];
	double iJ_1[6000];

	cudaMemcpy(rJ_1, (void*) d_int_rJ3_1, gridsize*sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(iJ_1, (void*) d_int_iJ3_1, gridsize*sizeof(double), cudaMemcpyDeviceToHost);

/*	cudaMemcpy(rJ_1, (void*) d_rJ3, gridsize*sizeof(double)*Nstop, cudaMemcpyDeviceToHost);
  	cudaMemcpy(iJ_1, (void*) d_iJ3, gridsize*sizeof(double)*Nstop, cudaMemcpyDeviceToHost);

	for(int a = 0; a < Nstop; a++) fprintf(file, "%g,%g,%g\n", a*dz, rJ_1[a]*sin(Pi/La*a*dz), iJ_1[a]*sin(Pi/La*a*dz));*/

	for(int i = 0; i < gridsize; i++)
	{
		rJ += rJ_1[i]; iJ += iJ_1[i];
	//	printf("rJ[%i] = %g\tiJ[%i] = %g\n", i,rJ_1[i], i, iJ_1[i]);
	}

	
	//	printf("rJ = %g\tiJ = %g\n", rJ/gridsize*dz, iJ/gridsize*dz);


/*	delete [] energy;
	delete [] phase;
	fclose(res_enr);*/

//	fclose(file);
/*	delete[] reQ1;
	delete[] imQ1;*/


	std::complex<double> res = std::complex<double>(rJ/gridsize*dz, iJ/gridsize*dz);
	res *= exp(I*arg(reA + I*imA));

	return res;

}
std::complex<double>  Multiplier_SpCharge_2D::CurrentBSpaceCharge2D(double rB, double iB, double A, double enPrint)
{
	int GQ = Nq/NQ; int GS = Ns/NS; int GV = Nv;
	

//	int warpsize = NP*NQ*NS;
	int gridsize = GQ*GS*GV;

	ParamsM par = setPar();

	double K[4] = {0, 0.5, 0.5, 1.};

	dim3 motionwarp = dim3(NP, NQ, NS);
	dim3 motiongrid = dim3(GQ, GS, GV);

	cudaMemset(d_Qk, 0, 4*sizeof(double)*Np*Nq*Ns*Nv);	cudaMemset(d_W0, 0, sizeof(double)*Np*Nq*Ns*Nv);
	cudaMemset(d_Wk, 0, 4*sizeof(double)*Np*Nq*Ns*Nv);	cudaMemset(d_rAk, 0, sizeof(double)*4*NR);
	cudaMemset(d_Q0, 0, sizeof(double)*Np*Nq*Ns*Nv);	cudaMemset(d_iAk, 0, sizeof(double)*4*NR);
	cudaMemset(d_rAq1k, 0, sizeof(double)*4*NR);
	cudaMemset(d_iAq1k, 0, sizeof(double)*4*NR);
	

	cudaMemset(d_ai0, 0, sizeof(double)*NR);
	cudaMemset(d_ar0, 0, sizeof(double)*NR);
	cudaMemset(d_rAq1, 0, sizeof(double)*NR);	
	cudaMemset(d_iAq1, 0, sizeof(double)*NR);	

	cudaMemset(d_avEN, 0, sizeof(double));

	double La = Nperiods*period;
	double dz = Lmax/double(Nmax);

//	double *energy = new double [Np*Nq*Ns*Nv];
//	double *phase = new double [Np*Nq*Ns*Nv];
//	FILE *res_enr = fopen("F:\\Piotr\\bwo_Data\\enPhase.csv", "w");
	
	int Nstop = ceil((La + Lb + Ld)/dz);
	int nd = round(period/dz);

	double z, fA, fB;
//	double debug[3000];

	cudaMemcpy(d_par, &par, sizeof(PAR), cudaMemcpyHostToDevice);
/*
	double *reQ1 = new double [NR*Nstop];
	double *imQ1 = new double [NR*Nstop];	
	FILE *file = fopen("F:\\Piotr\\Multiplier_260GHz_Data\\q2d.txt","w");*/

		
	for(int i = 0; i < Nstop; i++)
	{
		z = double(i)*dz;
		//if(z < La*0.5) fA = 0.6*sin(2.*Pi/La*z) + 0.08; else sin(2.*Pi/La*(z - La*0.5))+0.08;
		fA = sin(Pi/La*z);
		fA *= A;
		if(z > La) fA = 0;
		if((z > La+Ld)&&(z < La+ Ld+ Lb)) fB = sin(Pi*(z - La - Ld)/Lb); else fB = 0;


		for(int k = 0; k < 4; k++)
		{
			motionstep_spacecharge_2d <<<motiongrid,motionwarp>>>(d_par, i, k, K[k], fA, rB*fB, iB*fB, 1);
			amplitudestep_spacecharge_2d  <<<1, Nq*Ns*Nv>>>(d_par, i, k, K[k], spchQ1, spchQ3);
		}
//		spacecharge_endstep<<<gridsize, warpsize>>>(d_par, i, fB);
		spacecharge_2d_endstep<<<motiongrid, motionwarp>>>(d_par, i, 3, fB);

/*		cudaMemcpy(reQ1 + i*NR, d_ar0, sizeof(double)*NR, cudaMemcpyDeviceToHost);
		cudaMemcpy(imQ1 + i*NR, d_ai0, sizeof(double)*NR, cudaMemcpyDeviceToHost);
		for(int a = 0; a < NR; a++) {double r = 3.*wall/double(NR)*(a) - wall; fprintf(file, "%g,%g,%g, %g\n", z, r, reQ1[a+i*NR], imQ1[a+i*NR]);}*/

/*		if((enPrint != 0)&&(i%nd == 0))
		{
			cudaMemcpy(energy, d_W0, sizeof(double)*Np*Nq*Ns*Nv, cudaMemcpyDeviceToHost);
			cudaMemcpy(phase, d_Q0, sizeof(double)*Np*Nq*Ns*Nv, cudaMemcpyDeviceToHost);
			for(int en = 0; en < Np*Nq*Ns*Nv; en++)
			{
				fprintf(res_enr, "%g,%g\n", fmod(phase[en]+ 0.5*Pi, 2.*Pi) , energy[en]);

			}

		}*/

	}
		
	double rJ = 0, iJ = 0;	
	double rJ_1[512];
	double iJ_1[512];

	cudaMemcpy(rJ_1, (void*) d_int_rJ3_1, gridsize*sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(iJ_1, (void*) d_int_iJ3_1, gridsize*sizeof(double), cudaMemcpyDeviceToHost);

	for(int i = 0; i < gridsize; i++)
	{
		rJ += rJ_1[i]; iJ += iJ_1[i];
	//	printf("rJ[%i] = %g\tiJ[%i] = %g\n", i,rJ_1[i], i, iJ_1[i]);
	}

	
	//	printf("rJ = %g\tiJ = %g\n", rJ/gridsize*dz, iJ/gridsize*dz);


/*	delete [] energy;
	delete [] phase;
	fclose(res_enr);*/
/*
	fclose(file);
	delete[] reQ1;
	delete[] imQ1;*/



	return std::complex<double>(rJ/gridsize*dz, iJ/gridsize*dz);
	


}
void  Multiplier_SpCharge_2D::PrintCurrentsSpaceCharge2D(char *filename1, char *filename2)
{
	int GQ = Nq/NQ; int GS = Ns/NS; int GV = Nv;


	double A = A_stat;
	double rB = real(B_stat);
	double iB = imag(B_stat);
	
//	A = 6.077e-5; rB = -0.000102469; iB = 3.1098e-5;

//	int warpsize = NP*NQ*NS;
	int gridsize = GQ*GS*GV;

	ParamsM par = setPar();

	double K[4] = {0, 0.5, 0.5, 1.};

	dim3 motionwarp = dim3(NP, NQ, NS);
	dim3 motiongrid = dim3(GQ, GS, GV);

	cudaMemset(d_Qk, 0, 4*sizeof(double)*Np*Nq*Ns*Nv);	cudaMemset(d_W0, 0, sizeof(double)*Np*Nq*Ns*Nv);
	cudaMemset(d_Wk, 0, 4*sizeof(double)*Np*Nq*Ns*Nv);	cudaMemset(d_rAk, 0, sizeof(double)*4*NR);
	cudaMemset(d_Q0, 0, sizeof(double)*Np*Nq*Ns*Nv);	cudaMemset(d_iAk, 0, sizeof(double)*4*NR);
	cudaMemset(d_rAq1k, 0, sizeof(double)*4*NR);
	cudaMemset(d_iAq1k, 0, sizeof(double)*4*NR);
	

	cudaMemset(d_ai0, 0, sizeof(double)*NR);
	cudaMemset(d_ar0, 0, sizeof(double)*NR);
	cudaMemset(d_rAq1, 0, sizeof(double)*NR);	
	cudaMemset(d_iAq1, 0, sizeof(double)*NR);	

	cudaMemset(d_avEN, 0, sizeof(double));

	double La = Nperiods*period;
	double dz = Lmax/double(Nmax);

	int Nstop = ceil((La + Lb + Ld)/dz);
	int nd = round(period/dz);

	double z, fA, fB;

	cudaMemcpy(d_par, &par, sizeof(PAR), cudaMemcpyHostToDevice);

	double *reQ1 = new double [NR*Nstop];
	double *imQ1 = new double [NR*Nstop];	

	double *reQ3 = new double [NR*Nstop];
	double *imQ3 = new double [NR*Nstop];

	double *rJ3 = new double [GQ*GS*GV*Nmax];
	double *iJ3 = new double [GQ*GS*GV*Nmax];

	FILE *file1 = fopen(filename1,"w");
	FILE *file2 = fopen(filename2,"w");

	double rJ = 0, iJ = 0;	
//	double rJ_1[1024];
//	double iJ_1[1024];

	
	for(int i = 0; i < Nstop; i++)
	{
		z = double(i)*dz;
		fA = sin(Pi/La*z);
		fA *= A;
		if(z > La) fA = 0;
		if((z > La+Ld)&&(z < La+ Ld+ Lb)) fB = sin(Pi*(z - La - Ld)/Lb); else fB = 0;


		for(int k = 0; k < 4; k++)
		{
			motionstep_spacecharge_2d <<<motiongrid,motionwarp>>>(d_par, i, k, K[k], fA, rB*fB, iB*fB, 1);
			amplitudestep_spacecharge_2d  <<<1, Nq*Ns*Nv>>>(d_par, i, k, K[k], spchQ1, spchQ3);
		}
		spacecharge_2d_endstep<<<motiongrid, motionwarp>>>(d_par, i, 3, fB);

/*		cudaMemcpy(reQ1 + i*NR, d_rAq1, sizeof(double)*NR, cudaMemcpyDeviceToHost);
		cudaMemcpy(imQ1 + i*NR, d_iAq1, sizeof(double)*NR, cudaMemcpyDeviceToHost);

		cudaMemcpy(reQ3 + i*NR, d_ar0, sizeof(double)*NR, cudaMemcpyDeviceToHost);
		cudaMemcpy(imQ3 + i*NR, d_ai0, sizeof(double)*NR, cudaMemcpyDeviceToHost); */

	/*	for(int a = 0; a < NR; a++) 
		{
			double r = 3.*wall/double(NR)*(a) - wall; 
			fprintf(file1, "%g,%g,%g, %g\n", z, r, reQ1[a+i*NR], imQ1[a+i*NR]);
			fprintf(file2, "%g,%g,%g, %g\n", z, r, reQ3[a+i*NR], imQ3[a+i*NR]);
		}*/

		
	}

	printf("Loop over\n");
	printf("J3 size = %i\n", GQ*GS*GV*Nmax);

	printf("memcpy %i\n", cudaMemcpy(rJ3, (void*) d_rJ3, GQ*GS*GV*sizeof(double)*Nmax, cudaMemcpyDeviceToHost));
  	printf("memcpy %i\n", cudaMemcpy(iJ3, (void*) d_iJ3, GQ*GS*GV*sizeof(double)*Nmax, cudaMemcpyDeviceToHost));
	z = 0;
	for(int i = 0; i < Nstop; i++)
	{
		rJ = 0; iJ = 0;
		z = (double) i*dz;
		for(int a = 0; a < gridsize; a++)
		{
			if(Nmax*a + i >= GQ*GS*GV*Nmax) printf("i = %i, a = %i, adr = %i;", i, a, Nmax*a + i);
			rJ += rJ3[Nmax*a + i]; 
			iJ += iJ3[Nmax*a + i];
//			printf("rJ[%i] = %g\tiJ[%i] = %g\n", i,rJ_1[i], i, iJ_1[i]);
//			fprintf(file2, "%g,%i,%g,%g\n", i*dz, a, rJ3[Nmax*a + i]/*sin(Pi/La*a*dz)*/, iJ3[Nmax*a + i]/*sin(Pi/La*a*dz)*/);
		}
			
		fprintf(file2, "%g,%g,%g\n", i*dz, rJ/(double)gridsize/*sin(Pi/La*a*dz)*/, iJ/(double)gridsize/*sin(Pi/La*a*dz)*/);
	}	
	printf("\n");
		
	fclose(file1);
	fclose(file2);
	delete[] reQ1;
	delete[] imQ1;
	delete[] reQ3;
	delete[] imQ3;
	delete[] rJ3;
	delete[] iJ3;

	


}
double Multiplier_SpCharge_2D::ElectronsDeltaEnergy(double _A)
{
	return DeltaEnergySpaceCharge2D(_A);
}
cplx Multiplier_SpCharge_2D::ElectronCurrentA(double reA, double imA)
{
	return CurrentASpaceCharge2D(reA, imA);
}

cplx Multiplier_SpCharge_2D::HfElectronCurrent(double  _reB, double _imB, double _A)
{
	return CurrentBSpaceCharge2D(_reB, _imB, _A);
}
bool Multiplier_SpCharge_2D::initMultiplierSolver(int nz, double lsolver, double groupSpeedCoeff, char *_solverName)
{
	Device::initSolver(nz, lsolver, groupSpeedCoeff, _solverName);
//	int GQ = Nq / NQ; int GS = Ns / NS; int GV = Nv; //grid dimensions
	printf("Allocate space charge arrays ...\n<font color = \"red\">");
	gpuErrChk(cudaMalloc((void**)&d_Qk, 4 * sizeof(double)*Np*Nq*Ns*Nv))
	gpuErrChk(cudaMalloc((void**)&d_Wk, 4 * sizeof(double)*Np*Nq*Ns*Nv))
	gpuErrChk(cudaMalloc((void**)&d_Q0, sizeof(double)*Np*Nq*Ns*Nv))
	gpuErrChk(cudaMalloc((void**)&d_W0, sizeof(double)*Np*Nq*Ns*Nv))
	gpuErrChk(cudaMalloc((void**)&d_ifnotdestroyed, sizeof(int)*Np*Nq*Ns*Nv))


	gpuErrChk(cudaMalloc((void**)&d_ar0, sizeof(double)*NR))
	gpuErrChk(cudaMalloc((void**)&d_ai0, sizeof(double)*NR))

	gpuErrChk(cudaMalloc((void**)&d_int_rJ3, sizeof(double)*Nq*Ns*Nv))		//Третья гармоника пространственного заряда
	gpuErrChk(cudaMalloc((void**)&d_int_iJ3, sizeof(double)*Nq*Ns*Nv))		//(усреднённая только по фазе)

	gpuErrChk(cudaMalloc((void**)&d_rAq1, sizeof(double)*NR))
	gpuErrChk(cudaMalloc((void**)&d_iAq1, sizeof(double)*NR))

	gpuErrChk(cudaMalloc((void**)&d_int_rQ1, sizeof(double)*Nq*Ns*Nv))		//Первая гармоника пространственного заряда
	gpuErrChk(cudaMalloc((void**)&d_int_iQ1, sizeof(double)*Nq*Ns*Nv))		//(усреднённая только по фазе)

	gpuErrChk(cudaMalloc((void**)&d_rAk, 4 * NR*sizeof(double)))	 //Функция поля пространственного заряда
	gpuErrChk(cudaMalloc((void**)&d_iAk, 4 * NR*sizeof(double)))
	gpuErrChk(cudaMalloc((void**)&d_rAq1k, 4 * NR*sizeof(double)))
	gpuErrChk(cudaMalloc((void**)&d_iAq1k, 4 * NR*sizeof(double)))
	gpuErrChk(cudaMalloc((void**)&d_radii, Nq*Ns*Nv*sizeof(double)))

	printf("<\font>... end allocation\n");


	return 1;
}
