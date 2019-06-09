#include <stdio.h>
#include <complex>
#include "multiplier.h"
#include "multiplier_three_cavity.h"
#include "multiplier_multimode.h"
#include "cu_mult.h"
//#include "print.h"
  

#define PAR ParamsM
#ifndef dm_Pi
__device__ const double dm_Pi = 3.141592653589793;
#endif

int Nz = 0;
double Lsolver = 0;
double MultiplierGroupSpeedCoefficient = 0;

__device__  double grSpeedCoeff;
int *d_Nz;
double *d_Lsolver;

//cudaPitchedPtr  d2_rJ3, d2_iJ3, d2_int_rJ3, d2_int_iJ3, d2_W;
//cudaPitchedPtr  d1_rJ3, d1_iJ3, d1_int_rJ3, d1_int_iJ3, d1_W;

__device__ void biReduce(double *A, double *B, int p0, int datasize, int logsize)
{
	int stride = datasize;
	for (int q = 1; q <= logsize; q++)
	{
		stride = stride >> 1;
		if (p0 < stride)
		{
			A[p0] += A[p0 + stride];

		}
		else
		{
			if (p0 < 2 * stride)
			{

				B[p0 - stride] += B[p0];

			}
		}
		__syncthreads();
	}
}
__device__ double dh( double delta)
{
	return delta*delta*(grSpeedCoeff);
}
__device__ void funcA (double z, double *rA, double *iA, double2 *Amps, int Na)
{
	*rA = 0; 
	*iA = 0;
	double rF, iF;	
	for(int i = 0; i < Na; i++)
	{
		sincos(z*double(i - Na/2), &iF, &rF);
		*rA += Amps[i].x*rF - Amps[i].y*iF;
		*iA += Amps[i].x*iF + Amps[i].y*rF;
	}

}
__global__ void 
__launch_bounds__ (512, 2)
MotionEquationMultiplier(PAR *par, double Lstop, int Nharm, double A, double2 B)//Fixed Structure
{

	unsigned int p0 = threadIdx.x;			unsigned int Np = blockDim.x;
	unsigned int q0 = threadIdx.y;			unsigned int Nq = blockDim.y;
	unsigned int s0 = threadIdx.z;			unsigned int Ns = blockDim.z;
    
	unsigned int q_init =		 Nq*blockIdx.x + q0;		unsigned int Nq_max =		 Nq*gridDim.x;
	unsigned int s_init =		 Ns*blockIdx.y + s0;		unsigned int Ns_max =		 Ns*gridDim.y;
	unsigned int v_init =			blockIdx.z;				unsigned int Nv_max =		    gridDim.z;

	int warpsize = Np*Nq*Ns;
	int log2warpsize = round(log2((double)warpsize));
	int X = blockIdx.x + gridDim.x*blockIdx.y + gridDim.y*gridDim.x*blockIdx.z;

	double la, lb, ld, h, k1, voltage, g1, g3;

	__shared__ double avEN, int_rJ3, int_iJ3;
		
	int N; 

	double dz;

	N = par->Nz;

	la = par->la; lb = par->lb; ld = par->ld;
	h  = par->h; k1 = par->k1;
	g1 = par->g1; g3 = par->g3;
	voltage = par->voltage;

	
	double *rJ3 = par->rJ3; 
	double *iJ3 = par->iJ3;

	dz = par->L/(double)N;


	double z;
	int ifinal = floor(Lstop/dz);

	
	double Q, Qk1, Qk2, Qk3, Qk4;
	double W, Wk1, Wk2, Wk3, Wk4;
	double fA, fB, rA, r;
	double Wmax, Wmin;
	double ifnotdestroyed = 1;

	double R_cyclotron, R_center, kappa_cyclotron, phase_cyclotron, initial_angle, angle_spread_factor;
	double wall = par->wall;

	R_center =  0.5*wall + wall*((double)q_init-0.5*(double)Nq_max)/(double)(Nq_max);
	initial_angle = 0;//(0.0810194 - 2.05972*R_center +  28.0433*R_center*R_center);
	R_cyclotron =	(0.568*initial_angle + 0*0.035156*((double)v_init)/double(Nv_max));
	kappa_cyclotron = 1.758;
	phase_cyclotron = 2.*dm_Pi*(double)s_init/(double)Ns_max;

	/*R_center =  0.5*wall + wall*((double)q_init-0.5*(double)Nq_max)/(double)(Nq_max);
	initial_angle = 0.087*((double)v_init - 0.5*Nv_max)/double(Nv_max);
	R_cyclotron =	0.744*initial_angle;//для f = 86.6 коэфф. 0,568; для f = 95.5: 0.744
	kappa_cyclotron = 1.344; //для f = 86.6 коэфф. 1.758; для f = 95.5: 1.344
	phase_cyclotron = 2.*dm_Pi*(double)s_init/(double)Ns_max;*/

	double en0 = 1. + voltage/511.;	
	angle_spread_factor = 1./sqrt(1. + initial_angle*initial_angle);
																		// Вызывает недопустимую операцию при достаточно больших wall. Наверное из-за резкости параболы в initial_angle
	Q = 2.*dm_Pi/double(Np)*double(p0);
	W = 0; 

	__shared__ double sh_sinQ[NS*NQ*NP];
	__shared__ double sh_cosQ[NS*NQ*NP];

/*	__shared__ double    shQ[NS][NQ][NP];
	__shared__ double    shW[NS][NQ][NP];
	__shared__ double d2_rJ3[NQ][NP];
	__shared__ double d2_iJ3[NQ][NP];
	__shared__ double d1_rJ3[NP];
	__shared__ double d1_iJ3[NP];*/

	


	double PH, EN, cosPH, sinPH, cosPS, sinPS, rB, iB;

	
	double H = h;//+dh(delta);

	if(p0+q_init+s_init + v_init  == 0) 
	{
		rJ3[0]  = 0;
		iJ3[0]  = 0;
	}

	if(p0+q0+s0== 0)
	{
		par->int_rJ3[gridDim.x*gridDim.y*blockIdx.z + gridDim.x*blockIdx.y+blockIdx.x] = 0;
		par->int_iJ3[gridDim.x*gridDim.y*blockIdx.z + gridDim.x*blockIdx.y+blockIdx.x] = 0;
		int_rJ3 = 0;
		int_iJ3 = 0;

		avEN    = 0;
	}

	
//	if(s0+p0+q0 == 0) printf("la = %g, ld = %g, lb = %g \n", la, ld, lb);
	
	int i = 0;
	for(i = 1; i < N; i++)
	{


/////////////////
		z = (double)i*dz; 
		r = (R_center + R_cyclotron*cos(kappa_cyclotron*z + phase_cyclotron));
		 

		ifnotdestroyed *= (r > -wall)? 1. : 0.;										///!!!!!!
	//	ifnotdestroyed = 1;

		PH = Q;
		EN = W + en0;

//		if((s0+p0+q0 == 0)) printf("%g\t%g\n", z, fB);

 
		fA = ((z<la)?sin(dm_Pi/la*z)*                             exp(-g1*r):0);
		fB = (((z>la+ld)&&(z < la+ld+lb))?sin(dm_Pi/lb*(z-la-ld))*exp(-g3*r):0);
		rA = A*fA; 
		rB = B.x*fB;
		iB = B.y*fB;
		
		sincos(PH, &sinPH, &cosPH);
		sincos(3.*PH, &sinPS, &cosPS);

		Qk1 = dz*(H - k1*EN*angle_spread_factor/sqrt(EN*EN-1.));
		Wk1 = -dz*((rA*cosPH)+(rB*cosPS-iB*sinPS))*ifnotdestroyed;

/////////////////
		z = ((double)i+0.5)*dz;
		r = (R_center + R_cyclotron*cos(kappa_cyclotron*z + phase_cyclotron));
	//	ifnotdestroyed *= (r > -wall)? 1. : 0.;
	//	ifnotdestroyed = 1;

		PH = Q + 0.5*Qk1; 
		EN = W + 0.5*Wk1 + en0;

		fA = ((z<la)?sin(dm_Pi/la*z) *                           exp(-g1*r):0);
		fB = ((z>la+ld)&&(z < la+ld+lb))?sin(dm_Pi/lb*(z-la-ld))*exp(-g3*r):0; 
		rA = A*fA;
		rB = B.x*fB;
		iB = B.y*fB;
		
		sincos(PH, &sinPH, &cosPH); 
		sincos(3.*PH, &sinPS, &cosPS);

		Qk2 = dz*(H - k1*EN*angle_spread_factor/sqrt(EN*EN-1.));
		Wk2 = -dz*((rA*cosPH)+(rB*cosPS-iB*sinPS))*ifnotdestroyed;	
/////////////////
		PH = Q + 0.5*Qk2; 
		EN = W + 0.5*Wk2 + en0;

		sincos(PH, &sinPH, &cosPH);	
		sincos(3.*PH, &sinPS, &cosPS);

		Qk3 = dz*(H - k1*EN*angle_spread_factor/sqrt(EN*EN-1.));
		Wk3 = -dz*((rA*cosPH)+(rB*cosPS-iB*sinPS))*ifnotdestroyed;				
/////////////////	
		z = ((double)i+1.)*dz;
		r = (R_center + R_cyclotron*cos(kappa_cyclotron*z + phase_cyclotron));
//		ifnotdestroyed *= (r > -wall)? 1. : 0.;
//		ifnotdestroyed = 1;

		PH = Q + Qk3; 
		EN = W + Wk3 + en0;

		fA = ((z<la)?					  sin(dm_Pi/la*z)*		  exp(-g1*r):0);
		fB = (((z>la+ld)&&(z < la+ld+lb))?sin(dm_Pi/lb*(z-la-ld))*exp(-g3*r):0); 
		rA = A*fA;
		rB = B.x*fB;
		iB = B.y*fB;

		
		sincos(PH, &sinPH, &cosPH);
		sincos(3.*PH, &sinPS, &cosPS);

		Qk4= dz*(H - k1*EN*angle_spread_factor/sqrt(EN*EN-1.));
		Wk4= -dz*((rA*cosPH)+(rB*cosPS-iB*sinPS))*ifnotdestroyed;	
///////////////		
			
		Q +=  1./6.*(Qk1+2.*Qk2+2.*Qk3+Qk4);

		W +=  1./6.*(Wk1+2.*Wk2+2.*Wk3+Wk4);

	/*	shQ[s0][q0][p0] = Q;
		shW[s0][q0][p0] = W;*/

		__syncthreads();
		

		sincos(double(Nharm)*Q, &sinPH, &cosPH);

		if(Nharm == 1) 
			fB = ((z<la)?sin(dm_Pi/la*z):0)*exp(-g1*r);			//fB используется как множитель при интегрировании тока вдоль продольной координаты ВНИМАНИЕ!! fB зависит от q0, s0
		else 
			fB = (((z>la+ld)&&(z < la+ld+lb))?sin(dm_Pi/lb*(z-la-ld)):0)*exp(-g3*r); 

		fB *= ifnotdestroyed;

		int xi = p0 + Np*q0 + Np*Nq*s0;



		sh_sinQ[xi] = sinPH*fB;
		sh_cosQ[xi] = cosPH*fB;
		
		unsigned int stride = warpsize;

		__syncthreads();
	
		for(int q = 1; q <= log2warpsize; q++)
		{
			stride = stride >> 1;//roundf(powf(2., q));
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
		
//		if((i == 1300)) printf("%g\n", Q);

		if(xi == 0)
		{
			rJ3[X*N+i] =  sh_cosQ[0];	
			iJ3[X*N+i] =  -sh_sinQ[0];

//			if((i == 1300)) printf("\n%g\n", sh_cosQ[0]);

			int_rJ3 += sh_cosQ[0];
			int_iJ3 += -sh_sinQ[0];
		}


	/*	 //////// усреднение Nharm гармоники
		if(s0 == 0)
		{
			 double tmp_rJ3 = 0, tmp_iJ3 = 0, tmpPhCycl = 0;
			 

			for(int ii = 0; ii < Ns; ii++)
			{

				int ii_init =  Ns*blockIdx.y + ii;
				tmpPhCycl = 2.*dm_Pi*(double)ii_init/(double)Ns_max;
				r = (R_center + R_cyclotron*cos(kappa_cyclotron*z + tmpPhCycl));

								
				if(Nharm == 1) 
					fB = ((z<la)?sin(dm_Pi/la*z):0)*exp(-g1*r);			//fB используется как множитель при интегрировании тока вдоль продольной координаты ВНИМАНИЕ!! fB зависит от q0, s0
				else 
					fB = (((z>la+ld)&&(z < la+ld+lb))?sin(dm_Pi/lb*(z-la-ld)):0)*exp(-g3*r); 

				fB *= ifnotdestroyed;

				PH = shQ[ii][q0][p0];
				sincos((double)Nharm*PH, &sinPS, &cosPS);

				tmp_rJ3 += cosPS*fB;		
				tmp_iJ3 -= sinPS*fB;

				
				
			}
			d2_rJ3[q0][p0] = tmp_rJ3;
			d2_iJ3[q0][p0] = tmp_iJ3;

		
		}


		__threadfence();
		__syncthreads();
		if(s0 + q0 == 0)
		{
			double tmp_rJ3 = 0, tmp_iJ3 = 0;

			for(int ii = 0; ii < Nq; ii++)
			{
				tmp_rJ3 += d2_rJ3[ii][p0];
				tmp_iJ3 += d2_iJ3[ii][p0];

			}
			d1_rJ3[p0] = tmp_rJ3;
			d1_iJ3[p0] = tmp_iJ3;
		
		}
		__threadfence();
		__syncthreads();
		if(p0 + q0 +s0 == 0)
		{
			double tmp_rJ3 = 0, tmp_iJ3 = 0;

			for(int ii = 0; ii < Np; ii++)
			{
				tmp_rJ3 += d1_rJ3[ii];
				tmp_iJ3 += d1_iJ3[ii];
			}
			rJ3[i] = tmp_rJ3;
			iJ3[i] = tmp_iJ3;

			int_rJ3 += tmp_rJ3;
			int_iJ3 += tmp_iJ3;
//			if(i == ifinal) printf("<< %g  %g\n", A*int_rJ3, A*tmp_rJ3);
		}
		__threadfence();
		__syncthreads();
		//////////////////// конец усреднения Nharm гармоники
*/
		if(i == ifinal) 
		{	
			sh_sinQ[xi] = W;
			__syncthreads();
			stride = warpsize;
			for(int q = 1; q <= log2warpsize; q++)
			{
				stride = stride >> 1;//warpsize/roundf(powf(2., q));
				if(xi < stride)
				{
					sh_sinQ[xi] += sh_sinQ[xi + stride];
				}
				__syncthreads();
			}

			if(xi == 0)
			{	
				avEN = sh_sinQ[0];
			}
			__syncthreads();

			sh_sinQ[xi] = W;
			sh_cosQ[xi] = W;
			stride = warpsize;
			for(int q = 1; q <= log2warpsize; q++)
			{
				stride = stride >> 1;// stride = warpsize/roundf(powf(2., q));
				if(xi < stride)
				{
					sh_sinQ[xi] = (sh_sinQ[xi] > sh_sinQ[xi + stride]) ? sh_sinQ[xi] : sh_sinQ[xi + stride];
				}
				else
				{
					if(xi < 2*stride)
					{
						sh_cosQ[xi - stride] = (sh_cosQ[xi - stride] < sh_cosQ[xi]) ? sh_cosQ[xi - stride] :  sh_cosQ[xi];
					}
				}
				__syncthreads();
			}

			if(xi == 0)
			{	
				Wmax = sh_sinQ[0];
				Wmin = sh_cosQ[0];
			}

/*
			if(s0 == 0)
			{
				double tmp_W = 0;

				for(int ii = 0; ii < Ns; ii++)
				{
					EN = shW[ii][q0][p0];
					tmp_W += EN;
				}
				d2_rJ3[q0][p0] = tmp_W;
			//	if((p0 == 0)) printf(" %g >>, \n", d2_rJ3[q0][p0]);
				
			}
			__threadfence();
			__syncthreads();
			if(s0 + q0 == 0)
			{
	
				double tmp_rJ3 = 0;
				for(int ii = 0; ii < Nq; ii++)
					tmp_rJ3 += d2_rJ3[ii][p0];
				d1_rJ3[p0] = tmp_rJ3;

			}
			__threadfence();
			__syncthreads();
			if(p0 + q0 +s0 == 0)
			{
				double tmp_rJ3 = 0;
				for(int ii = 0; ii < Np; ii++)
					tmp_rJ3 += d1_rJ3[ii];
				(avEN) += tmp_rJ3;
			}
*/
			__syncthreads();
			
			
		}
		__threadfence();
		__syncthreads();




		if(i > ifinal) break;
	

	}
//	printf("END\t");

	
//	if(p0 + s0 == 0) printf("(%i, %i, %i)...<%g, %g> =?= <%g>...\n", blockIdx.x, blockIdx.y, blockIdx.z,A*int_rJ3*dz, A*int_iJ3*dz, avEN);

	__syncthreads();
//	if(p0+q_init+s_init + v_init  == 0)
	if(p0+q0+s0 == 0)
	{
	/*	printf("%i, %i, %i\t (%g, %g)\n",
			blockIdx.x, blockIdx.y, blockIdx.z,
			     int_rJ3/double(Np*Nq_max*Ns_max*N)*(par->L), 
				 int_iJ3/double(Np*Nq_max*Ns_max*N)*(par->L)) ;*/
	    par->avEN[gridDim.x*gridDim.y*blockIdx.z + gridDim.x*blockIdx.y+blockIdx.x] = avEN;
		par->int_rJ3[gridDim.x*gridDim.y*blockIdx.z + gridDim.x*blockIdx.y+blockIdx.x] = int_rJ3;
		par->int_iJ3[gridDim.x*gridDim.y*blockIdx.z + gridDim.x*blockIdx.y+blockIdx.x] = int_iJ3;
		par->Wmax[gridDim.x*gridDim.y*blockIdx.z + gridDim.x*blockIdx.y+blockIdx.x] = Wmax;
		par->Wmin[gridDim.x*gridDim.y*blockIdx.z + gridDim.x*blockIdx.y+blockIdx.x] = Wmin;

	}

}

__global__ void 
__launch_bounds__ (512, 2)
MotionEquationMultiplierDoubleScheme(PAR par, double Lstop, int Nharm, double A, double2 A2, double2 B)//Fixed Structure
{

	unsigned int p0 = threadIdx.x;
	unsigned int q0 = threadIdx.y;
	unsigned int s0 = threadIdx.z;

    unsigned int Np = blockDim.x;
	unsigned int Nq = blockDim.y;
	unsigned int Ns = blockDim.z;

	unsigned int q_init =		 Nq*blockIdx.x + q0;
	unsigned int s_init =		 Ns*blockIdx.y + s0;
	unsigned int v_init =			blockIdx.z;

	unsigned int Nq_max =		 Nq*gridDim.x;
	unsigned int Ns_max =		 Ns*gridDim.y;
	unsigned int Nv_max =		    gridDim.z;

//	printf("Thread %i/%i, %i/%i started;    ", q0, Nq, p0, Np);

	double la1, la2, lb, ld1, ld2, h, k1, voltage, g1, g3;
	double fA2, rA2, iA2;

	__shared__ double avEN, int_rJ3, int_iJ3;
		
//	printf("Step0;    ");
	int N; 

	double dz;

	N = par.Nz;


	la1 = par.la1; la2 = par.la2;  lb = par.lb; 
	ld1 = par.ld1; ld2 = par.ld2;
	h  = par.h; k1 = par.k1;
	g1 = par.g1; g3 = par.g3;
	voltage = par.voltage;

	double la_tot = la1 + la2 + ld1;
	double ifnotdestroyed = 1;
	
	double *rJ3 = par.rJ3; 
	double *iJ3 = par.iJ3;

//	double *int_rJ3 = par.int_rJ3;
//	double *int_iJ3 = par.int_iJ3;


	dz = par.L/(double)N;

	double z;
	int ifinal = floor(Lstop/dz);
	
	double Q, Qk1, Qk2, Qk3, Qk4;
	double W, Wk1, Wk2, Wk3, Wk4;
	double fB;

	double R_cyclotron, R_center, kappa_cyclotron, phase_cyclotron, initial_angle;
	double wall = par.wall, r;
	
	R_center =  0.5*wall + wall*((double)q_init-0.5*(double)Nq_max)/(double)(Nq_max);
	initial_angle = (0.0810194 - 2.05972*R_center +  28.0433*R_center*R_center);
	R_cyclotron =	0.568*initial_angle + 0.035156*((double)v_init)/double(Nv_max);
	kappa_cyclotron = 1.758;
	phase_cyclotron = 2.*dm_Pi*(double)s_init/(double)Ns_max;

	double en0 = 1. + voltage/511.;	
	en0 -= 0.5*initial_angle*initial_angle*(en0*en0 - 1)*en0;

	
	Q = 2.*dm_Pi/double(Np)*double(p0);// + 1./(double)Nq*((double)q0 + (double)s0/(double)Ns));
	W = 0; 

	__shared__ double    shQ[NS][NQ][NP];
	__shared__ double    shW[NS][NQ][NP];
	__shared__ double d2_rJ3[NQ][NP];
	__shared__ double d2_iJ3[NQ][NP];
	__shared__ double d1_rJ3[NP];
	__shared__ double d1_iJ3[NP];

	


	double PH, EN, cosPH, sinPH, cosPS, sinPS, rA, rB, iB;

	// printf("Step3;    ");	

	
	double H = h;//+dh(delta);

	if(p0+q_init+s_init + v_init  == 0) 
	{
		rJ3[0]  = 0;
		iJ3[0]  = 0;
	//	printf("init:\t%i\t%i\t%i\t%i...........\n",p0, q_init,s_init,v_init);

	}

	if(p0+q0+s0== 0)
	{
		par.int_rJ3[gridDim.x*gridDim.y*blockIdx.z + gridDim.x*blockIdx.y+blockIdx.x] = 0;
		par.int_iJ3[gridDim.x*gridDim.y*blockIdx.z + gridDim.x*blockIdx.y+blockIdx.x] = 0;
		int_rJ3 = 0;
		int_iJ3 = 0;

		avEN    = 0;
	}

	int i = 0;
	for(i = 1; i < N; i++)
	{


/////////////////
		z = (double)i*dz; 
		r = (R_center + R_cyclotron*cos(kappa_cyclotron*z + phase_cyclotron));
		ifnotdestroyed *= (r > -wall)? 1. : 0.;

		PH = Q;
		EN = W + en0;

		rA = A  *((z<la1)				? sin(dm_Pi/la1*z)							*exp(-g1*r):0); 
		fA2= ( ((la1+ld1<z)&&(z<la_tot))? sin(dm_Pi/la2*(z-la1 - ld1))				*exp(-g1*r):0); 
		rA2 = A2.x*fA2;
		iA2 = A2.y*fA2;

		fB = ( ((z> la_tot+ld2)&&(z < la_tot+ld2+lb))?sin(dm_Pi/lb*(z-la_tot-ld2))  *exp(-g3*r):0);
		rB = B.x*fB;
		iB = B.y*fB;
		
		sincos(PH, &sinPH, &cosPH);
		sincos(3.*PH, &sinPS, &cosPS);

		Qk1 = dz*(H - k1*EN/sqrt(EN*EN-1.));
		Wk1 = -dz*((rA*cosPH)+(rA2*cosPH-iA2*sinPH)+(rB*cosPS-iB*sinPS))*ifnotdestroyed;

/////////////////
		z = ((double)i+0.5)*dz;
		r = (R_center + R_cyclotron*cos(kappa_cyclotron*z + phase_cyclotron));
		ifnotdestroyed *= (r > -wall)? 1. : 0.;

		PH = Q + 0.5*Qk1; 
		EN = W + 0.5*Wk1 + en0;

		rA = A  *((z<la1)				? sin(dm_Pi/la1*z)							*exp(-g1*r):0); 
		fA2= ( ((la1+ld1<z)&&(z<la_tot))? sin(dm_Pi/la2*(z-la1 - ld1))				*exp(-g1*r):0); 
		rA2 = A2.x*fA2;
		iA2 = A2.y*fA2;

		fB = ( ((z> la_tot+ld2)&&(z < la_tot+ld2+lb))?sin(dm_Pi/lb*(z-la_tot-ld2))   *exp(-g3*r):0);
		rB = B.x*fB;
		iB = B.y*fB;
		
		sincos(PH, &sinPH, &cosPH); 
		sincos(3.*PH, &sinPS, &cosPS);

		Qk2 = dz*(H - k1*EN/sqrt(EN*EN-1.));
		Wk2 = -dz*((rA*cosPH)+(rA2*cosPH-iA2*sinPH)+(rB*cosPS-iB*sinPS))*ifnotdestroyed;
/////////////////
		PH = Q + 0.5*Qk2; 
		EN = W + 0.5*Wk2 + en0;

		sincos(PH, &sinPH, &cosPH);	
		sincos(3.*PH, &sinPS, &cosPS);

		Qk3 = dz*(H - k1*EN/sqrt(EN*EN-1.));
		Wk3 = -dz*((rA*cosPH)+(rA2*cosPH-iA2*sinPH)+(rB*cosPS-iB*sinPS))*ifnotdestroyed;
/////////////////	
		z = ((double)i+1)*dz;
		r = (R_center + R_cyclotron*cos(kappa_cyclotron*z + phase_cyclotron));
		ifnotdestroyed *= (r > -wall)? 1. : 0.;
		PH = Q + Qk3; 
		EN = W + Wk3 + en0;

		rA = A  *((z<la1)				? sin(dm_Pi/la1*z)							*exp(-g1*r):0); 
		fA2= ( ((la1+ld1<z)&&(z<la_tot))? sin(dm_Pi/la2*(z-la1 - ld1))				*exp(-g1*r):0); 
		rA2 = A2.x*fA2;
		iA2 = A2.y*fA2;

		fB = ( ((z> la_tot+ld2)&&(z < la_tot+ld2+lb))?sin(dm_Pi/lb*(z-la_tot-ld2))   *exp(-g3*r):0);
		rB = B.x*fB;
		iB = B.y*fB;

		
		sincos(PH, &sinPH, &cosPH);
		sincos(3.*PH, &sinPS, &cosPS);

		Qk4= dz*(H - k1*EN/sqrt(EN*EN-1.));
		Wk4= -dz*((rA*cosPH)+(rA2*cosPH-iA2*sinPH)+(rB*cosPS-iB*sinPS))*ifnotdestroyed;
///////////////		
			
		Q +=  1./6.*(Qk1+2.*Qk2+2.*Qk3+Qk4);

		W +=  1./6.*(Wk1+2.*Wk2+2.*Wk3+Wk4);

		
//		printf("#< %i -> (%i, %i, %i)>\t", i, s0, q0, p0);
//		printf("#<%i>", Np_Q*Nq*(Ns*i+s0) + Np_Q*q0 + p0);

//		if(q0+p0+s0 == 0) printf("%i", i);

		shQ[s0][q0][p0] = Q;
				   	   
		shW[s0][q0][p0] = W;

			
		__threadfence();
		__syncthreads();

		//////// усреднение какой-то гармоники
		if(s0 == 0)
		{
			 double tmp_rJ3 = 0, tmp_iJ3 = 0;

			for(int ii = 0; ii < Ns; ii++)
			{
				PH = shQ[ii][q0][p0];
				sincos((double)Nharm*PH, &sinPS, &cosPS);

				tmp_rJ3 += cosPS;
				tmp_iJ3 -= sinPS;
			}
			d2_rJ3[q0][p0] = tmp_rJ3;
			d2_iJ3[q0][p0] = tmp_iJ3;

			
		}

		

		__threadfence();
		__syncthreads();
		if(s0 + q0 == 0)
		{
			double tmp_rJ3 = 0, tmp_iJ3 = 0;

			for(int ii = 0; ii < Nq; ii++)
			{
				tmp_rJ3 += d2_rJ3[ii][p0];
				tmp_iJ3 += d2_iJ3[ii][p0];
			}
			d1_rJ3[p0] = tmp_rJ3;
			d1_iJ3[p0] = tmp_iJ3;
		}
		__threadfence();
		__syncthreads();
		if(p0 + q0 +s0 == 0)
		{
			double tmp_rJ3 = 0, tmp_iJ3 = 0;

			for(int ii = 0; ii < Np; ii++)
			{
				tmp_rJ3 += d1_rJ3[ii];
				tmp_iJ3 += d1_iJ3[ii];
			}
			rJ3[i] = tmp_rJ3;
			iJ3[i] = tmp_iJ3;
			int_rJ3 += tmp_rJ3*((Nharm == 3)?fB:fA2);
			int_iJ3 += tmp_iJ3*((Nharm == 3)?fB:fA2);
		}
		__threadfence();
		__syncthreads();
		//////////////////// конец усреднения какой-то гармоники


//		if((q0+p0 == 0)&&(s0 == 0)) printf("%i\t%g\t%g\n", q0, PH, EN);

//		if(q0+p0+s0 == 0) printf("....%i\t", i);

		/////////////////////// усреднение энергии
		if(i ==  ifinal) 
		{
			if(s0 == 0)
			{
				double tmp_W = 0;

				for(int ii = 0; ii < Ns; ii++)
				{
					EN = shW[ii][q0][p0];
					tmp_W += EN;

				}
				d2_rJ3[q0][p0] = W;
			}
			__threadfence();
			__syncthreads();
			if(s0 + q0 == 0)
			{
	
				double tmp_rJ3 = 0;
				for(int ii = 0; ii < Nq; ii++)
					tmp_rJ3 += d2_rJ3[ii][p0];
				d1_rJ3[p0] = tmp_rJ3;

			}
			__threadfence();
			__syncthreads();
			if(p0 + q0 +s0 == 0)
			{
				double tmp_rJ3 = 0;
				for(int ii = 0; ii < Np; ii++)
					tmp_rJ3 += d1_rJ3[ii];
				(avEN) += tmp_rJ3;

			}
		} 
		///////////////// конец усреднения энергии
		__threadfence();
		__syncthreads();






		if(i > ifinal) break;
	

	}

	__syncthreads();
	if(p0+q0+s0 == 0)
	{
		*par.avEN = avEN;
		par.int_rJ3[gridDim.x*gridDim.y*blockIdx.z + gridDim.x*blockIdx.y+blockIdx.x] = int_rJ3;
		par.int_iJ3[gridDim.x*gridDim.y*blockIdx.z + gridDim.x*blockIdx.y+blockIdx.x] = int_iJ3;

	}

}

__global__ void 
__launch_bounds__ (512, 2)
MotionEquationMultiplierMultiModes(PAR par, double Lstop, int Nharm, int Na, double2 B)//Fixed Structure
{

	unsigned int p0 = threadIdx.x;
	unsigned int q0 = threadIdx.y;
	unsigned int s0 = threadIdx.z;

    unsigned int Np = blockDim.x;
	unsigned int Nq = blockDim.y;
	unsigned int Ns = blockDim.z;

	unsigned int q_init =		 Nq*blockIdx.x + q0;
	unsigned int s_init =		 Ns*blockIdx.y + s0;
	unsigned int v_init =			blockIdx.z;

	unsigned int Nq_max =		 Nq*gridDim.x;
	unsigned int Ns_max =		 Ns*gridDim.y;
//	unsigned int Nv_max =		    gridDim.z;

	int warpsize = Np*Nq*Ns;
	int log2warpsize = round(log2((double)warpsize));

	double la,  lb, ld,  h, k1, voltage, g1, g3;
	double rA1, iA1;

	__shared__ double avEN, int_rJ3, int_iJ3, int_rJ3_1, int_iJ3_1;

	int N; 

	double dz;

	N = par.Nz;

	la = par.la1;  lb = par.lb; 
	ld = par.ld; 
	h  = par.h; k1 = par.k1;
	g1 = par.g1; g3 = par.g3;
	voltage = par.voltage;

	double ifnotdestroyed = 1;
	
	double *rJ3 = par.rJ3; 
	double *iJ3 = par.iJ3;
	double2 *Amps = (double2 *)par.Amps;

//	double *int_rJ3 = par.int_rJ3;
//	double *int_iJ3 = par.int_iJ3;

	dz = par.L/(double)N;

	double z;
	int ifinal = floor(Lstop/dz);
	
	double Q, Qk1, Qk2, Qk3, Qk4;
	double W, Wk1, Wk2, Wk3, Wk4;
	double fB;

	double R_cyclotron, R_center, kappa_cyclotron, phase_cyclotron, initial_angle;
	double wall = par.wall, r;
	
	R_center =  0.5*wall + wall*((double)q_init-0.5*(double)Nq_max)/(double)(Nq_max);
	initial_angle = 0;//(0.0810194 - 2.05972*R_center +  28.0433*R_center*R_center);
	R_cyclotron =	0;//0.568*initial_angle + 0.035156*((double)v_init)/double(Nv_max);
	kappa_cyclotron = 1.758;
	phase_cyclotron = 2.*dm_Pi*(double)s_init/(double)Ns_max;

	double en0 = 1. + voltage/511.;	
	en0 -= 0.5*initial_angle*initial_angle*(en0*en0 - 1)*en0;
	double beta0 = sqrt(en0*en0 - 1)/en0;
//	double Delta =  k1*dm_Pi/(la*beta0) ;// \delta f / f  = (k_0 \pi /L)/beta_ph


	Q = 2.*dm_Pi/double(Np)*double(p0);// + 1./(double)Nq*((double)q0 + (double)s0/(double)Ns));
	W = 0; 

	__shared__ double2 shAmps[NP];
	__shared__ double sh_sinQ[NS*NQ*NP];
	__shared__ double sh_cosQ[NS*NQ*NP];

	double PH, EN, cosPH, sinPH, cosPS, sinPS, rB, iB;

	double H = h;//+dh(delta);

	if(p0+q_init+s_init + v_init  == 0) 
	{
		rJ3[0]  = 0;
		iJ3[0]  = 0;
	}

	if(p0+q0+s0== 0)
	{
		par.int_rJ3[gridDim.x*gridDim.y*blockIdx.z + gridDim.x*blockIdx.y+blockIdx.x] = 0;
		par.int_iJ3[gridDim.x*gridDim.y*blockIdx.z + gridDim.x*blockIdx.y+blockIdx.x] = 0;
		par.int_rJ3_1[gridDim.x*gridDim.y*blockIdx.z + gridDim.x*blockIdx.y+blockIdx.x] = 0;
		par.int_iJ3_1[gridDim.x*gridDim.y*blockIdx.z + gridDim.x*blockIdx.y+blockIdx.x] = 0;
		int_rJ3 = 0;
		int_iJ3 = 0;
		int_rJ3_1 = 0;
		int_iJ3_1 = 0;

		avEN    = 0;
		
	}



	if((q0 + s0 == 0)&&(p0 < Na))
	{
		shAmps[p0] = Amps[p0];
	}


	
	__syncthreads();
	int i = 0;
	for(i = 1; i < N; i++)
	{


/////////////////
		z = (double)i*dz; 
		r = (R_center + R_cyclotron*cos(kappa_cyclotron*z + phase_cyclotron));
		ifnotdestroyed *= 1;//(r > -wall)? 1. : 0.;

		PH = Q;
		EN = W + en0;

		funcA(dm_Pi/la*z, &rA1, &iA1, shAmps, Na);
	

		if(z > la) {rA1 =0; iA1 = 0;}

		rA1 *= exp(-g1*r); iA1 *= exp(-g1*r);

		fB = ( ((z> la+ld)&&(z < la+ld+lb))?sin(dm_Pi/lb*(z-la-ld))  *exp(-g3*r):0);
		rB = 0;//B.x*fB;
		iB = 0;//B.y*fB;
		
		sincos(PH, &sinPH, &cosPH);
		sincos(3.*PH, &sinPS, &cosPS);

		Qk1 = dz*(H - k1*EN/sqrt(EN*EN-1.));
		Wk1 = -dz*((rA1*cosPH - iA1*sinPH)+(rB*cosPS-iB*sinPS))*ifnotdestroyed;

//		if(s0 + p0 + q0  == 0 && (i == 1)) printf("%g,%g,%g,%g\n", r, g1, Qk1, Wk1);

/////////////////
		z = ((double)i+0.5)*dz;
		r = (R_center + R_cyclotron*cos(kappa_cyclotron*z + phase_cyclotron));
	//	ifnotdestroyed *= (r > -wall)? 1. : 0.;

		PH = Q + 0.5*Qk1; 
		EN = W + 0.5*Wk1 + en0;

		if(z > la) {rA1 =0; iA1 = 0;}

		fB = ( ((z> la+ld)&&(z < la+ld+lb))?sin(dm_Pi/lb*(z-la-ld))  *exp(-g3*r):0);
		rB = 0;//B.x*fB;
		iB = 0;//B.y*fB;
		
		sincos(PH, &sinPH, &cosPH); 
		sincos(3.*PH, &sinPS, &cosPS);

		Qk2 = dz*(H - k1*EN/sqrt(EN*EN-1.));
		Wk2 = -dz*((rA1*cosPH - iA1*sinPH)+(rB*cosPS-iB*sinPS))*ifnotdestroyed;
/////////////////
		PH = Q + 0.5*Qk2; 
		EN = W + 0.5*Wk2 + en0;

		sincos(PH, &sinPH, &cosPH);	
		sincos(3.*PH, &sinPS, &cosPS);

		Qk3 = dz*(H - k1*EN/sqrt(EN*EN-1.));
		Wk3 = -dz*((rA1*cosPH - iA1*sinPH)+(rB*cosPS-iB*sinPS))*ifnotdestroyed;
/////////////////	
		z = ((double)i+1)*dz;
		r = (R_center + R_cyclotron*cos(kappa_cyclotron*z + phase_cyclotron));
	//	ifnotdestroyed *= (r > -wall)? 1. : 0.;
		PH = Q + Qk3; 
		EN = W + Wk3 + en0;

		if(z > la) {rA1 =0; iA1 = 0;}

		fB = ( ((z> la+ld)&&(z < la+ld+lb))?sin(dm_Pi/lb*(z-la-ld))  *exp(-g3*r):0);
		rB = 0;//B.x*fB;
		iB = 0;//B.y*fB;

		
		sincos(PH, &sinPH, &cosPH);
		sincos(3.*PH, &sinPS, &cosPS);

		Qk4= dz*(H - k1*EN/sqrt(EN*EN-1.));
		Wk4= -dz*((rA1*cosPH - iA1*sinPH)+(rB*cosPS-iB*sinPS))*ifnotdestroyed;
///////////////		
			
		Q +=  1./6.*(Qk1+2.*Qk2+2.*Qk3+Qk4);
		W +=  1./6.*(Wk1+2.*Wk2+2.*Wk3+Wk4);

		__syncthreads();

		sincos(double(Nharm)*Q, &sinPH, &cosPH);

		if(Nharm == 1) 
			fB = exp(-g1*r);			//fB используется как множитель при интегрировании тока вдоль продольной координаты ВНИМАНИЕ!! fB зависит от q0, s0
		else 
			fB = (((z>la+ld)&&(z < la+ld+lb))?sin(dm_Pi/lb*(z-la-ld)):0)*exp(-g3*r);  

		fB *= ifnotdestroyed;

		int xi = p0 + Np*q0 + Np*Nq*s0;
		int X = blockIdx.x + gridDim.x*blockIdx.y + gridDim.y*gridDim.x*blockIdx.z;

		sh_sinQ[xi] = sinPH*fB;
		sh_cosQ[xi] = cosPH*fB;
		

		__syncthreads();
		biReduce(sh_sinQ, sh_cosQ, xi, warpsize, log2warpsize);
	
		if(xi == 0)
		{
			rJ3[X*N+i] =  sh_cosQ[0];	
			iJ3[X*N+i] =  -sh_sinQ[0];


			int_rJ3 += sh_cosQ[0];
			int_iJ3 += -sh_sinQ[0];
		}


		/////////////////////// усреднение энергии
		if(i ==  ifinal) 
		{
			sh_sinQ[xi] = W;
			__syncthreads();
			biReduce(sh_sinQ, sh_cosQ, xi, warpsize, log2warpsize);
			if(xi == 0)
			{	
				avEN = sh_sinQ[0];
			}
			__syncthreads();	
		} 
		///////////////// конец усреднения энергии
		__threadfence();
		__syncthreads();

		if(i > ifinal) break;
	}

	__syncthreads();
	if(p0+q0+s0 == 0)
	{
		*par.avEN = avEN;
		par.int_rJ3[gridDim.x*gridDim.y*blockIdx.z + gridDim.x*blockIdx.y+blockIdx.x] = int_rJ3;
		par.int_iJ3[gridDim.x*gridDim.y*blockIdx.z + gridDim.x*blockIdx.y+blockIdx.x] = int_iJ3;
		par.int_rJ3_1[gridDim.x*gridDim.y*blockIdx.z + gridDim.x*blockIdx.y+blockIdx.x] = int_rJ3_1;
		par.int_iJ3_1[gridDim.x*gridDim.y*blockIdx.z + gridDim.x*blockIdx.y+blockIdx.x] = int_iJ3_1;

	}

}

std::complex<double> Multiplier::retriveBCurr()
{
	int GQ = Nq/NQ; int GS = Ns/NS; int GV = Nv;

	double t_deltaEn[512]; double t_deltaEn2[512];
	double reJ = 0, imJ = 0;

//	printf("memcpy: %i\t", cudaMemcpy((void *)&t_deltaEn,  d_int_rJ3, sizeof(double), cudaMemcpyDeviceToHost));
//	printf("memcpy: %i\n", cudaMemcpy((void *)&t_deltaEn2, d_int_iJ3, sizeof(double), cudaMemcpyDeviceToHost));

	cudaMemcpy((void *) t_deltaEn,  d_int_rJ3, sizeof(double)*GQ*GS*GV, cudaMemcpyDeviceToHost);
	cudaMemcpy((void *) t_deltaEn2, d_int_iJ3, sizeof(double)*GQ*GS*GV, cudaMemcpyDeviceToHost);

	for(int i = 0; i < GQ*GS*GV; i++){
		reJ += t_deltaEn[i]; imJ += t_deltaEn2[i];
	}


	double coeff = Lsolver/double(Nz*Np*Nq*Ns*Nv);

//	printf("re = %g,  im = %g\n", reJ*coeff, imJ*coeff);

	std::complex<double>	 res = std::complex<double>	 (reJ*coeff, imJ*coeff);
	return res;

}
void Multiplier::retriveBCurr(std::complex<double> *J1, std::complex<double> *J2)
{
	int GQ = Nq/NQ; int GS = Ns/NS; int GV = Nv;

	double t_Jre[512]; double t_Jim[512];
	double t_J2re[512]; double t_J2im[512];

	double reJ = 0, imJ = 0;
	double re2J = 0, im2J = 0;

//	printf("memcpy: %i\t", cudaMemcpy((void *)&t_deltaEn,  d_int_rJ3, sizeof(double), cudaMemcpyDeviceToHost));
//	printf("memcpy: %i\n", cudaMemcpy((void *)&t_deltaEn2, d_int_iJ3, sizeof(double), cudaMemcpyDeviceToHost));

	cudaMemcpy((void *) t_Jre,  d_int_rJ3,   sizeof(double)*GQ*GS*GV, cudaMemcpyDeviceToHost);
	cudaMemcpy((void *) t_Jim,  d_int_iJ3,   sizeof(double)*GQ*GS*GV, cudaMemcpyDeviceToHost);
	cudaMemcpy((void *) t_J2re, d_int_rJ3_1, sizeof(double)*GQ*GS*GV, cudaMemcpyDeviceToHost);
	cudaMemcpy((void *) t_J2im, d_int_iJ3_1, sizeof(double)*GQ*GS*GV, cudaMemcpyDeviceToHost);

	for(int i = 0; i < GQ*GS*GV; i++){
		reJ  += t_Jre[i];   imJ  += t_Jim[i];
		re2J += t_J2re[i];  im2J += t_J2im[i];
	
	}



	double coeff = Lsolver/double(Nz*Np*Nq*Ns*Nv);

//	printf("re = %g,  im = %g\n", reJ*coeff, imJ*coeff);

	std::complex<double>	 res1 = std::complex<double>	 (reJ*coeff, imJ*coeff);
	std::complex<double>	 res2 = std::complex<double>	 (re2J*coeff, im2J*coeff);
	*J1 = res1; *J2 = res2;
//	printf("J1 = %g, %g\tJ2 = %g, %g\n", *J1, *J2);


}
double Multiplier::retriveDeltaEnergy()
{
	int GQ = Nq/NQ; int GS = Ns/NS; int GV = Nv;

	double t_deltaEn[512];
	double t_wmax[512];
	double t_wmin[512];
	double averagedEn = 0, wmax = -99999, wmin = 99999;

	cudaMemcpy( t_deltaEn,  d_avEN, sizeof(double)*GQ*GS*GV, cudaMemcpyDeviceToHost);
	cudaMemcpy( t_wmax, d_Wmax, sizeof(double)*GQ*GS*GV, cudaMemcpyDeviceToHost);
	cudaMemcpy( t_wmin, d_Wmin, sizeof(double)*GQ*GS*GV, cudaMemcpyDeviceToHost);

	for(int i = 0; i < GQ*GS*GV; i++)	
	{
		wmax =(wmax > t_wmax[i]) ? wmax : t_wmax[i];
		wmin =(wmin < t_wmin[i]) ? wmin : t_wmin[i];
		averagedEn += t_deltaEn[i];
//		printf("%g\n", t_deltaEn[i]/double(NP*NQ*NS));
		
	}
	double coeff = 1./double(Np*Nq*Ns*Nv);

//	printf("deltaW + = %g \t deltaW - = %g\n", wmax*511000., wmin*511000.);

	

	return averagedEn*coeff;

}

bool Device::initSolver(int nz, double lsolver, double groupSpeedCoeff, char *_solverName)
{
	Nz = nz;
	Lsolver = lsolver;
	Lmax = lsolver;
	solverName = _solverName;
	Nmax = nz;
	MultiplierGroupSpeedCoefficient = groupSpeedCoeff;

	printf("The %s solver is intialized\n", solverName);

	int GQ = Nq/NQ; int GS = Ns/NS; int GV = Nv;

//	printf(" Nq %i, Ns %i, Nv %i \t GQ %i, GS %i, GV %i \n",Nq, Ns, Nv, GQ, GS, GV);

	printf("Nz, Lsolver, grSpeed, %i, %g, %g\n", Nz, Lsolver,MultiplierGroupSpeedCoefficient);

	gpuErrChk(cudaMalloc((void**)&d_rJ3, Nz*GQ*GS*GV*sizeof(double)));
	gpuErrChk(cudaMalloc((void**)&d_iJ3, Nz*GQ*GS*GV*sizeof(double)));

	gpuErrChk(cudaMalloc((void**)&d_Nz, sizeof(int)));
	gpuErrChk(cudaMalloc((void**)&d_Lsolver, sizeof(double)));
	gpuErrChk(cudaMalloc((void**)&d_avEN, sizeof(double)*GQ*GS*GV));
		
	gpuErrChk(cudaMalloc((void**)&d_int_rJ3_1, sizeof(double)*GQ*GS*GV));
	gpuErrChk(cudaMalloc((void**)&d_int_iJ3_1, sizeof(double)*GQ*GS*GV));

	gpuErrChk(cudaMalloc((void**)&d_Amps, sizeof(cplx) * 30));
	
	if(strcmp(solverName,"multiplier_spcharge_2d") != 0)
	{		
		gpuErrChk(cudaMalloc((void**)&d_int_rJ3, sizeof(double)*GQ*GS*GV));
		gpuErrChk(cudaMalloc((void**)&d_int_iJ3, sizeof(double)*GQ*GS*GV));
	}

	gpuErrChk(cudaMalloc((void**)&d_Wmax, sizeof(double)*GQ*GS*GV));
	gpuErrChk(cudaMalloc((void**)&d_Wmin, sizeof(double)*GQ*GS*GV));


	gpuErrChk(cudaMalloc((void**)&d_par, sizeof(PAR)));
	gpuErrChk(cudaMalloc((void**)&grSpeedCoeff, sizeof(double)));
	gpuErrChk(cudaMemcpy((void*)d_Nz, &Nz, sizeof(int), cudaMemcpyHostToDevice));
	gpuErrChk(cudaMemcpy((void*)&grSpeedCoeff, &MultiplierGroupSpeedCoefficient, sizeof(double), cudaMemcpyHostToDevice)); // TODO Here is a bug
	gpuErrChk(cudaMemcpy((void*)d_Lsolver, (void*)&Lsolver, sizeof(double), cudaMemcpyHostToDevice));

	return 1;
}

void Device::releaseDeviceMemory()
{

	cudaFree((void*)d_Nz);
	cudaFree((void*)d_Lsolver);
	cudaFree((void*)d_avEN);
	cudaFree((void*)d_int_rJ3);
	cudaFree((void*)d_int_iJ3);

	if(fieldLoaded)
	{
		cudaFree((void*) d_tAr);
		cudaFree((void*) d_tAi);
	}

}


double Multiplier::DeltaEnergy(double A)
{
	PAR par;
	double d = period;
	double h = 2.*Pi/d;
	double La = period*double(Nperiods);

	par.la = La; par.lb = Lb; par.ld = Ld; par.k1 = k1; par.h = h; par.voltage = voltage;
	par.Nz = Nz; par.L = Lsolver; par.wall = wall;
	par.g1 = g1; par.g3 = g3;
	par.Wmax = d_Wmax; par.Wmin = d_Wmin;

	par.rJ3 = d_rJ3;  par.iJ3 = d_iJ3; 

	par.avEN = d_avEN;
	par.int_rJ3 = d_int_rJ3;
	par.int_iJ3 = d_int_iJ3;




	double2 zero = {0,0};

//	cudaMemcpy(	   d_rJ3, &dzero, sizeof(double), cudaMemcpyHostToDevice);
//	cudaMemcpy(	   d_iJ3, &dzero, sizeof(double), cudaMemcpyHostToDevice);
//	cudaMemcpy(d_int_rJ3, &dzero, sizeof(double), cudaMemcpyHostToDevice);
//	cudaMemcpy(d_int_iJ3, &dzero, sizeof(double), cudaMemcpyHostToDevice);
//	cudaMemcpy(   d_avEN, &dzero, sizeof(double), cudaMemcpyHostToDevice);


	

	dim3 threadsPerBlock(NP, NQ, NS); 

	cudaMemcpy(d_par, &par, sizeof(PAR), cudaMemcpyHostToDevice);
	
	MotionEquationMultiplier<<<dim3((size_t) Nq/NQ,(size_t) Ns/NS,(size_t) Nv), threadsPerBlock>>>(d_par, La, 1, A, zero);

/*	double *debRe = new double [Nz];
	double *debIm = new double [Nz];

	cudaError copy1 = cudaMemcpy((void*) debRe, (void *)dm_rJq, sizeof(double)*Nz, cudaMemcpyDeviceToHost);
	printf("copy1 = %i \n", copy1);
	cudaError_t copy2 = cudaMemcpy((void*) debIm, (void *)dm_iJq, sizeof(double)*Nz, cudaMemcpyDeviceToHost);
	printf("copy2 = %i \n", copy2);

*/


	//printf("memcpy: %i \n", cudaMemcpy((void*) &t_deltaEn, d_avEN, sizeof(double), cudaMemcpyDeviceToHost));

	//printf("Energy delta = %g \n", t_deltaEn/double(NP*NQ*NS));

	double res = retriveDeltaEnergy();

//	printf("Retrieve returned: %g \n", res);

	return res;

//	delete[] debRe; delete[] debIm;


}
std::complex<double> Multiplier::CurrentB(double  reB, double imB, double A)
{
	PAR par;
	double d = period;
	double h = 2.*Pi/d; 
	double La = period*double(Nperiods);


	par.la = La; par.lb = Lb; par.ld = Ld; par.k1 = k1; par.h = h; par.voltage = voltage;
	par.Nz = Nz; par.L = Lsolver; par.wall = wall;
	par.g1 = g1; par.g3 = g3;
	par.Wmax = d_Wmax; par.Wmin = d_Wmin;

//	printf("CurrentB: %g, %g, %g \n", La, Ld, Lb);

	cudaMemset(d_rJ3, 0, sizeof(double)*Nz);
	cudaMemset(d_iJ3, 0, sizeof(double)*Nz);

	par.rJ3 = d_rJ3;  par.iJ3 = d_iJ3; 

	par.avEN = d_avEN;
	par.int_rJ3 = d_int_rJ3;
	par.int_iJ3 = d_int_iJ3;

	double2 B; B.x = reB; B.y = imB;

//	printf("\n B loop: %g\n", La+Ld+Lb );
//	printf("\n Threads: %i, %i, %i\n", threadsPerBlock.x, threadsPerBlock.y, threadsPerBlock.z );

	dim3 numblocks(Nq/NQ, Ns/NS, Nv);
	dim3 threadsPerBlock(NP, NQ, NS);

//	cudaMemcpy(	   d_rJ3, &dzero, sizeof(double), cudaMemcpyHostToDevice);
//	cudaMemcpy(	   d_iJ3, &dzero, sizeof(double), cudaMemcpyHostToDevice);
//	cudaMemcpy(d_int_rJ3, &dzero, sizeof(double), cudaMemcpyHostToDevice);
//	cudaMemcpy(d_int_iJ3, &dzero, sizeof(double), cudaMemcpyHostToDevice);
//	cudaMemcpy(   d_avEN, &dzero, sizeof(double), cudaMemcpyHostToDevice);


	cudaMemcpy(d_par, &par, sizeof(PAR), cudaMemcpyHostToDevice);

	MotionEquationMultiplier << <numblocks, threadsPerBlock >> >(d_par, La + Ld + Lb, 3, A, B);

	double *jr = new double [Nz];
	double *ji = new double [Nz];
	cudaMemcpy(jr, d_rJ3, sizeof(double)*Nz, cudaMemcpyDeviceToHost);
	cudaMemcpy(ji, d_iJ3, sizeof(double)*Nz, cudaMemcpyDeviceToHost);

	FILE *resamp_ar = fopen("F:\\Piotr\\bwo_Data\\mdebug_jr.csv", "w");
	FILE *resamp_ai = fopen("F:\\Piotr\\bwo_Data\\mdebug_ji.csv", "w");
	for(int j = 0; j < Nz; j++)
	{
		fprintf(resamp_ar, "%i,%g\n", j, jr[j]);
		fprintf(resamp_ai, "%i,%g\n", j, ji[j]);
	
	}
	fclose(resamp_ar);
	fclose(resamp_ai);

	delete []jr;
	delete []ji;


	return retriveBCurr();



}
std::complex<double> MultiplierThreeCavity::CurrentB2(double  reB, double imB, double A, cplx A2)
{
	PAR par;
	double d = period;
	double h = 2.*Pi/d;

	double La1 = period*(double)Nperiods;

//	printf("CurrentB2: %g, %g, %g, %g, %g \n", La1, Ld1, La2, Ld2, Lb);

	par.la1 = La1; par.lb = Lb; par.ld1 = Ld1; par.k1 = k1; par.h = h; par.voltage = voltage;
	par.la2 = La2; par.ld2 = Ld2;
	par.Nz = Nz; par.L = Lsolver;
	par.wall = wall; par.g1 = g1; par.g3 = g3;
	par.Wmax = d_Wmax; par.Wmin = d_Wmin;


	par.rJ3 = d_rJ3;  par.iJ3 = d_iJ3; 

	par.avEN = d_avEN;
	par.int_rJ3 = d_int_rJ3;
	par.int_iJ3 = d_int_iJ3;

	double2 B; B.x = reB; B.y = imB;
	double2 Astat2 ={A2.real(), A2.imag()};

	dim3 numblocks(Nq/NQ, Ns/NS, Nv);
	dim3 threadsPerBlock(NP, NQ, NS);

	MotionEquationMultiplierDoubleScheme << <numblocks, threadsPerBlock >> >(par, La1 + Ld1 + La2 + Ld2 + Lb, 3, A, Astat2, B);

	return retriveBCurr();
}

std::complex<double> Multiplier::CurrentA(double  reA, double imA)
{
	PAR par;

	double d = period;
	double h = 2.*Pi/d;
	double La = period*double(Nperiods);

	par.la = La; par.lb = 1.; par.ld = 1.; par.k1 = k1; par.h = h; par.voltage = voltage;
	par.Nz = Nz; par.L = Lsolver; par.wall = wall;
	par.g1 = g1; par.g3 = g3;
	par.Wmax = d_Wmax; par.Wmin = d_Wmin;


	par.rJ3 = d_rJ3;  par.iJ3 = d_iJ3; 

	par.avEN = d_avEN;
	par.int_rJ3 = d_int_rJ3;
	par.int_iJ3 = d_int_iJ3;

	double2 zero = {0,0};

	dim3 threadsPerBlock(NP, NQ, NS); 
	dim3 numblocks(Nq / NQ, Ns / NS, Nv);

	double A; A = sqrt(reA*reA + imA*imA);	

//	printf("\n B loop: %g\n", La+Ld+Lb );
//	printf("\n Threads: %i, %i, %i\n", threadsPerBlock.x, threadsPerBlock.y, threadsPerBlock.z );	

//	cudaMemcpy(	   d_rJ3, &dzero, sizeof(double), cudaMemcpyHostToDevice);
//	cudaMemcpy(	   d_iJ3, &dzero, sizeof(double), cudaMemcpyHostToDevice);
//	cudaMemcpy(d_int_rJ3, &dzero, sizeof(double), cudaMemcpyHostToDevice);
//	cudaMemcpy(d_int_iJ3, &dzero, sizeof(double), cudaMemcpyHostToDevice);
//	cudaMemcpy(   d_avEN, &dzero, sizeof(double), cudaMemcpyHostToDevice);

	gpuErrChk(cudaMemcpy(d_par, &par, sizeof(PAR), cudaMemcpyHostToDevice));

	MotionEquationMultiplier << <numblocks, threadsPerBlock >> >(d_par, La, 1, A, zero);



	return retriveBCurr()*exp(I*arg(reA + I*imA));



}
std::complex<double> MultiplierThreeCavity::CurrentA2(double A1, double  reA, double imA)
{
	PAR par;

	double d = period;
	double h = 2.*Pi/d;
	double La1 = period*double(Nperiods);

	par.la1 = La1;  par.la2 = La2; par.ld1 = Ld1;

	par.lb = 1.; par.ld = 1.; par.k1 = k1; par.h = h; par.voltage = voltage;
	
	par.Nz = Nz; par.L = Lsolver; par.wall = wall;
	
	par.g1 = g1; par.g3 = g3;

	par.Wmax = d_Wmax; par.Wmin = d_Wmin;


	par.rJ3 = d_rJ3;  par.iJ3 = d_iJ3; 

	par.avEN = d_avEN;
	par.int_rJ3 = d_int_rJ3;
	par.int_iJ3 = d_int_iJ3;

	double2 zero = {0,0};
	double2 A = {reA, imA};

	dim3 threadsPerBlock(NP, NQ, NS); 
	dim3 numblocks(Nq/NQ, Ns/NS, Nv);


	MotionEquationMultiplierDoubleScheme <<< numblocks, threadsPerBlock >>>(par, La1 + La2 + Ld1, 1, A1, A, zero);


	return retriveBCurr();


}

void MultiplierMultiModes::CurrentAMultiModes(std::complex<double>  *Amps, std::complex<double> * currs, double *buffRe, double *buffIm, int Na, cplx *J1, cplx *J2)
{
	PAR par;

	double d = period;
	double h = 2.*Pi/d;
	double La = period*double(Nperiods);
	int Nstop = La/dz;

	par.la1 = La;  par.ld = Ld;
	par.lb = 1.; par.k1 = k1; par.h = h; par.voltage = voltage;
	par.Nz = Nz; par.L = Lsolver; par.wall = wall;
	par.g1 = g1; par.g3 = g3;

	par.rJ3 = d_rJ3;  par.iJ3 = d_iJ3; 

	par.avEN = d_avEN;
	par.int_rJ3 = d_int_rJ3;
	par.int_iJ3 = d_int_iJ3;
	par.int_rJ3_1 = d_int_rJ3_1;
	par.int_iJ3_1 = d_int_iJ3_1;
	

	double2 zero = {0,0};

	dim3 threadsPerBlock(NP, NQ, NS); 
	dim3 numblocks(Nq/NQ, Ns/NS, Nv);

	par.Amps = (double2*) d_Amps;

	int ierr = cudaMemcpy(d_Amps, (void*) Amps, sizeof(double2)*Na, cudaMemcpyHostToDevice);

	MotionEquationMultiplierMultiModes <<< numblocks, threadsPerBlock >>>(par, La, 1, Na, zero);
	gpuErrChk(cudaPeekAtLastError());

	retriveACurrComplex((std::complex<double>*)Amps, currs, buffRe, buffIm, Namm, Nstop);

}
void  MultiplierMultiModes::retriveACurrComplex(std::complex<double>  *Amps, std::complex<double>  *currs, double *currsBuffRe, double *currsBuffIm, int Na, int Nstop)
{
	int GQ = Nq / NQ; int GS = Ns / NS; int GV = Nv;

	double reJ = 0, imJ = 0;
	double rF, iF, z;
	double La = period*double(Nperiods);
	std::complex<double> J;


	//	printf("memcpy: %i\t", cudaMemcpy((void *)&t_deltaEn,  d_int_rJ3, sizeof(double), cudaMemcpyDeviceToHost));
	//	printf("memcpy: %i\n", cudaMemcpy((void *)&t_deltaEn2, d_int_iJ3, sizeof(double), cudaMemcpyDeviceToHost));

	gpuErrChk(cudaMemcpy((void *)currsBuffRe, d_rJ3, sizeof(double)*GQ*GS*GV*Nmax, cudaMemcpyDeviceToHost))
	gpuErrChk(cudaMemcpy((void *)currsBuffIm, d_iJ3, sizeof(double)*GQ*GS*GV*Nmax, cudaMemcpyDeviceToHost))

	for (int a = 0; a < Na; a++)
	{
		currs[a] = 0;
	}

	//	FILE* debugfile = fopen("F:\\Piotr\\CalcData\\mm_orotron_Data\\debug.txt", "w");
	for (int j = 0; j < Nstop; j++)
	{

		reJ = 0;  imJ = 0;
		for (int i = 0; i < GQ*GS*GV; i++)
		{
			reJ += currsBuffRe[i*Nmax + j]; imJ += currsBuffIm[i*Nmax + j];
		}
		for (int a = 0; a < Na; a++)
		{
			z = (double)j * dz;
			sincos(Pi / La*z*double(a - Na / 2), &iF, &rF);
			J = cplx(reJ, imJ)*cplx(rF, -iF);
			currs[a] += (J);
			//			if(a == 1) fprintf(debugfile, "%g,%g,%g,%g,%g\n",z, real(J)/double(Np*Nq*Ns*Nv), imag(J)/double(Np*Nq*Ns*Nv), abs(J)/double(Np*Nq*Ns*Nv), arg(J) );
		}
	}

	double coeff = Lsolver / double(Nz*Np*Nq*Ns*Nv);
	for (int a = 0; a < Na; a++) currs[a] *= coeff;
	//	fclose(debugfile);

}





//////////////////////////////////

ParamsM Device::setPar()
{
	ParamsM par;
	int GQ = Nq/NQ; int GS = Ns/NS; int GV = Nv;
	int gridsize = GQ*GS*GV;

	

	double La = Nperiods*period; 
	double h = 2.*Pi/period;

	par.la = La;	 par.k1 = k1; par.h = h; par.voltage = voltage;
	par.Nz = Nmax;	 par.L = Lmax; par.wall = wall;
	par.g1 = g1;	 par.Ngrid = gridsize;
	par.ar0 = d_ar0; par.ai0 = d_ai0;
	par.rJ3 = d_rJ3;  par.iJ3 = d_iJ3; 
	par.delta = 0;
	par.Q0 = d_Q0;	par.W0 = d_W0;

	par.rAk = d_rAk; par.iAk = d_iAk;
	par.rAq1k = d_rAq1k; par.iAq1k = d_iAq1k;

	par.Qk = d_Qk;	par.Wk = d_Wk;
	par.ar0_t = d_ar0_t; par.ai0_t = d_ai0_t;
	par.int_rQ1 = d_int_rQ1;
	par.int_iQ1 = d_int_iQ1;
	par.ifnotdestroyed = d_ifnotdestroyed;
	par.g3 = g3;

	par.rAq1 =d_rAq1;
	par.iAq1 =d_iAq1;
	par.radii = d_radii;

	par.int_rJ3 = d_int_rJ3;
	par.int_iJ3 = d_int_iJ3;
	par.int_rJ3_1 = d_int_rJ3_1;
	par.int_iJ3_1 = d_int_iJ3_1;
	par.avEN = d_avEN;

	int *mass = new int [Np*Nq*Ns*Nv];
	for(int a = 0; a < Np*Nq*Ns*Nv; a++) mass[a] = 1;
	gpuErrChk(cudaMemcpy(d_ifnotdestroyed, mass, sizeof(int)*Np*Nq*Ns*Nv, cudaMemcpyHostToDevice));
	delete [] mass;


	return par;
}