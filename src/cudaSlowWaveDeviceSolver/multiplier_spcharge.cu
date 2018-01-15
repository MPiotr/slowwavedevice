#include <stdio.h>
#include <complex>
#include "cu_mult.h"
#include "multiplier_spcharge.h"
//#include "print.h"
  

#define PAR ParamsM
#ifndef dm_Pi
__device__ const double dm_Pi = 3.141592653589793;
#endif

__global__ void
motionstep_spacecharge(PAR *par, int i, int k, double K, double A, double rB, double iB, int Nharm, double coupling = 1)//Fixed Structure
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
	unsigned int log2warpsize = round(log2((float)warpsize));
	int xi = (p0 + Np*q0 + s0*Np*Nq);


	double h, k1,  voltage, g1, g3;
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
	
	R_center =  0.5*wall + wall*((double)q_init-0.5*(double)Nq_max)/(double)(Nq_max);
	initial_angle = 0;//(0.0810194 - 2.05972*R_center +  28.0433*R_center*R_center);
	R_cyclotron =	(0.568*initial_angle + 0*0.035156*((double)v_init)/double(Nv_max));
	kappa_cyclotron = 1.758;
	phase_cyclotron = 2.*dm_Pi*(double)s_init/(double)Ns_max;

	double en0 = 1. + voltage/511.;	
	double angle_spread_factor = 1./sqrt(1. + initial_angle*initial_angle);




	__shared__ double    sh_sinQ[NS*NQ*NP];
	__shared__ double    sh_cosQ[NS*NQ*NP];

	double PH, EN, cosPH, sinPH, sinPS, cosPS;

	double H = h;//+dh(delta);
/*	double NU = 1./(en0*(en0*en0 - 1.));
	double DELTA = 1.-k1*en0/sqrt(en0*en0 - 1.);*/

	double *d_Qk = par->Qk;	double *Q0 = par->Q0;
	double *d_Wk = par->Wk;	double *W0 = par->W0;
	double *rAq3k = par->rAk;	double *iAq3k = par->iAk;
	double *rAq1k = par->rAq1k;	double *iAq1k = par->iAq1k;

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
		Q = 2.*dm_Pi/double(Np)*double(p0);// + 1./(double)Nq*((double)q0 + (double)s0/(double)Ns));
		W = 0; 
		Q0[gstride + xi] = Q;
		W0[gstride + xi] = W;

	}
	if( (i == 0)&&(xi == 0))
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

	PH = Q + K*Qk;
	EN = W + K*Wk + en0;

	double rak = (k == 0)?0:rAq3k[k-1];
	double iak = (k == 0)?0:iAq3k[k-1];

	rQ3 = (*par->ar0 + K*rak);
	iQ3 = (*par->ai0 + K*iak);

	rak = (k == 0)?0:rAq1k[k-1];
	iak = (k == 0)?0:iAq1k[k-1];

	iQ1 = (*(par->iAq1) + K*iak);
	rQ1 = (*(par->rAq1) + K*rak);

//	if((i < 2)&&(xi == 0)) printf("%i,%i\t%g\t%g\n",i, k, rQ1, iQ1);
	


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

	unsigned int stride;

	sh_sinQ[xi] = DQ*sinPS*double(ifnotdestroyed); //Амлитуда пр. заряда (3-я гармоника)
	sh_cosQ[xi] = DQ*cosPS*double(ifnotdestroyed);
	
	stride = warpsize;
	__syncthreads();
	for(int q = 1; q <= log2warpsize; q++)
	{
		stride =  stride >> 1;
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


	if(xi == 0)
	{
		par->int_rJ3[NX*NY*Z + NX*Y+X] =  -3.*sh_sinQ[0]/warpsize;
		par->int_iJ3[NX*NY*Z + NX*Y+X] =  -3.*sh_cosQ[0]/warpsize; 
	
	}
	__syncthreads();

	sh_sinQ[xi] = DQ*sinPH*double(ifnotdestroyed); //Амлитуда пр. заряда (1-я гармоника)
	sh_cosQ[xi] = DQ*cosPH*double(ifnotdestroyed);
	
	stride = warpsize;
	__syncthreads();
	for(int q = 1; q <= log2warpsize; q++) 
	{
		stride =  stride >> 1;
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


	if(xi == 0)
	{
		par->int_rQ1[NX*NY*Z + NX*Y+X] = -sh_sinQ[0]/warpsize;//??
		par->int_iQ1[NX*NY*Z + NX*Y+X] = -sh_cosQ[0]/warpsize; 
	
	}
	__syncthreads();



}
__global__ void
amplitudestep_spacecharge(PAR *par, int i,  int k, double K, double spchQ1, double spchQ3)
{
	unsigned int xi = threadIdx.x;

	unsigned int warpsize = blockDim.x;
	unsigned int log2warpsize = round(log2((float)warpsize));

	double *rAk = par->rAk;
	double *iAk = par->iAk;
	double *rAq1k = par->rAq1k;
	double *iAq1k = par->iAq1k;

	double *rJ = par->int_rJ3;
	double *iJ = par->int_iJ3;
	double *rQ1 = par->int_rQ1;
	double *iQ1 = par->int_iQ1;

	int N = par->Nz;
	double dz = par->L/(double)N;
//	double rak, iak;


	__shared__ double   shJr[256];
	__shared__ double   shJi[256];

	shJr[xi] = rJ[xi]/warpsize; //чтобы не делить после усреднения
	shJi[xi] = iJ[xi]/warpsize;


	unsigned int stride = warpsize;

	__syncthreads();

	for(int q = 1; q <= log2warpsize; q++)
	{
		stride = stride >> 1;
		if(xi < stride)
		{
			shJr[xi] += shJr[xi + stride];
		}
		else
		{
			if(xi < 2*stride)
			{
				shJi[xi - stride] += shJi[xi];
			}
		}
		__syncthreads();
	}
	

	if(xi == 0)
	{
		rAk[k] = dz*spchQ3*shJr[0];
		iAk[k] = dz*spchQ3*shJi[0];	
	}
	
	__syncthreads(); 

	shJr[xi] = rQ1[xi]/warpsize;
	shJi[xi] = iQ1[xi]/warpsize;

	__syncthreads();

	stride = warpsize;

	for(int q = 1; q <= log2warpsize; q++)
	{
		stride = stride >> 1;
		if(xi < stride)
		{
			shJr[xi] += shJr[xi + stride];
		}
		else
		{
			if(xi < 2*stride)
			{
				shJi[xi - stride] += shJi[xi];
			}
		}
		__syncthreads();
	}
	

	if(xi == 0)
	{
		
		rAq1k[k] = dz*spchQ1*shJr[0];
		iAq1k[k] = dz*spchQ1*shJi[0]; 
	
	}



}
__global__ void
spacecharge_endstep(PAR *par, int i, int N_averaging, double fB = 0)
{
	unsigned int p0 = threadIdx.x;	    unsigned int Np = blockDim.x; 	unsigned int X_ = blockIdx.x;		unsigned int NX = gridDim.x;
	unsigned int q0 = threadIdx.y;		unsigned int Nq = blockDim.y;  	unsigned int Y_ = blockIdx.y;		unsigned int NY = gridDim.y;
	unsigned int s0 = threadIdx.z;		unsigned int Ns = blockDim.z;  	unsigned int Z_ = blockIdx.z;	//	unsigned int NZ = gridDim.z;

	unsigned int q_init =		 Nq*X_ + q0;		unsigned int Nq_max =		 Nq*gridDim.x;
	unsigned int s_init =		 Ns*Y_ + s0;		unsigned int Ns_max =		 Ns*gridDim.y;
	unsigned int v_init =			Z_;				unsigned int Nv_max =		    gridDim.z;

	double wall = par->wall; 

	double R_center =  0.5*wall + wall*((double)q_init-0.5*(double)Nq_max)/(double)(Nq_max);
	double initial_angle = 0;//(0.0810194 - 2.05972*R_center +  28.0433*R_center*R_center);
	double R_cyclotron =	(0.568*initial_angle + 0*0.035156*((double)v_init)/double(Nv_max));
	double kappa_cyclotron = 1.758;
	double phase_cyclotron = 2.*dm_Pi*(double)s_init/(double)Ns_max;
	
	double z = par->L/double(par->Nz)*double(i);
	double r = (R_center + R_cyclotron*cos(kappa_cyclotron*z + phase_cyclotron));
	double g3 = par->g3;

	unsigned int warpsize = Ns*Nq*Np;
	unsigned int log2warpsize = round(log2((float)warpsize));
	int xi = (p0 + Np*q0 + s0*Np*Nq);

	
//	unsigned int xi = threadIdx.x;	    unsigned int Np = blockDim.x;
	unsigned int X = X_+NX*Y_+NX*NY*Z_;	

//	unsigned int warpsize = Np;
//	unsigned int gridsize = NX*NY*NZ;
	int Nz = par->Nz;
	

	double *d_Qk = par->Qk; double *Q0 = par->Q0;
	double *d_Wk = par->Wk; double *W0 = par->W0;
	double *rAq3k = par->rAk; double *iAq3k = par->iAk;
	double *rAq1k = par->rAq1k; double *iAq1k = par->iAq1k;
	double *rJ3 = par->rJ3; double *iJ3 = par->iJ3;

	unsigned int stride = 4*X*warpsize;
	unsigned int gstride = X*warpsize;

	if(r < -wall) par->ifnotdestroyed[gstride + xi] *= 0;
	int ifnotdestroyed = par->ifnotdestroyed[gstride + xi]; 
	
	Q0[gstride + xi] += (d_Qk[stride + xi] + 2.*d_Qk[stride+ xi + warpsize] + 2.*d_Qk[stride+ xi + 2*warpsize] + d_Qk[stride+ xi+ 3*warpsize])/6.;
	W0[gstride + xi] += (d_Wk[stride + xi] + 2.*d_Wk[stride +xi + warpsize] + 2.*d_Wk[stride +xi + 2*warpsize] + d_Wk[stride +xi+ 3*warpsize])/6.;


	
	if(X + xi == 0)
	{
		*par->ar0 += (rAq3k[0] + 2.*rAq3k[1] + 2.*rAq3k[2] + rAq3k[3])/6.;
		*par->ai0 += (iAq3k[0] + 2.*iAq3k[1] + 2.*iAq3k[2] + iAq3k[3])/6.;

		*par->rAq1+= (rAq1k[0] + 2.*rAq1k[1] + 2.*rAq1k[2] + rAq1k[3])/6.;
		*par->iAq1+= (iAq1k[0] + 2.*iAq1k[1] + 2.*iAq1k[2] + iAq1k[3])/6.;

	}
	
	__shared__ double sh_sinQ[NS*NQ*NP];
	__shared__ double sh_cosQ[NS*NQ*NP];

	double sinPS, cosPS;

	sincos(3.*Q0[gstride + xi], &sinPS, &cosPS);

	sh_sinQ[xi] = sinPS*exp(-g3*r)*double(ifnotdestroyed);
	sh_cosQ[xi] = cosPS*exp(-g3*r)*double(ifnotdestroyed);

	__syncthreads();	
	stride = warpsize;

	for(int q = 1; q <= log2warpsize; q++) //Это амплитуда ВЧ тока (?)
	{
		stride =  stride >> 1;
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

//	if((i == 1300)) printf("%g\t%g\n", g3, Q0[gstride + xi]);
	

	if(xi == 0)
	{
		rJ3[X*Nz + i] =  sh_cosQ[0]/warpsize;
		iJ3[X*Nz + i] = -sh_sinQ[0]/warpsize;

	//	if((i == 1300)) printf("\n%g\n", sh_cosQ[0]);

		par->int_rJ3_1[X] +=     sh_cosQ[0]/warpsize*fB;
		par->int_iJ3_1[X] +=    -sh_sinQ[0]/warpsize*fB; 

	}
	__syncthreads();

	if(i == N_averaging)
	{
		sh_sinQ[xi] = W0[gstride + xi];

		stride = warpsize;
		for(int q = 1; q <= log2warpsize; q++)
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

Multipler_SpCharge::Multipler_SpCharge()
{
	;
}
Multipler_SpCharge::Multipler_SpCharge(QDomDocument *doc) : Multiplier(doc)
{
	;
}
double  Multipler_SpCharge::DeltaEnergySpaceCharge(double A)
{
	int GQ = Nq/NQ; int GS = Ns/NS; int GV = Nv;

//	int warpsize = NQ*NS*NP;
	int gridsize = GQ*GS*GV;

	ParamsM par = setPar();

	double K[4] = {0, 0.5, 0.5, 1.};

	dim3 motionwarp = dim3(NP, NQ, NS);
	dim3 motiongrid = dim3(GQ, GS, GV);

	cudaMemset(d_Qk, 0, 4*sizeof(double)*Np*Nq*Ns*Nv);	cudaMemset(d_W0, 0, sizeof(double)*Np*Nq*Ns*Nv);
	cudaMemset(d_Wk, 0, 4*sizeof(double)*Np*Nq*Ns*Nv);	cudaMemset(d_rAk, 0, sizeof(double)*4);
	cudaMemset(d_Q0, 0, sizeof(double)*Np*Nq*Ns*Nv);	cudaMemset(d_iAk, 0, sizeof(double)*4);
	
	cudaMemset(d_ai0, 0, sizeof(double));
	cudaMemset(d_ar0, 0, sizeof(double));
	cudaMemset(d_rAq1, 0, sizeof(double));
	cudaMemset(d_iAq1, 0, sizeof(double));

	cudaMemset(d_avEN, 0, sizeof(double));


	double la = Nperiods*period;
	int Nz = Nmax;
	double dz = Lmax/double(Nz);

	int Nstop = floor(la/dz);
	double z, rB = 0, iB = 0;
	double res = 0; 
	double avEN[512];	

	cudaMemcpy(d_par, &par, sizeof(PAR), cudaMemcpyHostToDevice);

	double tmpReRho, tmpImRho, tmpReRho3, tmpImRho3;
	FILE *file =  fopen("F:\\Piotr\\Multiplier_260GHz_Data\\rhoDeb.txt", "w");
		
	for(int i = 0; i < Nstop; i++)
	{
		z = double(i)*dz;

		for(int k = 0; k < 4; k++)
		{
			motionstep_spacecharge <<<motiongrid,motionwarp>>>(d_par, i, k, K[k],  A*sin(Pi/la*(z+K[k]*dz)), rB, iB, 1);
			amplitudestep_spacecharge  <<<1, gridsize>>>(d_par, i, k, K[k], spchQ1, spchQ3);
		}

//		spacecharge_endstep<<<gridsize, warpsize>>>(d_par, i);
		spacecharge_endstep<<<motiongrid, motionwarp>>>(d_par, i, Nstop-1);
		cudaMemcpy(&tmpReRho, d_rAq1, sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(&tmpImRho, d_iAq1, sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(&tmpReRho3, d_rJ3 + i, sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(&tmpImRho3, d_iJ3 + i, sizeof(double), cudaMemcpyDeviceToHost);
		fprintf(file, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", z, tmpReRho, tmpImRho, sqrt(tmpReRho*tmpReRho+tmpImRho*tmpImRho), tmpReRho3*dz, tmpImRho3*dz, dz*sqrt(tmpReRho3*tmpReRho3+tmpImRho3*tmpImRho3));
	}
	fclose(file);
	res = 0;
	cudaMemcpy(&avEN, (void*) d_avEN, sizeof(double)*gridsize, cudaMemcpyDeviceToHost);

//	printf("\n\n");
	for(int q = 0; q < gridsize; q++) 
	{
		res += avEN[q];
//		printf("%i\t%g\n",q, avEN[q]);
	}
	


		
	res /= gridsize;
	




	return res;
	


}
std::complex<double>  Multipler_SpCharge::CurrentBSpaceCharge(double rB, double iB, double A, double enPrint)
{
	int GQ = Nq/NQ; int GS = Ns/NS; int GV = Nv;

//	int warpsize = NP*NQ*NS;
	int gridsize = GQ*GS*GV;

	ParamsM par = setPar();

	double K[4] = {0, 0.5, 0.5, 1.};

	dim3 motionwarp = dim3(NP, NQ, NS);
	dim3 motiongrid = dim3(GQ, GS, GV);

	cudaMemset(d_Qk, 0, 4*sizeof(double)*Np*Nq*Ns*Nv);	cudaMemset(d_W0, 0, sizeof(double)*Np*Nq*Ns*Nv);
	cudaMemset(d_Wk, 0, 4*sizeof(double)*Np*Nq*Ns*Nv);	cudaMemset(d_rAk, 0, sizeof(double)*4);
	cudaMemset(d_Q0, 0, sizeof(double)*Np*Nq*Ns*Nv);	cudaMemset(d_iAk, 0, sizeof(double)*4);
	
	cudaMemset(d_ai0, 0, sizeof(double));
	cudaMemset(d_ar0, 0, sizeof(double));
	cudaMemset(d_rAq1, 0, sizeof(double));
	cudaMemset(d_iAq1, 0, sizeof(double));

	cudaMemset(d_avEN, 0, sizeof(double));


	double La = Nperiods*period;
	double dz = Lmax/double(Nmax);

/*	double *energy = new double [Np*Nq*Ns*Nv];
	double *phase = new double [Np*Nq*Ns*Nv];
	FILE *res_enr = fopen("F:\\Piotr\\bwo_Data\\enPhase.csv", "w");*/
	
	
	int Nstop = ceil((La + Lb + Ld)/dz);
	int nd = round(period/dz);

	double z, fA, fB;

	cudaMemcpy(d_par, &par, sizeof(PAR), cudaMemcpyHostToDevice);

	double tmpReRho, tmpImRho, tmpReRho3, tmpImRho3;
	FILE *file =  fopen("F:\\Piotr\\CalcData\\twt_data\\rhoDeb.txt", "w");
		
	for(int i = 0; i < Nstop; i++)
	{

		for(int k = 0; k < 4; k++)
		{
			z = double(i)*dz + K[k]*dz;
			//if(z < La*0.5) fA = 0.6*sin(2.*Pi/La*z) + 0.08; else sin(2.*Pi/La*(z - La*0.5))+0.08;
			fA = sin(Pi/La*z);
			fA *= A;
			if(z > La) fA = 0;
			if((z > La+Ld)&&(z < La+ Ld+ Lb)) fB = sin(Pi*(z - La - Ld)/Lb); else fB = 0;
			motionstep_spacecharge <<<motiongrid,motionwarp>>>(d_par, i, k, K[k], fA, rB*fB, iB*fB, 1);
			amplitudestep_spacecharge<<<1, gridsize>>>(d_par, i, k, K[k], spchQ1, spchQ3);
		}
//		spacecharge_endstep<<<gridsize, warpsize>>>(d_par, i, fB);
		spacecharge_endstep<<<motiongrid, motionwarp>>>(d_par, i, 0, fB);

		cudaMemcpy(&tmpReRho, d_rAq1, sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(&tmpImRho, d_iAq1, sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(&tmpReRho3, d_rJ3 + i, sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(&tmpImRho3, d_iJ3 + i, sizeof(double), cudaMemcpyDeviceToHost);
		fprintf(file, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", z, tmpReRho, tmpImRho, sqrt(tmpReRho*tmpReRho+tmpImRho*tmpImRho), tmpReRho3, tmpImRho3, sqrt(tmpReRho3*tmpReRho3+tmpImRho3*tmpImRho3));
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

	fclose(file);
		
	double rJ = 0, iJ = 0;	
	double rJ_1[512];
	double iJ_1[512];

	cudaMemcpy(rJ_1, (void*) d_int_rJ3_1, gridsize*sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(iJ_1, (void*) d_int_iJ3_1, gridsize*sizeof(double), cudaMemcpyDeviceToHost);

	for(int i = 0; i < gridsize; i++) {rJ += rJ_1[i]; iJ += iJ_1[i];}
//	printf("rJ = %g\tiJ = %g\n", rJ/gridsize*dz, iJ/gridsize*dz);

/*	double *jr = new double [Nz];
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
	delete []ji;*/

/*	delete [] energy;
	delete [] phase;
	fclose(res_enr);*/



	return std::complex<double>(rJ/gridsize*dz, iJ/gridsize*dz);
	


}
double Multipler_SpCharge::ElectronsDeltaEnergy(double A)
{
	return DeltaEnergySpaceCharge(A);
}
cplx Multipler_SpCharge::HfElectronCurrent(double  rB, double iB, double Astat)
{
	return CurrentBSpaceCharge(rB, iB, Astat);
}
cplx Multipler_SpCharge::ElectronCurrentA(double reA, double imA)
{
	return 0;
}
bool Multipler_SpCharge::initMultiplierSolver(int nz, double lsolver, double groupSpeedCoeff, char *_solverName)
{
	Multiplier::initMultiplierSolver(nz, lsolver, groupSpeedCoeff, _solverName);
	int GQ = Nq / NQ; int GS = Ns / NS; int GV = Nv;
	printf("Allocate space charge arrays ...\n<font color = \"red\">");
	gpuErrChk(cudaMalloc((void**)&d_Qk, 4 * sizeof(double)*Np*Nq*Ns*Nv))
	gpuErrChk(cudaMalloc((void**)&d_Wk, 4 * sizeof(double)*Np*Nq*Ns*Nv))
	gpuErrChk(cudaMalloc((void**)&d_Q0, sizeof(double)*Np*Nq*Ns*Nv))
	gpuErrChk(cudaMalloc((void**)&d_W0, sizeof(double)*Np*Nq*Ns*Nv))
	gpuErrChk(cudaMalloc((void**)&d_ifnotdestroyed, sizeof(int)*Np*Nq*Ns*Nv))

	gpuErrChk(cudaMalloc((void**)&d_rAq1, sizeof(double)))
	gpuErrChk(cudaMalloc((void**)&d_iAq1, sizeof(double)))
	gpuErrChk(cudaMalloc((void**)&d_ar0, sizeof(double)))
	gpuErrChk(cudaMalloc((void**)&d_ai0, sizeof(double)))
	gpuErrChk(cudaMalloc((void**)&d_int_rQ1, sizeof(double)*GQ*GS*GV))
	gpuErrChk(cudaMalloc((void**)&d_int_iQ1, sizeof(double)*GQ*GS*GV))
	gpuErrChk( cudaMalloc((void**)&d_ar0_t, Nt*sizeof(double)))
	gpuErrChk( cudaMalloc((void**)&d_ai0_t, Nt*sizeof(double)))
	gpuErrChk( cudaMalloc((void**)&d_rAk, 4 * sizeof(double)))
	gpuErrChk( cudaMalloc((void**)&d_iAk, 4 * sizeof(double)))
	gpuErrChk(cudaMalloc((void**)&d_rAq1k, 4 * sizeof(double)))
	gpuErrChk(cudaMalloc((void**)&d_iAq1k, 4 * sizeof(double)))
	printf("<\font>... end allocation\n");


	return 1;
}