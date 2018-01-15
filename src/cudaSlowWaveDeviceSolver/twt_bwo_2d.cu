#include <stdio.h>
#include <complex>
#include "Multiplier.h"
#include "twt_2d.h"

//#include "cu_mult.h"

#define PAR ParamsM
#ifndef dm_Pi
__device__ const double dm_Pi = 3.141592653589793;
#endif

//extern __device__ void biAverage(double *A, double *B, int p0, int datasize, int logsize);
static __device__ void biAverage(double *A, double *B, int p0, int datasize, int logsize)
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


__global__ void particleStep(PAR *par, double *Ar, double *Ai, int i, int k, double K, double coupling = 1)
{
	unsigned int p0 = threadIdx.x;			unsigned int Np = blockDim.x;    //x and y -- are coordinates of beam x-section
	unsigned int x0 = threadIdx.y;			unsigned int Nx = blockDim.y;    //p -- phase
	unsigned int y0 = threadIdx.z;		    unsigned int Ny = blockDim.z;    
																			 
	unsigned int P = blockIdx.x;		    unsigned int NP_= gridDim.x;
	unsigned int X = blockIdx.y;			unsigned int NX = gridDim.y;
	unsigned int Y = blockIdx.z;			unsigned int NY = gridDim.z;

	unsigned int p =		 Np*P + p0;		unsigned int Np_max = Np*NP_;
	unsigned int x =		 Nx*X + x0;		unsigned int Nx_max = Nx*NX;
	unsigned int y =		 Ny*Y + y0;		unsigned int Ny_max = Ny*NY;

	unsigned int warpsize = Nx*Ny*Np;
	unsigned int log2warpsize = round(log2((double)warpsize));
	int xi = (p0 + Np*x0 + y0*Np*Nx);           //развернутый индекс
	int gstride = (NX*NY*P + NX*Y + X)*warpsize;

//...... load problem parameters ......................... 
	double h, k1, voltage;
	h  = par->h; 
	k1 = par->k1; 
	voltage = par->voltage;
	double en0 = 1. + voltage / 511.;
	int N = par->Nz;
	double dz = par->L / (double)N;
// .... get acess to values arrays..................
	double *d_Qk = par->Qk;   double *d_Wk = par->Wk;   //k - компоненты фазы и энергии
	double *Q0 =   par->Q0;   double *W0 =   par->W0;   //фаза and энергия
	double rA = par->ar0[i];  double iA = par->ai0[i];  //амплитуда
// ..... init shared memory ...........................

	__shared__ double    sh_imJ[NS*NQ*NP];
	__shared__ double    sh_reJ[NS*NQ*NP];
	
	double Q, Qk;
	double W, Wk;
	double PH, EN, z, sinPH, cosPH;
	double frA = Ar[x + y*Nx_max]*coupling;
	double fiA = Ai[x + y*Nx_max]*coupling;
// ...... init phase and energy ............	
	if(i != 0)
	{
		Q = Q0[gstride + xi];
		W = W0[gstride + xi];
	}
	else
	{ 
		Q =  2.*dm_Pi / double(Np)*double(p0);
		W = 0; 
		Q0[gstride + xi] = Q;
		W0[gstride + xi] = W;

	}
// ....... at the left end of the space, set HF current to zero
	if( (i == 0)&&(xi ==  0))
	{
		par->avEN[NX*Y+X] = 0; 
		par->int_rJ3[NX*Y+X] = 0;
		par->int_iJ3[NX*Y+X] = 0;
	}
	__syncthreads();
// .......init k coefficient of phase and energy 

   	Wk = (k > 0) ? d_Wk[4*gstride + warpsize*(k-1) + xi] : 0;
	Qk = (k > 0) ? d_Qk[4*gstride + warpsize*(k-1) + xi] : 0;

// --------- end of initialization phase ----------

// --------- make a step -------------------------

	z = ((double)i+ K)*dz; 

	PH = Q + K*Qk;
	EN = W + K*Wk + en0;

//	if((i < 2)&&(p0 == 0)&&(q0 < 2)) printf("%i,%i,%i\t%g\t%g\t%g\t%g\n",i, k, q0, rQ1, iQ1,  iq0, iak);

	sincos(PH, &sinPH, &cosPH);


//...........Уравнения движения определены здесь...............................
	double DQ = (1. - (k1/h)*EN/sqrt(EN*EN-1.));	
	double reF = frA*cosPH - fiA*sinPH;
	double imF = frA*sinPH + fiA*cosPH;
	Qk = dz*DQ;
	Wk = -dz*(rA*reF - iA*imF);
//....................................................................

	d_Qk[4*gstride + warpsize*k + xi] = Qk;
	d_Wk[4*gstride + warpsize*k + xi] = Wk;

// ===========Усреднение=======================	
	
	sh_imJ[xi] =  -imF / double(warpsize);
	sh_reJ[xi] =  reF / double(warpsize);

	__syncthreads();
	

	biAverage(sh_imJ, sh_reJ, xi, warpsize, log2warpsize); // averaged values are stored in zero indexes
	
//========== Запись=====================

	if(xi == 0)
	{
		par->int_iJ3[NX*NY*P + NX*Y+X] =  sh_imJ[0];
		par->int_rJ3[NX*NY*P + NX*Y+X] =  sh_reJ[0]; 	
	}
	__syncthreads();

}

__global__ void amplitudeStep(PAR *par, int i, int k, int prevK, double kappa)
{
	unsigned int xi = threadIdx.x;

	unsigned int warpsize = blockDim.x;
	unsigned int log2warpsize = round(log2((double)warpsize));

	double *rAk = par->rAk;
	double *iAk = par->iAk;

	double *rJ = par->int_rJ3;
	double *iJ = par->int_iJ3;
	double G = par->G;
//	double kappa = par->lossKappa;
	double delta = par->delta;

	int N = par->Nz;
	double dz = par->L/(double)N;
//	double rak, iak;


	__shared__ double   shJr[256];
	__shared__ double   shJi[256];

	shJr[xi] = rJ[xi] / warpsize; //чтобы не делить после усреднения
	shJi[xi] = iJ[xi] / warpsize;


//	unsigned int stride = warpsize;

	__syncthreads();
//....... усреднение по блоку (по варпу усреднено на пред. шаге)................
	biAverage(shJr, shJi, xi, warpsize, log2warpsize);
//..............................................................................
	
//.............. уравнения на амплитуду определены здесь .........................
	if(xi == 0)
	{
		par->rJ3[0] = shJr[0];
		par->iJ3[0] = shJi[0];
		double arr = par->ar0[i];
		double aii = par->ai0[i];
		if (k > 0)
		{
			rAk[k] = dz*(G*shJr[0] - kappa*(arr+ prevK*rAk[k - 1]) + delta*(aii + prevK*iAk[k-1]) );
			iAk[k] = dz*(G*shJi[0] - kappa*(aii+ prevK*iAk[k - 1]) - delta*(arr + prevK*rAk[k-1]) );
		}
		else
		{
			rAk[k] = dz*(G*shJr[0] - kappa*arr + delta*aii);
			iAk[k] = dz*(G*shJi[0] - kappa*aii - delta*arr);
		}
	}
//.....................................................
	
	__syncthreads(); 

}

__global__ void endstep(PAR *par, int i)
{
	// make an addition to the phase, energy, and amlitude in the end of one Runge-Kutta (4,4) step.
	unsigned int p0 = threadIdx.x;			unsigned int Np = blockDim.x;
	unsigned int x0 = threadIdx.y;			unsigned int Nx = blockDim.y;
	unsigned int y0 = threadIdx.z;		    unsigned int Ny = blockDim.z;

	unsigned int P = blockIdx.x;		unsigned int NP_ = gridDim.x;
	unsigned int X = blockIdx.y;		unsigned int NX = gridDim.y;
	unsigned int Y = blockIdx.z;		unsigned int NY = gridDim.z;

/*	unsigned int p =		 Np*P + p0;		    unsigned int Np_max =		 Np*NP_;
	unsigned int x =		 Nx*X + x0;		    unsigned int Nx_max =		 Nx*NX;
	unsigned int y =		 Ny*Y + y0;			unsigned int Ny_max =		 Ny*NY;*/

	unsigned int warpsize = Nx*Ny*Np;
	unsigned int log2warpsize = round(log2((double)warpsize));
	int xi = (p0 + Np*x0 + y0*Np*Nx);
	unsigned int gridsize = NX*NY*NP_;
	unsigned int XI = (NX*NY*P + NX*Y + X);
// ............ init problem parameters ...................................

//	unsigned int xi = threadIdx.x;	    unsigned int Np = blockDim.x;
//	unsigned int warpsize = Np;

// .................... get access to phase, energy and current  k-coefficients
	double *d_Qk = par->Qk;		double *Q0 = par->Q0;
	double *d_Wk = par->Wk;		double *W0 = par->W0;
	double *rAk = par->rAk;    double *iAk = par->iAk;   //k - компоненты амплитуды

	unsigned int stride = 4*XI*warpsize; 
	unsigned int gstride = XI*warpsize;

// ............... calculate next phase and energy .....................
	
	Q0[gstride + xi] += (d_Qk[stride + xi] + 2.*d_Qk[stride+ xi + warpsize] + 2.*d_Qk[stride+ xi + 2*warpsize] + d_Qk[stride+ xi+ 3*warpsize])/6.;
	W0[gstride + xi] += (d_Wk[stride + xi] + 2.*d_Wk[stride +xi + warpsize] + 2.*d_Wk[stride +xi + 2*warpsize] + d_Wk[stride +xi+ 3*warpsize])/6.;

// ............... calculate next amplitude ............................	
	if(X + xi == 0)
	{
		par->ar0[i + 1] = par->ar0[i] +(rAk[0] + 2.*rAk[1] + 2.*rAk[2] + rAk[3]) / 6.;
		par->ai0[i + 1] = par->ai0[i] +(iAk[0] + 2.*iAk[1] + 2.*iAk[2] + iAk[3]) / 6.;
	}
	
}

std::complex<double>  TWT_2D::solveTWT_2d(std::complex<double>  *A, double *Ar, double *Ai, double inputAmp, double lossKappa, double delta,
											double *fieldStructureRe, double *fieldStructureIm, double G, double enPrint, bool printField, double *lStrRe, double *lStrIm, double *qStr)
{
	int GX = Nq / NQ; int GY = Ns / NS; int GP = Nv;

	int warpsize = NP*NQ*NS;
	int gridsize = GP*GX*GY;

	double K[4] = { 0, 0.5, 0.5, 1. };

	dim3 motionwarp = dim3(NP, NQ, NS);
	dim3 motiongrid = dim3(GP, GX, GY);

	gpuErrChk(cudaMemset(d_Qk, 0, 4 * sizeof(double)*Np*Nq*Ns*Nv));	gpuErrChk(cudaMemset(d_W0, 0, sizeof(double)*Np*Nq*Ns*Nv));
	gpuErrChk(cudaMemset(d_Wk, 0, 4 * sizeof(double)*Np*Nq*Ns*Nv));	gpuErrChk(cudaMemset(d_rAk, 0, sizeof(double) * 4 * gridsize));
	gpuErrChk(cudaMemset(d_Q0, 0, sizeof(double)*Np*Nq*Ns*Nv));	    gpuErrChk(cudaMemset(d_iAk, 0, sizeof(double) * 4 * gridsize));

	gpuErrChk(cudaMemset(d_ai0, 0, sizeof(double)*Nmax));
	gpuErrChk(cudaMemset(d_ar0, 0, sizeof(double)*Nmax));
	gpuErrChk(cudaMemset(d_avEN, 0, sizeof(double)));

	gpuErrChk(cudaMemcpy(d_fAr, fieldStructureRe, Nq*Ns*sizeof(double), cudaMemcpyHostToDevice));
	gpuErrChk(cudaMemcpy(d_fAi, fieldStructureIm, Nq*Ns*sizeof(double), cudaMemcpyHostToDevice));

	double h = 2 * Pi / period*(synch_angle / 360.);
	double La = Nperiods*period*h;
	double dz = Lmax / double(Nmax);

	ParamsM par = setPar();
	par.G = G;
	par.h = h;
	par.delta = delta;
	par.lossKappa = lossKappa;

	double *energy = new double [Np*Nq*Ns*Nv];
	double *phase = new double [Np*Nq*Ns*Nv];
	
	int Nstop = ceil(La / dz);
	if (Nstop >= Nmax) Nstop = Nmax - 1;
	int nd = round(period / dz);

	double z, coupl, loss;
	
	Ar[0] = inputAmp; Ai[0] = 0;

	cudaMemcpy(d_ai0, Ai, sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_ar0, Ar, sizeof(double), cudaMemcpyHostToDevice);

	cudaMemcpy(d_par, &par, sizeof(PAR), cudaMemcpyHostToDevice);


	FILE *file;
	FILE *res_enr;

 	if (printField) {
		file = fopen("F:\\Piotr\\CalcData\\twt_data\\twt_debug.txt", "w");
		res_enr = fopen("F:\\Piotr\\CalcData\\twt_data\\enPhase.csv", "w");
	}

	for (int i = 0; i < Nstop; i++)
	{

		for (int k = 0; k < 4; k++)
		{
			z = double(i)*dz + K[k] * dz;
			if (lStrRe) coupl = (1. - K[k])*lStrRe[i] + K[k] * lStrRe[i + 1]; else coupl = 1;
			if (qStr) loss = lossKappa/((1. - K[k])*qStr[i] + K[k] * qStr[i + 1]); else loss = lossKappa;
			if (!isfinite(loss)) loss = lossKappa;
			particleStep  <<<motiongrid, motionwarp >>>(d_par, d_fAr, d_fAi, i, k, K[k], coupl);
			amplitudeStep <<<1, gridsize >>>(d_par, i, k, K[k], loss);

			//cudaMemcpy(&rJ, d_rJ3, sizeof(double), cudaMemcpyDeviceToHost);
			//cudaMemcpy(&iJ, d_iJ3, sizeof(double), cudaMemcpyDeviceToHost);
			//fprintf(res_enr, "%g,%g,%g\n", z, rJ, iJ);
		}
		endstep <<<motiongrid, motionwarp >>>(d_par, i);
		
	}
	gpuErrChk(cudaMemcpy(Ar, d_ar0, sizeof(double)*Nmax, cudaMemcpyDeviceToHost));
	gpuErrChk(cudaMemcpy(Ai, d_ai0, sizeof(double)*Nmax, cudaMemcpyDeviceToHost));

	for (int i = 0; i < Nstop; i++)
	{
		z = double(i)*dz;
		A[i] = complex<double>(Ar[i], Ai[i]);
		if (printField)
			fprintf(file, "%g,%g,%g,%g,%g\n", z, Ar[i], Ai[i], abs(A[i]));
	}


	if (printField){
		fclose(file);
		fclose(res_enr);
	}


	return std::complex<double>(Ar[Nstop-1], Ai[Nstop-1]);



}

bool TWT::initSolver(int nz, double lsolver)
{
	int Nz = nz;
	double Lsolver = lsolver;
	
	printf("\nInitializing TWT solver....\n<font color = \"red\">");

	int GX = Nq / NQ; int GY = Ns / NS; int GP = Nv;
	int gridSize = GX*GY*GP;

	printf("Nz = %i, Lsolver = %g\n", Nz, Lsolver);

	gpuErrChk(cudaMalloc((void**)&d_avEN, sizeof(double)*gridSize))
	gpuErrChk(cudaMalloc((void**)&d_rJ3, sizeof(double)*gridSize))
	gpuErrChk(cudaMalloc((void**)&d_iJ3, sizeof(double)*gridSize))
	gpuErrChk(cudaMalloc((void**)&d_par, sizeof(PAR)))
			
	gpuErrChk(cudaMalloc((void**)&d_Qk, 4 * sizeof(double)*Np*Nq*Ns*Nv))
	gpuErrChk(cudaMalloc((void**)&d_Wk, 4 * sizeof(double)*Np*Nq*Ns*Nv))
	gpuErrChk(cudaMalloc((void**)&d_Q0, sizeof(double)*Np*Nq*Ns*Nv))
	gpuErrChk(cudaMalloc((void**)&d_W0, sizeof(double)*Np*Nq*Ns*Nv))
	gpuErrChk(cudaMalloc((void**)&d_ifnotdestroyed, sizeof(int)*Np*Nq*Ns*Nv))

	gpuErrChk(cudaMalloc((void**)&d_ar0, (Nmax + 1)*sizeof(double)))
	gpuErrChk(cudaMalloc((void**)&d_ai0, (Nmax + 1)*sizeof(double)))
	gpuErrChk(cudaMalloc((void**)&d_int_rJ3, sizeof(double)*Nq*Ns*Nv))
	gpuErrChk(cudaMalloc((void**)&d_int_iJ3, sizeof(double)*Nq*Ns*Nv))

	gpuErrChk(cudaMalloc((void**)&d_rAk, 4 * gridSize*sizeof(double)))
	gpuErrChk(cudaMalloc((void**)&d_iAk, 4 * gridSize*sizeof(double)))
	
	printf("</font>....... End initialization\n\n");
		
	return 1;
}
bool TWT_2D::initSolver(int nz, double lsolver)
{
//	TWT::initSolver(nz, lsolver);

	gpuErrChk(cudaMalloc((void**)&d_fAr, Nq*Ns*sizeof(double)))
	gpuErrChk(cudaMalloc((void**)&d_fAi, Nq*Ns*sizeof(double)))

	printf("</font>....... End initialization\n\n");
	return 1;
}