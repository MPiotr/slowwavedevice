#include "cu_mult.h"
#include "twt_1d.h"

#define PAR ParamsM
#ifndef dm_Pi
__device__ const double dm_Pi = 3.141592653589793;
#endif



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
__device__ void findPosition(double *mesh, double x, double* xplus, double* xminus, int *nplus, int *nminus, int N)
{
	if (x >= mesh[N - 1])
	{
		*xplus = mesh[N - 1];
		*xminus = mesh[N - 2];
		*nplus = N - 1;
		*nminus = N - 2;
		return;
	}
	else if (x <= mesh[0])
	{
		*xplus = mesh[1];
		*xminus = mesh[0];
		*nplus = 1;
		*nminus = 0;
		return;
	}
	else
	{
		int nmax = N-1;
		int nmin = 0;
		for (int i = 0; i < round(log2((double)N)); i++)
		{
			if (mesh[(nmax - nmin) / 2 + nmin] == x) { nmax = (nmax - nmin) / 2 + nmin; nmin = nmax - 1; break; }
			if (mesh[(nmax - nmin) / 2 + nmin] > x) nmax = (nmax - nmin) / 2 + nmin;
			else nmin = (nmax - nmin) / 2 + nmin;
			if (nmax - nmin == 1) break;
		}
		*nplus = nmax;
		*nminus = nmin;
		*xplus = mesh[nmax];
		*xminus = mesh[nmin];
	}

}

__device__ double fieldAtX(double *F, double *mesh, double x, int Nmesh)
{
	double xplus, xminus;
	int nplus, nminus;

	findPosition(mesh, x, &xplus, &xminus, &nplus, &nminus, Nmesh);

	double coefMinus = (xplus - x) / (mesh[nplus] - mesh[nminus]);
	double coefPlus = (x - xminus) / (mesh[nplus] - mesh[nminus]);
	return coefPlus*F[nplus] + coefMinus*F[nminus];
}

__global__ void particleStep(PAR *par, double *Ar, double *Ai, double *mesh, int NumMesh, int i, int k, double K, double coupling = 1)
{
	unsigned int p0 = threadIdx.x;			unsigned int Np = blockDim.x;    //x coordinates of beam x-section, ph - cyclotron phase spread
	unsigned int x0 = threadIdx.y;			unsigned int Nx = blockDim.y;    //p -- wave phase
	unsigned int ph0 = threadIdx.z;		    unsigned int Nph = blockDim.z;

	unsigned int P  = blockIdx.x;		    unsigned int NP_ = gridDim.x;  //P-grid index is for v_trans spread
	unsigned int X  = blockIdx.y;			unsigned int NX =  gridDim.y;
	unsigned int PC = blockIdx.z;			unsigned int NPC = gridDim.z;

//	unsigned int p = Np*P + p0;	  	    unsigned int Np_max = Np*NP_;
	unsigned int x = Nx*X + x0;		    unsigned int Nx_max = Nx*NX;
	unsigned int ph = NPC*PC + ph0;	//	unsigned int Nph_max = Nph*NPC; 

	//index ph warps all over spreads, for now: cyclotron phase, transversal velocity.
	unsigned int NcyclPhase = par->Na1;
	unsigned int NvTrans = NP_;
	unsigned int ph_cyclPhase = ph;
	unsigned int ph_v_trans = P;

	unsigned int warpsize = Nx*Nph*Np;
	unsigned int log2warpsize = round(log2((double)warpsize));
	int xi = (p0 + Np*x0 + ph0*Np*Nx);           //развернутый индекс
	int gstride = (NX*NPC*P + NX*PC + X)*warpsize;

	//...... load problem parameters ......................... 
	double h, k1, voltage, v_trans_max;
	h = par->h;
	k1 = par->k1;
	voltage = par->voltage;
	v_trans_max = par->v_trans_max;
	// cyclotron parameters
	double h_cycl = par->h_cycl; double omega_cycl = par->omega_cycl; // циклотронное волновое число и частота
	double ph_cyl = 2 * dm_Pi / (double) NcyclPhase*ph_cyclPhase;
	double v_trans = v_trans_max / (double)NvTrans*ph_v_trans;
	double R_cycl = (omega_cycl != 0) ? v_trans * 299.8e9 / omega_cycl : 0; //R_cycl is in mm, so 2.998e10 -> 299.8e9 
	// mesh parameters
	double beamThick = par->beamThickX;
    //TODO: make an beam position in input. So far, zero of mesh, corresponds to the beam center
	double r0 = (int(x) - int(Nx_max/ 2) )*beamThick / (double)Nx_max;  //initial positions of the particles
	int N = par->Nz;
	double dz = par->L / (double)N;
	double en0 = sqrt(pow(1. + voltage / 511., 2) - v_trans*v_trans);
	// .... get acess to values arrays..................
	double *d_Qk = par->Qk;   double *d_Wk = par->Wk;   //k - компоненты фазы и энергии
	double *Q0 = par->Q0;   double *W0 = par->W0;   //фаза and энергия
	double rA = par->ar0[i];  double iA = par->ai0[i];  //амплитуда (огибающая по z) 	
		// ..... init shared memory ...........................

	__shared__ double    sh_imJ[NS*NQ*NP];
	__shared__ double    sh_reJ[NS*NQ*NP];
	__shared__ double    sh_mesh[NMESHMAX];
	__shared__ double    sh_Ar[NMESHMAX];
	__shared__ double    sh_Ai[NMESHMAX];

	double Q, Qk;
	double W, Wk;
	double PH, EN, z, sinPH, cosPH;

	// ...... init/load phase and energy ............	
	if (i != 0) //load
	{
		Q = Q0[gstride + xi];
		W = W0[gstride + xi];
	}
	else  // or init
	{
		Q = 2.*dm_Pi / double(Np)*double(p0);
		W = 0;
		Q0[gstride + xi] = Q;
		W0[gstride + xi] = W;

	}
	if (xi < NumMesh)
	{
		sh_mesh[xi] = mesh[xi];
		sh_Ar[xi] = Ar[xi];
		sh_Ai[xi] = Ai[xi];
	}
	// ....... at the left end of the space, set HF current to zero
	if ((i == 0) && (xi == 0))
	{
		par->avEN[NX*PC + X] = 0;
		par->int_rJ3[NX*PC + X] = 0;
		par->int_iJ3[NX*PC + X] = 0;
	}
	__syncthreads();
	// .......init k coefficient of phase and energy 

	Wk = (k > 0) ? d_Wk[4 * gstride + warpsize*(k - 1) + xi] : 0;
	Qk = (k > 0) ? d_Qk[4 * gstride + warpsize*(k - 1) + xi] : 0;

	// --------- end of initialization phase ----------

	// --------- make a step -------------------------

	z = ((double)i + K)*dz;
	
	double r = R_cycl*cos(h_cycl*z + ph_cyl) + r0;
	double frA = fieldAtX(sh_Ar, sh_mesh, r, NumMesh);
	double fiA = fieldAtX(sh_Ai, sh_mesh, r, NumMesh);

	frA *= coupling;
	fiA *= coupling;

	PH = Q + K*Qk;
	EN = W + K*Wk + en0;

	//	if((i < 2)&&(p0 == 0)&&(q0 < 2)) printf("%i,%i,%i\t%g\t%g\t%g\t%g\n",i, k, q0, rQ1, iQ1,  iq0, iak);

	sincos(PH, &sinPH, &cosPH);


	//...........Уравнения движения определены здесь...............................
	double DQ = (1. - (k1 / h)*EN / sqrt(EN*EN - 1.));
	double reF = frA*cosPH - fiA*sinPH;
	double imF = frA*sinPH + fiA*cosPH;
	Qk = dz*DQ;
	Wk = -dz*(rA*reF - iA*imF);
	//....................................................................

	d_Qk[4 * gstride + warpsize*k + xi] = Qk;
	d_Wk[4 * gstride + warpsize*k + xi] = Wk;

	// ===========Усреднение=======================	

	sh_imJ[xi] = -imF / double(warpsize);
	sh_reJ[xi] = reF / double(warpsize);

	__syncthreads();


	biAverage(sh_imJ, sh_reJ, xi, warpsize, log2warpsize); // averaged values are stored in zero indexes

	//========== Запись=====================

	if (xi == 0)
	{
		par->int_iJ3[NX*NPC*P + NX*PC + X] = sh_imJ[0];
		par->int_rJ3[NX*NPC*P + NX*PC + X] = sh_reJ[0];
	}
	__syncthreads();

}

static __global__ void amplitudeStep(PAR *par, int i, int k, int prevK, double kappa)
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
	double dz = par->L / (double)N;
	//	double rak, iak;


	__shared__ double   shJr[SHARRAYSIZE];
	__shared__ double   shJi[SHARRAYSIZE];

	shJr[xi] = rJ[xi] / warpsize; //чтобы не делить после усреднения
	shJi[xi] = iJ[xi] / warpsize;


	//	unsigned int stride = warpsize;

	__syncthreads();
	//....... усреднение по блоку (по варпу усреднено на пред. шаге)................
	biAverage(shJr, shJi, xi, warpsize, log2warpsize);
	//..............................................................................

	//.............. уравнения на амплитуду определены здесь .........................
	if (xi == 0)
	{
		par->rJ3[0] = shJr[0];
		par->iJ3[0] = shJi[0];
		double arr = par->ar0[i];
		double aii = par->ai0[i];
		if (k > 0)
		{
			rAk[k] = dz*(G*shJr[0] - kappa*(arr + prevK*rAk[k - 1]) + delta*(aii + prevK*iAk[k - 1]));
			iAk[k] = dz*(G*shJi[0] - kappa*(aii + prevK*iAk[k - 1]) - delta*(arr + prevK*rAk[k - 1]));
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

static __global__ void endstep(PAR *par, int i)
{
	// make an addition to the phase, energy, and amlitude in the end of one Runge-Kutta (4,4) step.
	unsigned int p0 = threadIdx.x;			unsigned int Np = blockDim.x;
	unsigned int x0 = threadIdx.y;			unsigned int Nx = blockDim.y;
	unsigned int y0 = threadIdx.z;		    unsigned int Ny = blockDim.z;

	unsigned int P = blockIdx.x;	  //unsigned int NP_ = gridDim.x;
	unsigned int X = blockIdx.y;		unsigned int NX = gridDim.y;
	unsigned int Y = blockIdx.z;		unsigned int NY = gridDim.z;
	
/*	unsigned int p = Np*P + p0;		    unsigned int Np_max = Np*NP_;
	unsigned int x = Nx*X + x0;		    unsigned int Nx_max = Nx*NX;
	unsigned int y = Ny*Y + y0;			unsigned int Ny_max = Ny*NY;*/

	unsigned int warpsize = Nx*Ny*Np;
	unsigned int log2warpsize = round(log2((double)warpsize));
	int xi = (p0 + Np*x0 + y0*Np*Nx);
//	unsigned int gridsize = NX*NY*NP_;
	unsigned int XI = (NX*NY*P + NX*Y + X);

	// ............ init problem parameters ...................................
	//	unsigned int xi = threadIdx.x;	    unsigned int Np = blockDim.x;
	//	unsigned int warpsize = Np;

	// .................... get access to phase, energy and current  k-coefficients
	double *d_Qk = par->Qk;		double *Q0 = par->Q0;
	double *d_Wk = par->Wk;		double *W0 = par->W0;
	double *rAk = par->rAk;    double *iAk = par->iAk;   //k - компоненты амплитуды

	unsigned int stride = 4 * XI*warpsize;
	unsigned int gstride = XI*warpsize;

	// ............... calculate next phase and energy .....................

	Q0[gstride + xi] += (d_Qk[stride + xi] + 2.*d_Qk[stride + xi + warpsize] + 2.*d_Qk[stride + xi + 2 * warpsize] + d_Qk[stride + xi + 3 * warpsize]) / 6.;
	W0[gstride + xi] += (d_Wk[stride + xi] + 2.*d_Wk[stride + xi + warpsize] + 2.*d_Wk[stride + xi + 2 * warpsize] + d_Wk[stride + xi + 3 * warpsize]) / 6.;


	// ............... calculate next amplitude ............................	
	if (X + xi == 0)
	{
		par->ar0[i + 1] = par->ar0[i] + (rAk[0] + 2.*rAk[1] + 2.*rAk[2] + rAk[3]) / 6.;
		par->ai0[i + 1] = par->ai0[i] + (iAk[0] + 2.*iAk[1] + 2.*iAk[2] + iAk[3]) / 6.;
	}

}

std::complex<double>  TWT_1D::solveTWT_1d(std::complex<double>  *A, double *Ar, double *Ai, double inputAmp, double lossKappa, double delta,
	double *fieldStructureRe, double *fieldStructureIm, double *mesh, double G, double enPrint, bool printField, double *lStrRe, double *lStrIm, double *qStr)
{
	int GX = Nq / NQ; int GY = Ns / NS; int GP = Nv;

	int warpsize = NP*NQ*NS;
	int gridsize = GP*GX*GY;

	double K[4] = { 0, 0.5, 0.5, 1. };

	dim3 motionwarp = dim3(NP, NQ, NS);
	dim3 motiongrid = dim3(GP, GX, GY);

	if (firstRun) {
		printf("warp setup:\n %i (wave phase) x %i (space spread) x  %i other (cycl. phase) = %i\n", NP, NQ, NS, warpsize);
		printf("grid setup:\n %i (trans. velocity) x %i (space spread) x %i (cycl. phase) = %i\n", GP, GX, GY, gridsize);
		printf("total number of particles\n, %i (wave phase), %i (space spread), %i (cyclotron spread), %i (trans_velocities)\n", NP, NQ*GX /*=Nq*/, NS*GY /*=Ns*/, GP);
		firstRun = false;
	}


	gpuErrChk(cudaMemset(d_Qk, 0, 4 * sizeof(double)*Np*Nq*Ns*Nv));	gpuErrChk(cudaMemset(d_W0, 0, sizeof(double)*Np*Nq*Ns*Nv));
	gpuErrChk(cudaMemset(d_Wk, 0, 4 * sizeof(double)*Np*Nq*Ns*Nv));	gpuErrChk(cudaMemset(d_rAk, 0, sizeof(double) * 4 * gridsize));
	gpuErrChk(cudaMemset(d_Q0, 0, sizeof(double)*Np*Nq*Ns*Nv));	    gpuErrChk(cudaMemset(d_iAk, 0, sizeof(double) * 4 * gridsize));

	gpuErrChk(cudaMemset(d_ai0, 0, sizeof(double)*Nmax));
	gpuErrChk(cudaMemset(d_ar0, 0, sizeof(double)*Nmax));
	gpuErrChk(cudaMemset(d_avEN, 0, sizeof(double)));

	gpuErrChk(cudaMemcpy(d_fAr, fieldStructureRe, NumMesh*sizeof(double), cudaMemcpyHostToDevice));
	gpuErrChk(cudaMemcpy(d_fAi, fieldStructureIm, NumMesh*sizeof(double), cudaMemcpyHostToDevice));
	gpuErrChk(cudaMemcpy(d_mesh, mesh, NumMesh*sizeof(double), cudaMemcpyHostToDevice));

	double h = 2 * Pi / period*(synch_angle / 360.);
	double La = Nperiods*period*h;
	double dz = Lmax / double(Nmax);

//TODO make separate functions for copying twt_1d params
	ParamsM par = setPar();
	par.G = G;
	par.h = h;
	par.delta = delta;
	par.lossKappa = lossKappa;
	par.beamThickX = beamWidth;
	par.beamThickY = beamHeight;
	par.omega_cycl = omega_cycl;
	par.h_cycl = h_cycl / h;
	par.v_trans_max = v_trans_max;

//	double *energy = new double[Np*Nq*Ns*Nv];
//	double *phase = new double[Np*Nq*Ns*Nv];

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

	//TODO:  do something with filenames
	if (printField) {
		file = fopen("twt_debug.csv", "w");
		res_enr = fopen("enPhase.csv", "w");
	}

	if (NumMesh > NMESHMAX)
	{
		printf("Length of the field table is larger than NMESHMAX. Either change NMESHMAX and recompile of use sparser amplitude table\n");
		return 0;
	}
	if (gridsize > SHARRAYSIZE)
	{
		printf("The size of the grid is larger than the shared memory array used for averaging\n");
		//TODO get rid of this limitations : add grid averaging for amplitudeStep
		return 0;
	}

	for (int i = 0; i < Nstop; i++)
	{

		for (int k = 0; k < 4; k++)
		{
			z = double(i)*dz + K[k] * dz;
			if (lStrRe) coupl = (1. - K[k])*lStrRe[i] + K[k] * lStrRe[i + 1]; else coupl = 1;
			if (qStr) loss = lossKappa / ((1. - K[k])*qStr[i] + K[k] * qStr[i + 1]); else loss = lossKappa;
			if (!isfinite(loss)) loss = lossKappa;
			// TODO: put d_mesh and NumMesh inside par
			particleStep << <motiongrid, motionwarp >> >(d_par, d_fAr, d_fAi, d_mesh, NumMesh, i, k, K[k], coupl);
			amplitudeStep << <1, gridsize >> >(d_par, i, k, K[k], loss);

	/*		double rJ, iJ;
			cudaMemcpy(&rJ, d_rJ3, sizeof(double), cudaMemcpyDeviceToHost);
			cudaMemcpy(&iJ, d_iJ3, sizeof(double), cudaMemcpyDeviceToHost);
		//	fprintf(res_enr, "%g,%g,%g\n", z, rJ, iJ);
			printf("%g,%g,%g\n", z, rJ, iJ);*/
		}
		endstep << <motiongrid, motionwarp >> >(d_par, i);

	}
	gpuErrChk(cudaMemcpy(Ar, d_ar0, sizeof(double)*Nstop, cudaMemcpyDeviceToHost));
	gpuErrChk(cudaMemcpy(Ai, d_ai0, sizeof(double)*Nstop, cudaMemcpyDeviceToHost));

	for (int i = 0; i < Nstop; i++)
	{
		z = double(i)*dz;
		A[i] = cplx(Ar[i], Ai[i]);
		if (printField)
			fprintf(file, "%g,%g,%g,%g\n", z, Ar[i], Ai[i], abs(A[i]));
	}


	if (printField){
		fclose(file);
		fclose(res_enr);
	}

	return cplx (Ar[Nstop - 1], Ai[Nstop - 1]);

}
bool TWT_1D::initSolver(int nz, double lsolver)
{
//	TWT::initSolver(nz, lsolver);

	gpuErrChk(cudaMalloc((void**)&d_fAr, NMESHMAX*sizeof(double)))
	gpuErrChk(cudaMalloc((void**)&d_fAi, NMESHMAX*sizeof(double)))
	gpuErrChk(cudaMalloc((void**)&d_mesh, NMESHMAX*sizeof(double)))

	printf("</font>....... End initialization\n\n");
	return 1;
}
	
