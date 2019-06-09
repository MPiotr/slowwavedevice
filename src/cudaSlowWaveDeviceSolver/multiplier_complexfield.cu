#include "multiplier_complexfield.h"
#include "cu_mult.h"

#define PAR ParamsM
#ifndef dm_Pi
__device__ const double dm_Pi = 3.141592653589793;
#endif

extern int Nz;
extern double Lsolver;

__global__ void
__launch_bounds__(512, 2)
MotionEquationMultiplierComplexField(PAR par, double Lstop, int Nharm, double A, double2 B)//Fixed Structure
{

	unsigned int p0 = threadIdx.x;     unsigned int Np = blockDim.x;
	unsigned int q0 = threadIdx.y;	   unsigned int Nq = blockDim.y;
	unsigned int s0 = threadIdx.z;	   unsigned int Ns = blockDim.z;

	unsigned int q_init = Nq*blockIdx.x + q0;  unsigned int Nq_max = Nq*gridDim.x;
	unsigned int s_init = Ns*blockIdx.y + s0;  unsigned int Ns_max = Ns*gridDim.y;
	unsigned int v_init = blockIdx.z;		   unsigned int Nv_max = gridDim.z;	

	//	printf("Thread %i/%i, %i/%i started;    ", q0, Nq, p0, Np);

	double la, lb, ld, h, k1, delta, voltage, g1, g3, *tAr, *tAi;

	__shared__ double avEN, int_rJ3, int_iJ3;

	//	printf("Step0;    ");
	int N;

	double dz;

	N = par.Nz;


	la = par.la; lb = par.lb; ld = par.ld;
	h = par.h; k1 = par.k1; delta = par.delta;
	g1 = par.g1; g3 = par.g3;
	tAr = par.tAr;  tAi = par.tAi;
	voltage = par.voltage;

	double *rJ3 = par.rJ3;
	double *iJ3 = par.iJ3;

	//	double *int_rJ3 = par.int_rJ3;
	//	double *int_iJ3 = par.int_iJ3;

	dz = par.L / (double)N;

	double z;
	int ifinal = floor(Lstop / dz);

	double Q, Qk1, Qk2, Qk3, Qk4;
	double W, Wk1, Wk2, Wk3, Wk4;
	double fAr, fAi, fB, rAp, rApp, iAp, iApp, r;
	double ifnotdestroyed = 1;

	double R_cyclotron, R_center, kappa_cyclotron, phase_cyclotron, initial_angle, tmpPhCycl;
	double wall = par.wall;

	R_center = 0.5*wall + wall*((double)q_init - 0.5*(double)Nq_max) / (double)(Nq_max);
	initial_angle = (0.0810194 - 2.05972*R_center + 28.0433*R_center*R_center);
	R_cyclotron = (0.568*initial_angle + 0.035156*((double)v_init) / double(Nv_max));
	kappa_cyclotron = 1.758;
	phase_cyclotron = 2.*dm_Pi*(double)s_init / (double)Ns_max;

	double en0 = 1. + voltage / 511.;
	en0 -= 0.5*initial_angle*initial_angle*(en0*en0 - 1)*en0;

	Q = 2.*dm_Pi / double(Np)*double(p0);// + 1./(double)Nq*((double)q0 + (double)s0/(double)Ns));
	W = 0;


	__shared__ double    shQ[NS][NQ][NP];
	__shared__ double    shW[NS][NQ][NP];
	__shared__ double d2_rJ3[NQ][NP];
	__shared__ double d2_iJ3[NQ][NP];
	__shared__ double d1_rJ3[NP];
	__shared__ double d1_iJ3[NP];

	double PH, EN, cosPH, sinPH, cosPS, sinPS, rB, iB;

	double H = h;//+dh(delta);

	if (p0 + q_init + s_init + v_init == 0)
	{
		rJ3[0] = 0;
		iJ3[0] = 0;
	}

	if (p0 + q0 + s0 == 0)
	{
		par.int_rJ3[gridDim.x*gridDim.y*blockIdx.z + gridDim.x*blockIdx.y + blockIdx.x] = 0;
		par.int_iJ3[gridDim.x*gridDim.y*blockIdx.z + gridDim.x*blockIdx.y + blockIdx.x] = 0;
		int_rJ3 = 0;
		int_iJ3 = 0;

		avEN = 0;
	}


	int NAstop = floor(la / dz) - 1;
	//	if(q0 + p0 + s0 == 0) printf("NAstop: %i\n" , NAstop);
	int i = 0;  double rdr;
	for (i = 1; i < N; i++)
	{


		/////////////////
		z = (double)i*dz;
		r = (R_center + R_cyclotron*cos(kappa_cyclotron*z + phase_cyclotron)); rdr = (R_center + R_cyclotron*cos(kappa_cyclotron*(z + dz) + phase_cyclotron));
		ifnotdestroyed *= (r > -wall) ? 1. : 0.;
		rAp = (i < NAstop) ? tAr[i] * exp(-g1*r) : 0; rApp = (i < NAstop) ? tAr[i + 1] * exp(-g1*r) : 0;
		iAp = (i < NAstop) ? tAi[i] * exp(-g1*r) : 0; iApp = (i < NAstop) ? tAi[i + 1] * exp(-g1*r) : 0;
		//	rAp =  (z < la)? sin(dm_Pi/la*z)*exp(-g1*r) : 0; 
		//	rApp = (z < la)? sin(dm_Pi/la*(z+dz))*exp(-g1*rdr) : 0;
		iAp = 0; iApp = 0;

		//	if((q0 + p0 + s0 == 0)&&(i < 10)) printf("kern: %g, %g, %g, %g\n", rAp, iAp, rApp, iApp);
		//	ifnotdestroyed = 1;

		PH = Q;
		EN = W + en0;

		fAr = rAp; fAi = iAp;
		fB = (((z>la + ld) && (z < la + ld + lb)) ? sin(dm_Pi / lb*(z - la - ld))*exp(-g3*r) : 0);

		rB = B.x*fB;
		iB = B.y*fB;

		sincos(PH, &sinPH, &cosPH);
		sincos(3.*PH, &sinPS, &cosPS);

		Qk1 = dz*(H - k1*EN / sqrt(EN*EN - 1.));
		Wk1 = -dz*(A*(fAr*cosPH - fAi*sinPH) + (rB*cosPS - iB*sinPS))*ifnotdestroyed;

		/////////////////
		z = ((double)i + 0.5)*dz;
		r = (R_center + R_cyclotron*cos(kappa_cyclotron*z + phase_cyclotron));
		ifnotdestroyed *= (r > -wall) ? 1. : 0.;
		//	ifnotdestroyed = 1;

		PH = Q + 0.5*Qk1;
		EN = W + 0.5*Wk1 + en0;

		fAr += rApp; fAr *= 0.5;
		fAi += iApp; fAi *= 0.5;
		fB = ((z>la + ld) && (z < la + ld + lb)) ? sin(dm_Pi / lb*(z - la - ld))*exp(-g3*r) : 0;
		rB = B.x*fB;
		iB = B.y*fB;

		sincos(PH, &sinPH, &cosPH);
		sincos(3.*PH, &sinPS, &cosPS);

		Qk2 = dz*(H - k1*EN / sqrt(EN*EN - 1.));
		Wk2 = -dz*(A*(fAr*cosPH - fAi*sinPH) + (rB*cosPS - iB*sinPS))*ifnotdestroyed;
		/////////////////
		PH = Q + 0.5*Qk2;
		EN = W + 0.5*Wk2 + en0;

		sincos(PH, &sinPH, &cosPH);
		sincos(3.*PH, &sinPS, &cosPS);

		Qk3 = dz*(H - k1*EN / sqrt(EN*EN - 1.));
		Wk3 = -dz*(A*(fAr*cosPH - fAi*sinPH) + (rB*cosPS - iB*sinPS))*ifnotdestroyed;
		/////////////////	
		z = ((double)i + 1)*dz;
		r = (R_center + R_cyclotron*cos(kappa_cyclotron*z + phase_cyclotron));
		ifnotdestroyed *= (r > -wall) ? 1. : 0.;
		//		ifnotdestroyed = 1;

		PH = Q + Qk3;
		EN = W + Wk3 + en0;

		fAr -= 0.5*rAp; fAr *= 2.;
		fAi -= 0.5*iAp; fAi *= 2.;
		fB = (((z>la + ld) && (z < la + ld + lb)) ? sin(dm_Pi / lb*(z - la - ld))*exp(-g3*r) : 0);
		rB = B.x*fB;
		iB = B.y*fB;

		sincos(PH, &sinPH, &cosPH);
		sincos(3.*PH, &sinPS, &cosPS);

		Qk4 = dz*(H - k1*EN / sqrt(EN*EN - 1.));
		Wk4 = -dz*(A*(fAr*cosPH - fAi*sinPH) + (rB*cosPS - iB*sinPS))*ifnotdestroyed;
		///////////////		

		Q += 1. / 6.*(Qk1 + 2.*Qk2 + 2.*Qk3 + Qk4);

		W += 1. / 6.*(Wk1 + 2.*Wk2 + 2.*Wk3 + Wk4);


		//		printf("#< %i -> (%i, %i, %i)>\t", i, s0, q0, p0);
		//		printf("#<%i>", Np_Q*Nq*(Ns*i+s0) + Np_Q*q0 + p0);

		//		if(q0+p0+s0 == 0) printf("%i", i);

		shQ[s0][q0][p0] = Q;

		shW[s0][q0][p0] = W;



		__threadfence();
		__syncthreads();

		//////// усреднение какой-то гармоники
		if (s0 == 0)
		{
			double tmp_rJ3 = 0, tmp_iJ3 = 0;

			for (int ii = 0; ii < Ns; ii++)
			{
				int ii_init = Ns*blockIdx.y + ii;
				tmpPhCycl = 2.*dm_Pi*(double)ii_init / (double)Ns_max;
				r = (R_center + R_cyclotron*cos(kappa_cyclotron*z + tmpPhCycl));





				if (Nharm == 1)
				{
					double tm_exp = exp(-g1*r);			//rAp + I iAp используется как множитель при интегрировании тока вдоль продольной координаты ВНИМАНИЕ!! fB зависит от q0, s0
					rAp = (i < NAstop) ? tAr[i] * tm_exp : 0;
					iAp = (i < NAstop) ? tAi[i] * tm_exp : 0;

				}
				else
				{
					rAp = (((z>la + ld) && (z < la + ld + lb)) ? sin(dm_Pi / lb*(z - la - ld)) : 0)*exp(-g3*r);
					iAp = 0;
				}

				fB *= ifnotdestroyed;


				PH = shQ[ii][q0][p0];
				sincos((double)Nharm*PH, &sinPS, &cosPS);

				tmp_rJ3 += (cosPS*rAp + sinPS*iAp)*ifnotdestroyed;
				tmp_iJ3 += (-sinPS*rAp + cosPS*iAp)*ifnotdestroyed;
			}
			d2_rJ3[q0][p0] = tmp_rJ3;
			d2_iJ3[q0][p0] = tmp_iJ3;


		}



		__threadfence();
		__syncthreads();
		if (s0 + q0 == 0)
		{
			double tmp_rJ3 = 0, tmp_iJ3 = 0;

			for (int ii = 0; ii < Nq; ii++)
			{
				tmp_rJ3 += d2_rJ3[ii][p0];
				tmp_iJ3 += d2_iJ3[ii][p0];
			}
			d1_rJ3[p0] = tmp_rJ3;
			d1_iJ3[p0] = tmp_iJ3;
		}
		__threadfence();
		__syncthreads();
		if (p0 + q0 + s0 == 0)
		{
			double tmp_rJ3 = 0, tmp_iJ3 = 0;

			for (int ii = 0; ii < Np; ii++)
			{
				tmp_rJ3 += d1_rJ3[ii];
				tmp_iJ3 += d1_iJ3[ii];
			}
			rJ3[i] = tmp_rJ3;
			iJ3[i] = tmp_iJ3;
			if (Nharm == 3)
			{
				int_rJ3 += tmp_rJ3*fB;
				int_iJ3 += tmp_iJ3*fB;
			}
			else
			{
				int_rJ3 += tmp_rJ3*rAp - tmp_iJ3*iAp;
				int_iJ3 += tmp_rJ3*iAp + tmp_iJ3*rAp;
			}
		}
		__threadfence();
		__syncthreads();
		//////////////////// конец усреднения какой-то гармоники

		//		if((q0+p0 == 0)&&(s0 == 0)) printf("%i\t%g\t%g\n", q0, PH, EN);

		//		if(q0+p0+s0 == 0) printf("....%i\t", i);


		if (i == ifinal)
		{
			if (s0 == 0)
			{
				double tmp_W = 0;

				for (int ii = 0; ii < Ns; ii++)
				{
					EN = shW[ii][q0][p0];
					tmp_W += EN;
				}
				d2_rJ3[q0][p0] = tmp_W;
			}
			__threadfence();
			__syncthreads();
			if (s0 + q0 == 0)
			{

				double tmp_rJ3 = 0;
				for (int ii = 0; ii < Nq; ii++)
					tmp_rJ3 += d2_rJ3[ii][p0];
				d1_rJ3[p0] = tmp_rJ3;

			}
			__threadfence();
			__syncthreads();
			if (p0 + q0 + s0 == 0)
			{
				double tmp_rJ3 = 0;
				for (int ii = 0; ii < Np; ii++)
					tmp_rJ3 += d1_rJ3[ii];
				(avEN) += tmp_rJ3;
				//				if(i < 320)  printf("%g, %g\n", z, tmp_rJ3);


			}
		}
		__threadfence();
		__syncthreads();






		if (i > ifinal) break;


	}
	//	printf("END\t");

	//	if(p0 + q0 +s0 == 0) printf("%g, %g>__<%g, %g>...", rJ3[i-1]/double(16*16*2), iJ3[i-1]/double(16*16*2), int_rJ3/double(16*16*2*N)*z, int_iJ3/double(16*16*2*N)*z);

	__syncthreads();
	//	if(p0+q_init+s_init + v_init  == 0)
	if (p0 + q0 + s0 == 0)
	{
		/*	printf("%i, %i, %i\t (%g, %g)\n",
		blockIdx.x, blockIdx.y, blockIdx.z,
		int_rJ3/double(Np*Nq_max*Ns_max*N)*(par.L),
		int_iJ3/double(Np*Nq_max*Ns_max*N)*(par.L)) ;*/
		par.avEN[gridDim.x*gridDim.y*blockIdx.z + gridDim.x*blockIdx.y + blockIdx.x] = avEN;
		par.int_rJ3[gridDim.x*gridDim.y*blockIdx.z + gridDim.x*blockIdx.y + blockIdx.x] = int_rJ3;
		par.int_iJ3[gridDim.x*gridDim.y*blockIdx.z + gridDim.x*blockIdx.y + blockIdx.x] = int_iJ3;

	}

}

void Multiplier_ComplexField::loadField(char *filename)
{
	double dz = Lmax/double(Nmax);
	double La = (double)Nperiods*period;
	int Nstop = ceil(La/dz) + 1;

	FILE *file;
	fopen_s(&file, filename, "r");
//	fseek(file, 0L, SEEK_SET );

	float X[6000];
	float fieldRe[6000];
	float fieldIm[6000];

	tAr = new double [Nstop];
	tAi = new double [Nstop];

	float x1 = 1, x2 = 1, x3 = 1;

	
	int q = 0;
	int res;
	do
	{
		res = fscanf_s(file, "%g %g %g", &x1, &x2, &x3);
		X[q] = x1; fieldRe[q] = x2; fieldIm[q] = x3;


		q++;
	
	}while((q < 6000)&&(res > 0));
	


	int s = 0;
	double x;
	for(int i = 0; i < Nstop; i++)
	{
		x = dz * i;
		for(int j = s; j < q - 2; j++)
		{
			if((X[j] < x) && (X[j + 1] > x)) 
			{
				tAr[i]  = fieldRe[j] + (fieldRe[j+1] - fieldRe[j])/(X[j+1] - X[j])*(x - X[j]);
				tAi[i]  = fieldIm[j] + (fieldIm[j+1] - fieldIm[j])/(X[j+1] - X[j])*(x - X[j]);
				s = j;
				break;
			}
			if(X[j] == x)
			{
				tAr[i] = fieldRe[j];
				tAi[i] = fieldIm[j];
				s = j;
				break;
			}
		}
		
	}
	tAr[Nstop - 1] = 0; 	tAi[Nstop - 1] = 0;

	fclose(file);

	double absMax = -1;
	int iMax = 0;
	for(int i = 0; i < Nstop; i++)
	{
		if(  (tAr[i]*tAr[i] + tAi[i]*tAi[i]) > absMax) {  absMax = (tAr[i]*tAr[i] + tAi[i]*tAi[i]) ; iMax = i;}
	}
	absMax = sqrt(tAr[iMax]*tAr[iMax] + tAi[iMax]*tAi[iMax]);
	for(int i = 0; i < Nstop; i++) {tAr[i] /= absMax; tAi[i] /= absMax;}
	double integral = 0;
	for(int i = 0; i < Nstop; i++) {integral +=  (tAr[i]*tAr[i] + tAi[i]*tAi[i])*dz;}


	Norma = Norma*integral/period;


//debug...
	FILE *debfile = fopen("F:\\Piotr\\FemProjects\\mult_proj\\debug.csv", "w");
	for(int i = 0; i < Nstop; i++)
	{
		fprintf(debfile, "%g,%g,%g\n", dz*i, tAr[i], tAi[i]);
	//	fprintf(debfile, "%g,%g,%g\n", X[i], fieldRe[i], fieldIm[i]);
	}
	fclose(debfile);
	//...debug
	
	printf("d_tAr alloc: %i;\t",cudaMalloc((void**) &d_tAr, Nstop*sizeof(double)));
	printf("d_tAi alloc: %i;\t",cudaMalloc((void**) &d_tAi, Nstop*sizeof(double)));

	printf("d_tAr copy: %i;\t",cudaMemcpy((void*) d_tAr, tAr,  Nstop*sizeof(double), cudaMemcpyHostToDevice));
	printf("d_tAi copy: %i;\t",cudaMemcpy((void*) d_tAi, tAi,  Nstop*sizeof(double), cudaMemcpyHostToDevice));

/*	double *debAr = new double [100];
	cudaMemcpy(debAr, d_tAr, 99*sizeof(double), cudaMemcpyDeviceToHost);
	for(int j = 0; j < 100; j++) printf("deb: %g\n", debAr[j]);
	delete [] debAr;*/

	fieldLoaded = true;
	
}
double Multiplier_ComplexField::ElectronsDeltaEnergy(double A)
{
	if (fieldLoaded) return DeltaEnergyComplexField(A);
	else return 0;
}
cplx Multiplier_ComplexField::HfElectronCurrent(double  rB, double iB, double Astat)
{
	if (fieldLoaded) return CurrentBComplexField(rB, iB, Astat);
	else return 0;
}
cplx Multiplier_ComplexField::ElectronCurrentA(double reA, double imA)
{
	if (fieldLoaded) return CurrentAComplexField(reA, imA);
	else return 0;
}
std::complex<double> Multiplier_ComplexField::CurrentBComplexField(double  reB, double imB, double A)
{
	PAR par;
	double d = period;
	double h = 2.*Pi / d;
	double La = period*double(Nperiods);


	par.tAr = d_tAr;  par.tAi = d_tAi;
	par.la = La; par.lb = Lb; par.ld = Ld; par.k1 = k1; par.h = h; par.voltage = voltage;
	par.Nz = Nz; par.L = Lsolver; par.wall = wall;
	par.g1 = g1; par.g3 = g3;

	//	printf("CurrentB: %g, %g, %g \n", La, Ld, Lb);

	par.rJ3 = d_rJ3;  par.iJ3 = d_iJ3;

	par.avEN = d_avEN;
	par.int_rJ3 = d_int_rJ3;
	par.int_iJ3 = d_int_iJ3;

	double dz = (La + Lb + Ld) / double(Nz);


	dim3 threadsPerBlock(NP, NQ, NS);


	double2 B; B.x = reB; B.y = imB;

	//	printf("\n B loop: %g\n", La+Ld+Lb );
	//	printf("\n Threads: %i, %i, %i\n", threadsPerBlock.x, threadsPerBlock.y, threadsPerBlock.z );

	dim3 numblocks = (Nq / NQ, Ns / NS, Nv);

	//	cudaMemcpy(	   d_rJ3, &dzero, sizeof(double), cudaMemcpyHostToDevice);
	//	cudaMemcpy(	   d_iJ3, &dzero, sizeof(double), cudaMemcpyHostToDevice);
	//	cudaMemcpy(d_int_rJ3, &dzero, sizeof(double), cudaMemcpyHostToDevice);
	//	cudaMemcpy(d_int_iJ3, &dzero, sizeof(double), cudaMemcpyHostToDevice);
	//	cudaMemcpy(   d_avEN, &dzero, sizeof(double), cudaMemcpyHostToDevice);



	MotionEquationMultiplierComplexField << <dim3((size_t)Nq / NQ, (size_t)Ns / NS, (size_t)Nv), threadsPerBlock >> >(par, La + Ld + Lb, 3, A, B);


	return retriveBCurr();



}
std::complex<double> Multiplier_ComplexField::CurrentAComplexField(double  reA, double imA)
{
	PAR par;

	double d = period;
	double h = 2.*Pi / d;
	double La = period*double(Nperiods);

	par.la = La; par.lb = 1.; par.ld = 1.; par.k1 = k1; par.h = h; par.voltage = voltage;
	par.Nz = Nz; par.L = Lsolver; par.wall = wall;
	par.g1 = g1; par.g3 = g3;


	par.rJ3 = d_rJ3;  par.iJ3 = d_iJ3;
	par.tAr = d_tAr;  par.tAi = d_tAi;

	par.avEN = d_avEN;
	par.int_rJ3 = d_int_rJ3;
	par.int_iJ3 = d_int_iJ3;

	double2 zero = { 0, 0 };

	dim3 threadsPerBlock(NP, NQ, NS);

	double A; A = sqrt(reA*reA + imA*imA);

	//	printf("\n B loop: %g\n", La+Ld+Lb );
	//	printf("\n Threads: %i, %i, %i\n", threadsPerBlock.x, threadsPerBlock.y, threadsPerBlock.z );

	dim3 numblocks = (Nq / NQ, Ns / NS, Nv);

	//	cudaMemcpy(	   d_rJ3, &dzero, sizeof(double), cudaMemcpyHostToDevice);
	//	cudaMemcpy(	   d_iJ3, &dzero, sizeof(double), cudaMemcpyHostToDevice);
	//	cudaMemcpy(d_int_rJ3, &dzero, sizeof(double), cudaMemcpyHostToDevice);
	//	cudaMemcpy(d_int_iJ3, &dzero, sizeof(double), cudaMemcpyHostToDevice);
	//	cudaMemcpy(   d_avEN, &dzero, sizeof(double), cudaMemcpyHostToDevice);



	MotionEquationMultiplierComplexField << <dim3((size_t)Nq / NQ, (size_t)Ns / NS, (size_t)Nv), threadsPerBlock >> >(par, La, 1, A, zero);


	return retriveBCurr()*exp(I*arg(reA + I*imA));



}
double Multiplier_ComplexField::DeltaEnergyComplexField(double A)
{
	PAR par;
	double	t_deltaEn = 0;
	double d = period;
	double h = 2.*Pi / d;
	double La = period*double(Nperiods);

	par.la = La; par.lb = Lb; par.ld = Ld; par.k1 = k1; par.h = h; par.voltage = voltage;
	par.Nz = Nz; par.L = Lsolver; par.wall = wall;
	par.g1 = g1; par.g3 = g3;

	par.rJ3 = d_rJ3;  par.iJ3 = d_iJ3;
	par.tAr = d_tAr;  par.tAi = d_tAi;





	par.avEN = d_avEN;
	par.int_rJ3 = d_int_rJ3;
	par.int_iJ3 = d_int_iJ3;




	double2 zero = { 0, 0 };

	//	cudaMemcpy(	   d_rJ3, &dzero, sizeof(double), cudaMemcpyHostToDevice);
	//	cudaMemcpy(	   d_iJ3, &dzero, sizeof(double), cudaMemcpyHostToDevice);
	//	cudaMemcpy(d_int_rJ3, &dzero, sizeof(double), cudaMemcpyHostToDevice);
	//	cudaMemcpy(d_int_iJ3, &dzero, sizeof(double), cudaMemcpyHostToDevice);
	//	cudaMemcpy(   d_avEN, &dzero, sizeof(double), cudaMemcpyHostToDevice);




	dim3 threadsPerBlock(NP, NQ, NS);

	MotionEquationMultiplierComplexField << <dim3((size_t)Nq / NQ, (size_t)Ns / NS, (size_t)Nv), threadsPerBlock >> >(par, La, 1, A, zero);

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


