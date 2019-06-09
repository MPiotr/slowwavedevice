#include "orotron_multimode.h"
OrotronMutlimode::OrotronMutlimode(QDomDocument *doc) : MultiplierMultiModes(doc)
{
	;
}

double OrotronMutlimode::solveMultiModeOrotron(double _Current_ampers, double _period, int _Nperiods, double _Ld, double _Lb, double _k1, double _Norma, double Norma1B, double _voltage,
	double _inputPower_watts, double _delta, double _Qa, double _Qb, double _wall, double Time, int _Namm, double initAmpMax, char*outputAmps, char *filename, char *comment, double a00)
{

	Current_ampers = _Current_ampers;   Ld = _Ld; Lb = _Lb;
	period = _period;				  	k1 = _k1; Norma = _Norma;
	Nperiods = _Nperiods;			  	voltage = _voltage;
	delta_freq = _delta;				Q = _Qa; Qb = _Qb;
	inputPower_watts = _inputPower_watts;
	Norma = _Norma*0.5*double(Nperiods);
	Norma1b = Norma1B;
	wall = _wall;
	Namm = _Namm;
	dz = Lmax / (double)Nmax;

	sprintf(difractionFlag, "Empiric");
	cplx *Amps = new cplx[Namm];

	if (A_stat == 0) PrintParams(filename, comment);

	int N_it = 15;

	double G = paramG(); //k1 β μμ^{-1}, Norma β μμ^3.
	double omega = k1*(10.*c);
	double h = 2.*Pi / period, d = period;
	g1 = sqrt(h*h - k1*k1); g3 = 3 * g1;



	//	if(!solverStatus) 	solverStatus = initMultiplierSolver(3000, 10.1, 0., Ns);

	int GQ = Nq / NQ; int GS = Ns / NS; int GV = Nv;
	cplx *kAmps = new cplx[Namm * 5];
	double *buffRe = new double[GQ*GS*GV*Nmax];
	double *buffIm = new double[GQ*GS*GV*Nmax];
	cplx *tmpAmps = new cplx[Namm];
	cplx *currs = new cplx[Namm];

	for (int i = 0; i <= Namm / 2; i++) Amps[i] = initAmpMax*cplx(double(rand()) / double(RAND_MAX), double(rand()) / double(RAND_MAX));
	for (int n = Namm / 2 + 1; n<Namm; n++) Amps[n] = Amps[Namm - n - 1];


	//...... RK method 4-4
	/*debug*/ G = 9e-4*h*h;
	/*debug*/ Q = 16;
	/*debug*/// Amps[0] = 0.1e-10*h; Amps[1] = 1e-10*h; Amps[2] = 0.5e-10*h;

	double coeffs[4] = { 0, 0.5, 0.5, 1 };
	int Niter = 10000;
	double dT = Time / double(Niter);

	FILE* ampsFile = fopen(outputAmps, "w");

	for (int iter = 0; iter < Niter; iter++)
	{
		fprintf(ampsFile, "%g,", dT*iter);
		for (int am = 0; am < Namm; am++) fprintf(ampsFile, "%g,%g,%g,%g,", Amps[am].real(), Amps[am].imag(), abs(Amps[am]) / h, arg(Amps[am]));
		fprintf(ampsFile, "\n");

		for (int k = 0; k < 4; k++)
		{
			for (int n = 0; n <= Namm / 2; n++)
				tmpAmps[n] = Amps[n] + coeffs[k] * kAmps[k*Namm + n];
			for (int n = Namm / 2 + 1; n<Namm; n++)
				tmpAmps[n] = tmpAmps[Namm - n - 1];


			CurrentAMultiModes(tmpAmps, currs, buffRe, buffIm);

			for (int n = 0; n <= Namm / 2; n++)
				kAmps[(k + 1)*Namm + n] = dT*(G*(currs[n] + currs[Namm - n - 1]) - 0.5*(tmpAmps[n]) / Q);
			for (int n = Namm / 2 + 1; n<Namm; n++)
				kAmps[(k + 1)*Namm + n] = kAmps[(k + 1)*Namm + Namm - n - 1];

		}

		for (int n = 0; n < Namm; n++)
			Amps[n] += (kAmps[Namm + n] + 2.*kAmps[2 * Namm + n] + 2.*kAmps[3 * Namm + n] + kAmps[4 * Namm + n]) / 6.;

		double deltaK = G*imag(currs[1] * exp(-I*arg(Amps[1]))) / abs(Amps[1]);
		if (fabs(deltaK) < 0.1) k1 -= deltaK*k1;



	}
	//...... End RK method 

	delete[] buffRe, buffIm, kAmps, tmpAmps, currs, Amps;
	return 0;


}