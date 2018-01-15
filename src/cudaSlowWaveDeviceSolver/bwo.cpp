#include <direct.h>
#include "bwo.h"
#include "twt_1d.h"
#include "twt_2d.h"

#include "xml_routines.h"

template class BWO<TWT_1D>;
template class BWO<TWT_2D>;

template  <class TWTsolver>
BWO<TWTsolver>::BWO(QDomDocument *doc) : TWTsolver(doc)
{
	delta = 0;
	a0 = 1e-9;
	parasiteIndex = 0;
	if (!setXMLEntry(doc, "initialAmp", &a0)) a0 = 1e-9;
	if (!setXMLEntry(doc, "initialDelta", &delta)) delta = -0.01;
	if (!setXMLEntry(doc, "boundaryReflection", &reflectionCoeff)) reflectionCoeff = 0;
	if (!setXMLEntry(doc, "boundaryReflectionPhase", &reflectionPhase)) reflectionPhase = 0;

}
template  <class TWTsolver>
void BWO<TWTsolver>::iterate(int paramsInd)
{
	for (vec::iterator s = iteratedParams.at(paramsInd).begin(); s != iteratedParams.at(paramsInd).end(); ++s)
	{
		if (paramsInd != iteratedParams.size() - 1)
		{
			changeParam(iteratedNames[paramsInd], *s);
			iterate(paramsInd + 1);
		}
		else
		{
			changeParam(iteratedNames[paramsInd], *s);
			if (syncwave){
				synch_angle = syncwave->wavenumber();
				group_speed = syncwave->groupSpeed();
				k1 = syncwave->frequency() * 2 * Pi / 299.8;
			}

			printf("periods =%g, voltage = %g, synch_angle = %g\n", Nperiods, voltage, synch_angle);

			double h = 2.*Pi / period*(synch_angle / 360.);
			double G = paramG(h);
			double lossK = lossKappa(h);
			solveBWO(a0, delta);
			printResults(results);
		}
	}
}

template  <class TWTsolver>
void BWO<TWTsolver>::printResults(FILE *file)
{
	double eff;

	double h = 2 * Pi / period*(synch_angle / 360.);
	double La = Nperiods*period*h;

	double PoutLeft = -power(a0, h);
	double PoutRight = -power(a0*exp(lossKappa(h)*La)*(1 - reflectionCoeff), h);
	double shiftedFre = (1. + delta*(k1 / h) / fabs(group_speed))*299.8*k1 / (2.*PI);
	eff = efficiency(a0, h);

	printCurrentParams(file);
	fprintf(file, "%15g,%15g,%15g,%15g,%15g,%15g\n", a0, PoutLeft, PoutRight, delta, eff, shiftedFre);
}
template  <class TWTsolver>
void BWO<TWTsolver>::printResultsStartCurrent(FILE *file)
{
	double h = 2 * Pi / period*(synch_angle / 360.);
	printCurrentParams(file);
	double tmpdelta = delta;
	double Lst = getNormalizedLength(Current_ampers, &tmpdelta);
	fprintf(file, "%15g,%15g,%15g,%15g,%15g\n", Current_ampers, delta, paramG(h), Lst, tmpdelta);
}
template  <class TWTsolver> 
void BWO<TWTsolver>::printDataHeader(FILE *file)
{
	printParamsHeader(file);
	printResultHeader(file);
}
template  <class TWTsolver> 
void BWO<TWTsolver>::printParamsHeader(FILE *file)
{
	fprintf(file, "кол-во периодов,угол[градусы],Ток[A],Hапряжение[кВ],групповая скорость,добротность (периода),входная мощность[Вт],");
}
template  <class TWTsolver> 
void BWO<TWTsolver>::printResultHeader(FILE *file)
{
	fprintf(file, "output amp[W], output power (cathode)[W], output power (collector) [W], delta, efficiency, shifted frequency\n");
}
template  <class TWTsolver> 
void BWO<TWTsolver>::printStCurrResultHeader(FILE* file)
{
	fprintf(file, "start current [A], delta, G start, normalized length, normalized delta\n");
}

template  <class TWTsolver>  
void BWO<TWTsolver>::solveBWO()
{
	printDataHeader(results);
	if (iteratedParams.size() > 0)
	{
		iterate(0);
	}
	else
	{
		solveBWO(a0, delta);
		printResults(results);
	}

}
template  <class TWTsolver>  
void BWO<TWTsolver>::solveBWO(double A0, double Delta0)
{
	if (syncwave){
		synch_angle = syncwave->wavenumber(voltage);
		double fre = syncwave->frequency(voltage);
		group_speed = syncwave->groupSpeed(voltage);
		k1 = 2 * Pi*fre / 299.8;
	}

	double h = 2.*Pi / period*(synch_angle / 360.);
	double G = paramG(h);
	double lossKappa = k1 / h*0.5 / (Q*group_speed);
	a0 = A0;
	delta = Delta0;

	printf("calculating BWO:\n\t freq. = %g, wavenum. = %g deg., voltage = %g kV, current = %g A, reflection = %g,\n\tbeta_gr = %g,  G = %g\n\n", 299.8*k1 / (2.*Pi), synch_angle, voltage, Current_ampers, reflectionCoeff, group_speed, G);


	double La = Nperiods*period*h;
	double la = Nperiods*period*h;
	double dz = Lmax / double(Nmax);
	int Nstop = ceil(La / dz);
	if (Nstop > Nmax) Nstop = Nmax - 1;

	if (longitudinalStructureRe) generateLongitudinalStructure(h);

	cplx F00, F01, F10, F11;
	double J11, J12, J21, J22;
	double iJ11, iJ12, iJ21, iJ22;
	printf("%10s\t%10s\t%10s\t%10s\t%10s\n", "det", "delta", "a0", "real(F00)", "imag(F00)");
	double eps = 1e-11;
	double tol = 1e-10;

	cplx boundaryA;
	cplx reflection = reflectionCoeff*exp(I*reflectionPhase / 180. * Pi);
	for (int i = 0; i < 150; i++)
	{
		solve(A, ar, ai, a0, lossKappa, delta, structReal, structIm, G, 0, false, longitudinalStructureRe, longitudinalStructureIm, qStructure, mesh);
		boundaryA = a0*exp(-I*delta*La + lossKappa*La); //NB: lossKappa is negative
		F00 = A[Nstop - 1] / boundaryA - reflection;

		a0 *= 1. + eps;
		solve(A, ar, ai, a0, lossKappa, delta, structReal, structIm, G, 0, false, longitudinalStructureRe, longitudinalStructureIm, qStructure, mesh);
		boundaryA = a0*exp(-I*delta*La + lossKappa*La);
		F01 = A[Nstop - 1] / boundaryA - reflection;
		a0 /= 1. + eps;


		delta *= 1. + eps;
		solve(A, ar, ai, a0, lossKappa, delta, structReal, structIm, G, 0, false, longitudinalStructureRe, longitudinalStructureIm, qStructure, mesh);
		boundaryA = a0*exp(-I*delta*La + lossKappa*La);
		F10 = A[Nstop - 1] / boundaryA - reflection;
		delta /= 1. + eps;

		J11 = real(F10 - F00) / (delta*eps);   J12 = real(F01 - F00) / (a0*eps);
		J21 = imag(F10 - F00) / (delta*eps);   J22 = imag(F01 - F00) / (a0*eps);

		double det = J11*J22 - J21*J12;
		if (!isfinite(det)) { printf("iterating error: infinite determinant\n "); break; }
		if (det == 0)       { printf("iterating error: determinant is zero \n ");   break; }

		iJ11 = J22 / det;	iJ12 = -J12 / det;
		iJ21 = -J21 / det;  iJ22 = J11 / det;

		double ddelta = real(F00)*iJ11 + imag(F00)*iJ12;
		double dA = real(F00)*iJ21 + imag(F00)*iJ22;

		if (!isfinite(ddelta)) { printf("iterating error: infinite ddelta\n"); break; }
		if (!isfinite(dA))    { printf("iterating error: infinite da0\n"); break; }

		delta -= ddelta;
		a0 -= dA;

		if (shMemoryCreated) printAbsAtoSharedMemory(Nstop);
		if (fabs(ddelta / delta) < tol && fabs(ddelta / delta) < tol) break;

		printf("%10g\t%10g\t%10g\t%10g\t%10g\n", det, delta, a0, real(F00), imag(F00));
		fflush(0);


	}

	if (parasiteIndex == 0) solve(A, ar, ai, a0, lossKappa, delta, structReal, structIm, G, 0, false, longitudinalStructureRe, longitudinalStructureIm, qStructure, mesh);
	if (shMemoryCreated) printAbsAtoSharedMemory(Nstop);
	printf("power (output at cathode side) = %g\npower (output at collector side) =%g\ndelta = %g\nfrequency=%g\n", -power(a0, h), -power(a0*exp(lossKappa*La)*(1. - reflectionCoeff), h),
		delta,
		(1. + delta*(k1 / h) / fabs(group_speed))*299.8*k1 / (2.*PI));
}
template  <class TWTsolver>  
double BWO<TWTsolver>::currentFromNormalizedLength(double normLength, double *Delta0)
{
	if (syncwave){
		synch_angle = syncwave->wavenumber(voltage);
		double fre = syncwave->frequency(voltage);
		k1 = 2 * Pi*fre / 299.8;
		group_speed = syncwave->groupSpeed(voltage);
	}

	double h = 2.*Pi / period*(synch_angle / 360.);
	double gamma = 1. + voltage / 511.;
	double beta = sqrt(gamma*gamma - 1) / gamma;
	double nu = (k1 / h) / pow(gamma*gamma - 1, 1.5);
	double L = Nperiods*h*period;
	double norm1 = 2 * Norma / double(Nperiods);

	double res = pow(normLength / L, 3)*2. / (nu / norm1)*(Current_ampers / (-paramG(h)));
	*Delta0 = *Delta0*normLength / L;

	return res;

}
template  <class TWTsolver>  
double BWO<TWTsolver>::getNormalizedLength(double curr, double *Delta0)
{
	double h = 2.*Pi / period*(synch_angle / 360.);
	double gamma = 1. + voltage / 511.;
	double beta = sqrt(gamma*gamma - 1) / gamma;
	double nu = (k1 / h) / pow(gamma*gamma - 1, 1.5);
	double L = Nperiods*h*period;
	double norm1 = 2 * Norma / double(Nperiods);
	if (norm1 < 0) printf("Norma is negative, error calculating dimensionless values. Probably, synch. field is set by profile function\n");
	double normLength = pow(paramG(h)/ Current_ampers*curr*nu / norm1*0.5, 1. / 3.)*L; // .5 под степенью чтобы совпадало с чёрной книжкой (там делиться не на 2пи, а на пи)
	//	double normLength = pow(-curr*nu*0.5, 1. / 3.)*L; // .5 под степенью чтобы совпадало с чёрной книжкой (там делиться не на 2пи, а на пи)


	printf("k1 = %g\nh = %g\ngamma = %g\nbeta = %g\n\nnu=%g\nL = %g\ndelta = %g\n", k1, h, gamma, beta, nu, normLength, *Delta0 * L / normLength);

	*Delta0 = *Delta0 * L / normLength;
	return normLength;
}
template  <class TWTsolver>  
double BWO<TWTsolver>::findStartCurrent()
{
	double res = findStartCurrent(Current_ampers, &delta);
	Current_ampers = res;
	double tmpdelta = delta;
	double h = 2.*Pi / period*(synch_angle / 360.);
	double Lst = getNormalizedLength(res, &tmpdelta);
	printf("\nCalculated start current:   %g\n\t start delta:     %g\n\t start G:     %g\n\tdim-less start length     %g\n\tnormalized delta     %g\n",
		res, delta, res*paramG(h) / Current_ampers, Lst, tmpdelta); fflush(0);

	char tmp[300]; _getcwd(tmp, 300);
	printf(tmp);
	printf("\n");
	FILE *startCurrResultFile = fopen("startCurrResult.dat", "w");
	printParamsHeader(startCurrResultFile);
	printStCurrResultHeader(startCurrResultFile);
	printResultsStartCurrent(startCurrResultFile);
	fclose(startCurrResultFile);
	return res;
}
template  <class TWTsolver>  
double BWO<TWTsolver>::findStartCurrent(double I0, double *Delta0)
{
	if (syncwave){
		synch_angle = syncwave->wavenumber(voltage);
		double fre = syncwave->frequency(voltage);
		k1 = 2 * Pi*fre / 299.8;
		group_speed = syncwave->groupSpeed(voltage);
	}


	double h = 2.*Pi / period*(synch_angle / 360.);
	double G = paramG(h)*I0 / Current_ampers;
	double lossKappa = k1 / h*0.5 / (Q*group_speed);
	a0 = 1e-10;
	delta = *Delta0;

	if (longitudinalStructureRe) generateLongitudinalStructure(h);

	printf("calculatring BWO start curent: \n\t frequency = %g, wavenumber = %g degree,\n\t group beta = %g, voltage = %gkV, reflection = %g\n\n", 299.8*k1 / (2.*Pi), synch_angle, group_speed, voltage, reflectionCoeff);

	double La = Nperiods*period*h;
	double la = Nperiods*period*h;
	double dz = Lmax / double(Nmax);
	int Nstop = ceil(La / dz);
	if (Nstop > Nmax) Nstop = Nmax - 1;

	cplx F00, F01, F10, F11;
	double J11, J12, J21, J22;
	double iJ11, iJ12, iJ21, iJ22;
	if (debug)	printf("Newton iterations\n%10s\t%10s\t%10s\t%10s\t%10s\n", "det", "delta", "G", "real(F00)", "imag(F00)");
	double eps = 1e-10;
	double tol = 1e-10;
	cplx boundaryA;
	cplx reflection = reflectionCoeff*exp(I*reflectionPhase / 180. * Pi);
	for (int i = 0; i < 150; i++)
	{
		solve(A, ar, ai, a0, lossKappa, delta, structReal, structIm, G, 0, false, longitudinalStructureRe, longitudinalStructureIm, qStructure, mesh);
		boundaryA = a0*exp(-I*delta*La - lossKappa*La);
		F00 = A[Nstop - 1] / boundaryA - reflection;


		G *= 1. + eps;
		solve(A, ar, ai, a0, lossKappa, delta, structReal, structIm, G, 0, false, longitudinalStructureRe, longitudinalStructureIm,  qStructure, mesh);
		boundaryA = a0*exp(-I*delta*La - lossKappa*La);
		F01 = A[Nstop - 1] / boundaryA - reflection;
		G /= 1. + eps;

		delta *= 1. + eps;
		solve(A, ar, ai, fabs(a0), lossKappa, delta, structReal, structIm, G, 0, false, longitudinalStructureRe, longitudinalStructureIm, qStructure, mesh);
		boundaryA = a0*exp(-I*delta*La - lossKappa*La);
		F10 = A[Nstop - 1] / boundaryA - reflection;
		delta /= 1. + eps;

		J11 = real(F10 - F00) / (delta*eps);   J12 = real(F01 - F00) / (G*eps);
		J21 = imag(F10 - F00) / (delta*eps);   J22 = imag(F01 - F00) / (G*eps);

		double det = J11*J22 - J21*J12;
		if (det == 0) break;

		iJ11 = J22 / det;	iJ12 = -J12 / det;
		iJ21 = -J21 / det;  iJ22 = J11 / det;

		double ddelta = real(F00)*iJ11 + imag(F00)*iJ12;
		double dG = real(F00)*iJ21 + imag(F00)*iJ22;


		if (!isfinite(ddelta)) { printf("iterating error: infinite ddelta\n"); break; }
		if (!isfinite(dG))    { printf("iterating error: infinite dG\n"); break; }

		delta -= ddelta;
		G -= dG;

		if (shMemoryCreated) printAbsAtoSharedMemory(Nstop);
		if (fabs(ddelta / delta) < tol && fabs(dG / G) < tol) break;
		if (debug) printf("%10g\t%10g\t%10g\t%10g\t%10g\n", det, delta, G, real(F00), imag(F00));
		fflush(0);
	}
	printf("used amplitude:    %g\n", abs(A[0]));
	if (parasiteIndex == 0) solve(A, ar, ai, fabs(a0), lossKappa, delta, structReal, structIm, G, 0, true, longitudinalStructureRe, longitudinalStructureIm, qStructure, mesh);
	*Delta0 = delta;
	return G / (paramG(h))*Current_ampers;

}
//@obsolete
template  <class TWTsolver>  
double BWO<TWTsolver>::findStartCurrent(double I0, double Delta0, double fre, double synch_angle, double gr_speed, double Q)
{
	double h = 2.*Pi / period*(synch_angle / 360.);
	k1 = 2.*Pi*fre / 299.8;
	group_speed = gr_speed;

	double G = paramG(h)*I0 / Current_ampers;
	double inputAmp = amplitude(inputPower_watts, h);
	double lossKappa = k1 / h*0.5 / (Q*group_speed);
	a0 = 1e-10;
	delta = Delta0;


	double La = Nperiods*period*h;
	double la = Nperiods*period*h;
	double dz = Lmax / double(Nmax);
	int Nstop = ceil(La / dz);
	if (Nstop > Nmax) Nstop = Nmax - 1;

	cplx F00, F01, F10, F11;
	double J11, J12, J21, J22;
	double iJ11, iJ12, iJ21, iJ22;
	if (debug) printf("%10s\t%10s\t%10s\t%10s\t%10s\n", "det", "delta", "G", "real(F00)", "imag(F00)");
	double eps = 1e-7;
	cplx boundaryA;
	cplx reflection = reflectionCoeff*exp(I*reflectionPhase / 180. * Pi);
	for (int i = 0; i < 150; i++)
	{
		solve(A, ar, ai, a0, lossKappa, delta, structReal, structIm, G, 0, false, longitudinalStructureRe, longitudinalStructureIm, qStructure, mesh);
		boundaryA = a0*exp(-I*delta*La - lossKappa*La);
		F00 = A[Nstop - 1] / boundaryA - reflection;

		G *= 1. + eps;
		solve(A, ar, ai, a0, lossKappa, delta, structReal, structIm, G, 0, false, longitudinalStructureRe, longitudinalStructureIm, qStructure, mesh);
		boundaryA = a0*exp(-I*delta*La - lossKappa*La);
		F01 = A[Nstop - 1] / boundaryA - reflection;
		G /= 1. + eps;

		delta *= 1. + eps;
		solve(A, ar, ai, fabs(a0), lossKappa, delta, structReal, structIm, G, 0, false, longitudinalStructureRe, longitudinalStructureIm, qStructure, mesh);
		boundaryA = a0*exp(-I*delta*La - lossKappa*La);
		F10 = A[Nstop - 1] / boundaryA - reflection;
		delta /= 1. + eps;

		J11 = real(F10 - F00) / (delta*eps);   J12 = real(F01 - F00) / (G*eps);
		J21 = imag(F10 - F00) / (delta*eps);   J22 = imag(F01 - F00) / (G*eps);

		double det = J11*J22 - J21*J12;
		if (det == 0) break;

		iJ11 = J22 / det;	iJ12 = -J12 / det;
		iJ21 = -J21 / det;  iJ22 = J11 / det;

		double ddelta = real(F00)*iJ11 + imag(F00)*iJ12;
		double dG = real(F00)*iJ21 + imag(F00)*iJ22;


		if (!isfinite(ddelta)) { printf("iterating error: infinite ddelta\n"); break; }
		if (!isfinite(dG))    { printf("iterating error: infinite dG\n"); break; }

		delta -= ddelta;
		G -= dG;

		if (fabs(ddelta / delta) < 1e-7 && fabs(dG / G) < 1e-7) break;
		if (debug) printf("%10g\t%10g\t%10g\t%10g\t%10g\n", det, delta, G, real(F00), imag(F00));


	}

	if (parasiteIndex == 0) solve(A, ar, ai, fabs(a0), lossKappa, delta, structReal, structIm, G, 0, true, longitudinalStructureRe, longitudinalStructureIm, qStructure, mesh);
	return G / (-paramG(h))*Current_ampers;

}
