#include "orotron.h"

/*Orotron::Orotron() :Multiplier()
{
	;
}*/
Orotron::Orotron(QDomDocument *doc) : Multiplier(doc)
{
	;
}
double Orotron::ElectronsDeltaEnergy(double A)
{
	return DeltaEnergy(A);
}
void Orotron::iterate(int paramsInd)
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

			printf("periods =%g, voltage = %g, current = %g\n", Nperiods, voltage, Current_ampers);

			double kpd = solveOrotron(Current_ampers, period, Nperiods, Ld, Lb, k1, Norma, Norma1b, voltage, 0, delta_freq, Q, Qb, wall, NULL, NULL, 0);
			printResultOrotron(results, kpd);
		}
	}
}
void Orotron::orotronStartCurrent()
{
	if (!initialized)
		printf("Error:  Not initialized\n"); return;

	double startCurr = orotronStartCurrent(Current_ampers, period, Nperiods, Ld, Lb, k1, Norma, Norma1b, voltage, 0, delta_freq, Q, Qb, wall, NULL, NULL, 0);
	printf("calculated start current = %g\n", startCurr);


}
double Orotron::orotronStartCurrent(double _Current_ampers, double _period, int _Nperiods, double _Ld, double _Lb, double _k1, double _Norma, double Norma1B, double _voltage,
double _inputPower_watts, double _delta, double _Qa, double _Qb, double _wall, char *filename, char *comment, double a00)
{

	Current_ampers = _Current_ampers;   Ld = _Ld; Lb = _Lb;
	period = _period;				  	k1 = _k1; Norma = _Norma;
	Nperiods = _Nperiods;			  	voltage = _voltage;
	delta_freq = _delta;				Q = _Qa; Qb = _Qb;
	inputPower_watts = _inputPower_watts;
	Norma = _Norma*0.5*double(Nperiods);
	Norma1b = Norma1B;
	wall = _wall;


	if (filename && comment && A_stat == 0) PrintParams(filename, comment);

	int N_it = 15;

	double G = paramG(); //k1 в мм^{-1}, Norma в мм^3.
	double omega = k1*(10.*c);
	double h = 2.*Pi / period, d = period;
	g1 = sqrt(h*h - k1*k1); g3 = 3 * g1;

	//	if(!solverStatus) 	solverStatus = initMultiplierSolver(3000, 10.1, 0., Ns);


	double nextA = 0.001;
	double A0 = 1e-9;
	double enDelta = ElectronsDeltaEnergy(A0);

	return  -(A0*A0 / (enDelta)) / (2 * Q*G);



}
void Orotron::printCurrentParamsOrotron(FILE *file)
{
	fprintf(file, "%g,%g,%g,%g,", Nperiods, Current_ampers, voltage, Q);
}
void Orotron::printDataHeaderOrotron(FILE *file)
{
	fprintf(file, "колво периодов,Ток[A],Hапряжение[кВ],добротность (периода),power,efficiency\n");
}
void Orotron::printResultOrotron(FILE *file, double kpd)
{
	fprintf(file, "%g,%g\n", power(kpd), kpd);
}
void Orotron::solveOrotron()
{
	if (!initialized)
		printf("Error:  Not initialized\n"); return;

	printDataHeaderOrotron(results);
	if (iteratedParams.size() > 0)
		iterate(0);
	else
	{
		printf("voltage = %g, Nperiod = %i, beta_gr = %g\n", voltage, Nperiods, group_speed);
		fflush(0);

		double kpd = solveOrotron(Current_ampers, period, Nperiods, Ld, Lb, k1, Norma, Norma1b, voltage, 0, delta_freq, Q, Qb, wall, NULL, NULL, 0);
		printResultOrotron(results, kpd);

	}
	return;
}
double Orotron::solveOrotron(double _Current_ampers, double _period, int _Nperiods, double _Ld, double _Lb, double _k1, double _Norma, double Norma1B, double _voltage,
double _inputPower_watts, double _delta, double _Qa, double _Qb, double _wall, char *filename, char *comment, double a00)
{

	Current_ampers = _Current_ampers;   Ld = _Ld; Lb = _Lb;
	period = _period;				  	k1 = _k1; Norma = _Norma;
	Nperiods = _Nperiods;			  	voltage = _voltage;
	delta_freq = _delta;				Q = _Qa; Qb = _Qb;
	inputPower_watts = _inputPower_watts;
	Norma = _Norma*0.5*double(Nperiods);
	Norma1b = Norma1B;
	wall = _wall;


	if (filename && comment && A_stat == 0) PrintParams(filename, comment);

	int N_it = 15;

	double G = paramG(); //k1 в мм^{-1}, Norma в мм^3.
	double omega = k1*(10.*c);
	double h = 2.*Pi / period, d = period;
	g1 = sqrt(h*h - k1*k1); g3 = 3 * g1;

	//	if(!solverStatus) 	solverStatus = initMultiplierSolver(3000, 10.1, 0., Ns);

	double nextA = 0.001;
	double A0 = 1e-6;
	if (A0*A0 + 2.*Q*G*ElectronsDeltaEnergy(A0) > 0) { printf("Start conditions are not satisfied\n"); return 0; }//Тест на стартовый ток (входной мощности нет) F > 0 -> старта нет
	double nextB = 0;
	if (a00 == 0)
		if (A_stat > 0.0001) { nextA = A_stat; N_it = 10; }
		else
			if (a00 > 0.00001) { nextA = a00; N_it = 10; }
	A_stat = findAstat(nextA, N_it, G);

	double deltaGamma = ElectronsDeltaEnergy(A_stat);

	double gamma = 1. + voltage / 511.;
	return deltaGamma / (gamma - 1);


}

