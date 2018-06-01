#include "multiplier.h"
#ifndef OROTRON__H
#define OROTRON__H
class Orotron : public Multiplier
{
	
	void printResultOrotron(FILE *file, double kpd);
	void printCurrentParamsOrotron(FILE *file);
	void printDataHeaderOrotron(FILE *file);
	virtual double ElectronsDeltaEnergy(double A);
	virtual void iterate(int);
public: 
//	Orotron();
	Orotron(QDomDocument *doc);
	void solveOrotron();
	double solveOrotron(double _Current_ampers, double _period, int _Nperiods, double _Ld, double _Lb, double _k1, double _Norma, double _NormaB, double _voltage, double _inputPower_watts, double _delta, double _Qa, double _Qb, double _wall, char *filename, char *comment, double a00 = 0);

	void orotronStartCurrent();
	double orotronStartCurrent(double _Current_ampers, double _period, int _Nperiods, double _Ld, double _Lb, double _k1, double _Norma, double _NormaB, double _voltage, double _inputPower_watts, double _delta, double _Qa, double _Qb, double _wall, char *filename, char *comment, double a00 = 0);
};

#endif