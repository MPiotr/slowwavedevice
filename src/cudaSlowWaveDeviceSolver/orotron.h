#include "device.h"
#ifndef OROTRON__H
#define OROTRON__H
class Orotron : public Device
{
	void iterate(int);
	void printResultOrotron(FILE *file, double kpd);
	void printCurrentParamsOrotron(FILE *file);
	void printDataHeaderOrotron(FILE *file);
	double ElectronsDeltaEnergy(double A);
public: 
	Orotron();
	Orotron(QDomDocument *doc);
	void solveOrotron();
	double solveOrotron(double _Current_ampers, double _period, int _Nperiods, double _Ld, double _Lb, double _k1, double _Norma, double _NormaB, double _voltage, double _inputPower_watts, double _delta, double _Qa, double _Qb, double _wall, char *filename, char *comment, double a00 = 0);

	void orotronStartCurrent();
	double orotronStartCurrent(double _Current_ampers, double _period, int _Nperiods, double _Ld, double _Lb, double _k1, double _Norma, double _NormaB, double _voltage, double _inputPower_watts, double _delta, double _Qa, double _Qb, double _wall, char *filename, char *comment, double a00 = 0);

	double solveMultiModeOrotron(double _Current_ampers, double _period, int _Nperiods, double _Ld, double _Lb, double _k1, double _Norma, double _NormaB, double _voltage, double _inputPower_watts, double _delta, double _Qa, double _Qb, double _wall, double Time, int _Namm, double initAmpMax, char* outputAmps, char *filename, char *comment, double a00 = 0);
};

#endif