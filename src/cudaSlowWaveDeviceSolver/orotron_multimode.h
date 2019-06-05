#ifndef __OROTRON_MULTIMODE_H
#define __OROTRON_MULTIMODE_H
#include "multiplier_multimode.h"
class OrotronMutlimode : public MultiplierMultiModes
{
public:
	OrotronMutlimode() : MultiplierMultiModes() { ; }
	OrotronMutlimode(QDomDocument *);
	double solveMultiModeOrotron(double _Current_ampers, double _period, int _Nperiods, double _Ld, double _Lb, double _k1, double _Norma, double _NormaB, double _voltage, double _inputPower_watts, double _delta, double _Qa, double _Qb, double _wall, double Time, int _Namm, double initAmpMax, char* outputAmps, char *filename, char *comment, double a00 = 0);
};


#endif  __OROTRON_MULTIMODE_H