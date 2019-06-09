#ifndef __MULTIPLIER_MULTIMODE_H
#define __MULTIPLIER_MULTIMODE_H
#include "multiplier.h"
class MultiplierMultiModes : public Multiplier {
protected:
	int Namm; // количество  продольных мод для мультимодового солвера
	double Ld1, Lb, Ld;
	double Qa1, Qa2, Qb;
	double Norma1, Norma1b;

	void PrintParamsMultiModes(char *filename, char *comment);

	void findAstatMultiModes(double* nextX, int N_it, double G1, double G2); //переписывает Х

	cplx CurrentAMultiModes(double *X, int Na = 3, cplx *J1 = 0, cplx * J2 = 0);				//НЧ ток (интеграл) 
	void CurrentAMultiModes(std::complex<double> *Amps, std::complex<double> * currs, double *buffRe, double *buffIm, int Na = 3, cplx *J1 = 0, cplx *J2 = 0); // Для многомодового решателя перегруженный
	void retriveACurrComplex(std::complex<double>  *Amps, std::complex<double>  *currs, double *currsBuffRe, double *currsBuffIm, int Na, int Nstop);

public:
	double getHFoutputPowerMultiModes(double _Current_ampers, double _period, int _Nperiods, double _Ld, double _Lb, double _k1, double _Norma, double Norma2, double _voltage, double _inputPower_watts, double _delta, double _Qa1, double _Qa2, double _Qb, double _wall, char *filename, char *comment);
	MultiplierMultiModes();
	MultiplierMultiModes(QDomDocument *doc);

};

#endif