#ifndef SYNCHRONISM_H
#define SYNCHRONISM_H
#include <vector>
#include "Interpolation.h"

class Synchronism
{
	int Nmodes;
	double _voltage, _frequency, _phaseSpeed, _d, _groupSpeed, _inputFreq;

	std::vector< Interpolation * > modes;
	std::vector<double> frequencies;
	std::vector<double> wavenumbers;
	
	void loadModes(char *modesFile);
	double findRootFixedFre(Interpolation* mode, double beta, double h0);
	double findRootFixedBeta(Interpolation* mode, double beta, double h0);

	
public:
	Synchronism(char *modesFile);
	Synchronism(char *modesFile, double d, double inputFre);
	Synchronism();

	void addMode(char *modesFile);
	double frequency();
	double frequency(double voltage);
	double frequency(double voltage, int n);
	double groupSpeed();
	double groupSpeed(double voltage);
	double groupSpeed(double voltage, int n);
	double wavenumber();
	double wavenumber(double voltage);
	double wavenumber(double voltage, int n);
	double voltage();

	void setVoltage(double voltage);
	void setFrequency(double fre, double h0 = 0);
	
	void unitTesting();
};


#endif