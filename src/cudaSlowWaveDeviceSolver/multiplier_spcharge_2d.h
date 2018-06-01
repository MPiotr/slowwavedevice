#include "multiplier.h"
#ifndef __MULTIPLIER_SPCHARGE_2D_H
#define __MULTIPLIER_SPCHARGE_2D_H
class Multiplier_SpCharge_2D : public Multiplier
{
	double ElectronsDeltaEnergy(double _A);
	cplx ElectronCurrentA(double reA, double imA);
	cplx HfElectronCurrent(double  _reB, double _imB, double _A);

	double DeltaEnergySpaceCharge2D(double A);
	cplx CurrentASpaceCharge2D(double rB, double iB);
	cplx CurrentBSpaceCharge2D(double  _reB, double _imB, double _A, double enPrint = 0);
public:
//	Multiplier_SpCharge_2D() : Device() { ; }
	Multiplier_SpCharge_2D(QDomDocument * doc) : Multiplier(doc) { ; }
	void  PrintCurrentsSpaceCharge2D(char *filename1, char *filename2);
	bool initMultiplierSolver(int nz, double lsolver, double groupSpeedCoeff, char *_solverName);
};
#endif