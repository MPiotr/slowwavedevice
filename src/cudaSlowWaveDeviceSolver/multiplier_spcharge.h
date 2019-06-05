#ifndef	__MULTIPLIER_SPCHARGE_H
#define __MULTIPLIER_SPCHARGE_H
#include "multiplier.h"
class Multipler_SpCharge : public Multiplier
{
	virtual double ElectronsDeltaEnergy(double _A);
	cplx ElectronCurrentA(double reA, double imA);
	double DeltaEnergySpaceCharge(double A);
	cplx HfElectronCurrent(double  _reB, double _imB, double _A);
	cplx CurrentBSpaceCharge(double  _reB, double _imB, double _A, double enPrint = 0);
public:
	Multipler_SpCharge();
	Multipler_SpCharge(QDomDocument *doc); 
	virtual bool initSolver(int nz, double lsolver, double groupSpeedCoeff, char *_solverName);

};
#endif