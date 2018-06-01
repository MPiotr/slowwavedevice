#include "multiplier.h"
#ifndef __MULTIPLIER_COMPLEXFIELD_H
#define __MULTIPLIER_COMPLEXFIELD_H
class Multiplier_ComplexField : public Multiplier
{
	double ElectronsDeltaEnergy(double _A);
	cplx ElectronCurrentA(double reA, double imA);
	cplx HfElectronCurrent(double  _reB, double _imB, double _A);

	double DeltaEnergyComplexField(double _A);
	cplx CurrentAComplexField(double _reA, double _imA);	//�� ��� ��� ���������� ��������� ������������ �����
	cplx CurrentBComplexField(double  _reB, double _imB, double _A);	//�� ��� ��� ���������� ��������� ������������ �����
	
public:
//	Multiplier_ComplexField() : Device() { ; }
	Multiplier_ComplexField(QDomDocument * doc) : Multiplier(doc) { ; }
	void loadField(char *filename);
};
#endif