#ifndef __MULTIPLIER_THREE_CAVITY_H
#define __MULTIPLIER_THREE_CAVITY_H
#include "multiplier.h"
class MultiplierTreeCavity : public Multiplier {
protected:
	double Lb;							    //����� �� ������ (�������)
	double Qb;								//��������� ����������� �� ������

	double La2, Ld1, Ld2;				    //����� �� ������������ � ��������� ������ � ��������������� �����
	double Norma1, Norma2;					//����� �� ������ (�������) � ��������������� �����
	double Norma1b;							//����� �� ������ (���� ������)
	double Qa1, Qa2;						//����������� ������ �� � ������ �� ������ (��������������� �����)
	
	void PrintParamsDoubleScheme(char *filename, char*comment);
	cplx findAstatDoubleScheme(double nextA, int Nit, double G1, double G2);  //���������� ������������ ��������� �� ������ �� ������, ��������� � ������ ������������� A_stat
	cplx CurrentA2(double A1, double  _reA, double _imA);		//� ������� �����, ��� �� ������ ������?
	
	cplx FindBstatDoubleScheme(cplx b0, int Nb_it, double Gb, double Astat, cplx A2_stat);
	cplx CurrentB2(double  _reB, double _imB, double _A, cplx _A2stat); //�� ��� � ������� �����

	double DeltaEnergyDoubleScheme();
public:
	MultiplierTreeCavity(QDomDocument* doc);
	double getHFoutputPowerDoubleScheme(double _Current_ampers, double _period, int _Nperiods, double _Ld1, double _Ld2, double _La2, double _Lb, double _k1, double _Norma, double Norma1B, double _voltage, double _inputPower_watts, double _delta, double _Qa1, double _Qa2, double _Qb, double _wall, char *filename, char *comment);
	double retriveBPower(cplx b, double omega, double Qb, double NormaB);
};
#endif