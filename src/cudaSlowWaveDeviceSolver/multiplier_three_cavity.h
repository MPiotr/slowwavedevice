#ifndef __MULTIPLIER_THREE_CAVITY_H
#define __MULTIPLIER_THREE_CAVITY_H
#include "multiplier.h"
class MultiplierTreeCavity : public Multiplier {
protected:
	double Lb;							    //длина ВЧ секции (третьей)
	double Qb;								//Омическая добротность ВЧ секции

	double La2, Ld1, Ld2;				    //Длины НЧ резонаторных и дрейфовых секций в трёхрезонаторной схеме
	double Norma1, Norma2;					//норма нч секций (целиком) в трёхрезонаторной схеме
	double Norma1b;							//норма вч секции (один период)
	double Qa1, Qa2;						//Добротность первой НЧ и второй НЧ секций (трёхрезонаторная схема)
	
	void PrintParamsDoubleScheme(char *filename, char*comment);
	cplx findAstatDoubleScheme(double nextA, int Nit, double G1, double G2);  //Возвращает стационарную амплитуду во второй НЧ секции, амплитуда в первой присваивается A_stat
	cplx CurrentA2(double A1, double  _reA, double _imA);		//В двойной схеме, ток во второй секции?
	
	cplx FindBstatDoubleScheme(cplx b0, int Nb_it, double Gb, double Astat, cplx A2_stat);
	cplx CurrentB2(double  _reB, double _imB, double _A, cplx _A2stat); //ВЧ ток в двойной схеме

	double DeltaEnergyDoubleScheme();
public:
	MultiplierTreeCavity(QDomDocument* doc);
	double getHFoutputPowerDoubleScheme(double _Current_ampers, double _period, int _Nperiods, double _Ld1, double _Ld2, double _La2, double _Lb, double _k1, double _Norma, double Norma1B, double _voltage, double _inputPower_watts, double _delta, double _Qa1, double _Qa2, double _Qb, double _wall, char *filename, char *comment);
	double retriveBPower(cplx b, double omega, double Qb, double NormaB);
};
#endif