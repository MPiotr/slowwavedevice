#ifndef __MULTIPLIER_H
#define __MULTIPLIER_H
#include "device.h"

class Multiplier : public Device
{
protected:
	cplx  B_stat;
	double Ld, Lb;							//длина дрейфовой секции, длина второй (ВЧ) секции
	double Norma1b;							//норма вч секции (один период)
	double Qb;								//Омическая добротность ВЧ секции
	double Q_difr;							//Дифракционная добротность ВЧ секции
	double spchQ3;							//Коэффициент перед пространственным зарядом (третья гармоника)	
	double avKPD;							//КПД (усреднённый)
	char difractionFlag[200];

	double diffractionQ(char *);				//Задаёт значение дифракционной добротности (проверяет соответствующий флаг)

	double ElectronsDeltaEnergy(double _A); // обёртка для изменения средней энергии электронов
	double DeltaEnergy(double _A);

	virtual cplx ElectronCurrentA(double reA, double imA);          //обёртка для ВЧ тока на первой гармонике
	cplx CurrentA(double  _reA, double _imA);				        //НЧ ток (интеграл)

	virtual cplx HfElectronCurrent(double  _reB, double _imB, double _A); //обёртка для ВЧ тока
	cplx CurrentB(double  _reB, double _imB, double _A);                  //ВЧ ток (интеграл)

	cplx FindBstat(cplx b0, int Nb_it, double Gb, double Astat);
	cplx FindBstatDetuned(cplx b0, int Nb_it, double Gb, double Astat);

	virtual void PrintParams(char *filename, char *comment);

	virtual double findAstat(double nextA, int Nit, double G, double *resS = 0);	
	virtual double findAstatDetuned(double nextA, int Nit, double G, double S = 0);
	virtual double findAstatDetuned_full(double nextA, int Nit, double G); //в отличие от предыдущей функции меняет амплитуду входящего тока

	cplx retriveBCurr();
	void retriveBCurr(cplx* J1, cplx* J2);						//Для решателя с двумя продольными гармониками, 
	double retriveDeltaEnergy();
public:
	Multiplier(QDomDocument* doc);

	double getHFoutputPower(double _Current_ampers, double _period, int _Nperiods, double _Ld, double _Lb, double _k1, double _Norma, double _NormaB, double _voltage, double _inputPower_watts, double _delta, double _Qa, double _Qb, double _wall, double _sQ1, double _sQ3, char *filename, char *comment, double *monitor = 0, char *difrFlag = 0);
	
	double retriveBPower(cplx b, double omega, double Qb, double NormaB);

};

#endif