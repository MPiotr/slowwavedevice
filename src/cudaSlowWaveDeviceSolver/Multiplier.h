#ifndef KLINOTRON_H
#define KLINOTRON_H

#include <cuda.h>
#include <complex>
#include <driver_types.h>
#include <cuda_runtime_api.h>
#include<QtXml/QDomDocument>
#include<QtCore/QFile>
//#include "Objects.h"	
#include "synchronism.h"

//#define NP 64

#define gpuErrChk(ans){ gpuAssert( (ans), __FILE__, __LINE__);}
inline void gpuAssert(cudaError_t code, char *file, int line, bool abort = false)
{
	if (code != CUDA_SUCCESS)
	{
		printf("GPU Assert: %s, %s, line %i\n", cudaGetErrorString(code), file, line);
		if (abort) exit(code);
	}
}

#define NP 64//64debug
#define NQ 4
#define NS 2 
#define NR 16 //Количество точек в табличном представлении функции пространственного заряда от радиуса. не должна быть больше 512/4 = 128, иначе не будет правильно обнуляться rAk, iAk ...
#define NMESHMAX  128  //Максимальное количество точек в таблице амплитуде в twt_1d
#define SHARRAYSIZE 512

//#define NP 32
//#define NQ 8
//#define NS 2	
using namespace std;
#ifndef vec
#define vec vector<double>
#endif

#ifndef cplx
#define cplx complex<double>
#endif

#ifndef I
const cplx I = cplx(0., 1.);
#endif

#ifndef Pi
#define Pi 3.1415926535897932384626433832795
#endif

#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif

#ifndef c
const double c = 2.998e10;
#endif

#ifndef m 
const double m = 1e-27;
#endif

#ifndef e
const double e = 4.8e-10;
#endif

#ifndef ParamsM
struct ParamsM
{
	double la, lb, ld, k1, la1, la2, ld1, ld2, delta, h, voltage, L, wall, g1, g3, lossKappa, h_cycl, omega_cycl = 0, v_trans_max = 0, beamThickX, beamThickY, *enMax, *enMin;

	double *rJ3, *iJ3, *avEN, *int_rJ3, *int_iJ3, *int_rJ3_1, *int_iJ3_1, *tAr, *tAi, *Wmax, *Wmin;

	void *Amps;

	double *Qk, *Wk, *Q0, *W0, *ar0, *ai0, *ar0_t, *ai0_t, *rAk, *iAk, G, *int_rQ1, *int_iQ1, *rAq1, *iAq1, *rAq1k, *iAq1k; //пространственный заряд

	double *radii;			//двухмерный пространственный заряд

	int *ifnotdestroyed;

	int Nz, Ngrid, Na1, Na2, Na3; // Na - additional dimensions for spreading parameters
};
#endif

class Multiplier
{
protected:
	const static int NV = 4;

	int N, Nt, Nq, Np, Ns, Nx, Nv;
	int Namm; // количество  продольных мод для мультимодового солвера
	int Nmax;
	int fileIndex;
	double dTime, dz, Delta;
	double A_stat;
	double k0;
	double Lmax;
	cplx  B_stat;

    double  grpSpeedCoeff;

	double *d_rJ3, *d_iJ3, *d_avEN, *d_int_rJ3, *d_int_iJ3, *d_int_rJ3_1, *d_int_iJ3_1;
	double *d_Qz, *d_Wz, *d_tAr, *d_tAi, *tAr, *tAi;
	double *d_ar0, *d_ai0, *d_ar0_t, *d_ai0_t, *d_rAk, *d_iAk, *d_radii;
	double *d_Qk, *d_Wk, *d_Q0, *d_W0;
	double *d_Wmax, *d_Wmin, *d_int_rQ1, *d_int_iQ1, *d_rAq1, *d_iAq1, *d_rAq1k, *d_iAq1k;
	cplx *d_Amps;
	int *d_ifnotdestroyed;
	ParamsM *d_par;
	unsigned int N_Wz, N_Qz;

	double avKPD;

	double g1, g3;

	vector< vec > iteratedParams;
	vector< string > iteratedNames;
	void changeParam(string, double);
	inline double paramG();
	double power(double kpd);
	double efficieny;
	
	cplx *A, *d_A, *b;
	cplx *CurrentData = NULL, *crnt = NULL;
	double *structReal = NULL, *structIm = NULL, *longitudinalStructureRe = NULL, *longitudinalStructureIm = NULL, *qStructure = NULL, *mesh = NULL;
	
	Interpolation *longStructRealRe = NULL;
	Interpolation *longStructRealIm = NULL;
	Interpolation *qStructRe = NULL;
	Interpolation *qStructIm = NULL;

	Synchronism* syncwave = NULL;

	bool solverStatus;
	bool fieldLoaded;
	bool printTablesFlag;
	bool debflag;
	bool initialized;						//проверка была ли проведена инициализация параметров.
	bool debug;							    //печать отладочной информации
	bool firstRun = true;					//для лучшего вывода логов

	ParamsM setPar();
	void PrintResults();


	FILE *results;

	double Current_ampers;					//Ток
	double period,  Nperiods;		    	//период, число периодов первой сек
    double Ld, Lb;							//длина дрейфовой секции, длина второй (ВЧ) секции
	double structureWidth, beamWidth, beamHeight;       //ширина структуры и ширина пучка (для ЛБВ решателя)
	double group_speed, synch_angle;		//групповая скорость (для ЛБВ/ЛОВ решателя) и угол синхронизма
	double La2, Ld1, Ld2;				    //Длины НЧ резонаторных и дрейфовых секций в трёхклистронной схеме
	double k1;								//волновое число
	double Norma;							//норма первой секции (вся секция)
	double Norma1, Norma2;					//норма нч секций (целиком) в трёхрезонаторной схеме
	double Norma1b;							//норма вч секции (один период)
	double inputPower_watts, delta_freq;	//входная мощность, и сдвиг частоты
	double voltage;							//Напряжение
	double Q, Qb;							//Добротность первой (НЧ) и омическая добротность второй (ВЧ) секций (двухрезонаторная схема)
	double Qa1, Qa2;						//Добротность первой НЧ и второй НЧ секций (трёхрезонаторная схема)
	double wall;							//расстояние от стенки и толщина пучка
	double spchQ1;							//Коэффициент перед пространственным зарядом (первая гармоника)
	double spchQ3;							//Коэффициент перед пространственным зарядом (третья гармоника)
	char *solverName;						//Выбор солвера
	char *fieldFileName;					//файл с поперечной структурой.
	char *dispersionFileName;				//файл с дисперсией (зависимость частоты ГГц от волнового числа в градусах)
	char difractionFlag[200];
	double Q_difr;


	double s_stat;
	double monitor_par;

	double diffractionQ(char *);				//Задаёт значение дифракционной добротности (проверяет соответствующий флаг)

	virtual double findAstat(double nextA, int Nit, double G, double *resS = 0);
	double findAstatDetuned(double nextA, int Nit, double G, double S = 0);
	double findAstatDetuned_full(double nextA, int Nit, double G); //в отличие от предыдущей функции меняет амплитуду входящего тока

	cplx FindBstat(cplx b0,			int Nb_it,  double Gb, double Astat);
	cplx FindBstatDetuned(cplx b0,	int Nb_it,  double Gb, double Astat);

	double ElectronsDeltaEnergy(double _A); //warper для изменения средней энергии электронов
	double DeltaEnergy(double _A);
	

	cplx retriveBCurr();
	void retriveBCurr(cplx* J1, cplx* J2);						//Для решателя с двумя продольными гармониками, 
	double retriveDeltaEnergy();

	virtual cplx ElectronCurrentA(double reA, double imA);          //warper для ВЧ тока на первой гармонике
	cplx CurrentA(double  _reA, double _imA);				//НЧ ток (интеграл)
	
	cplx CurrentAMultiModes(double *X, int Na = 3, cplx *J1 = 0, cplx * J2 = 0);				//НЧ ток (интеграл) 
	void CurrentAMultiModes(std::complex<double> *Amps, std::complex<double> * currs, double *buffRe, double *buffIm, int Na = 3, cplx *J1 = 0, cplx *J2 = 0); // Для многомодового решателя перегруженный
	void retriveACurrComplex(std::complex<double>  *Amps, std::complex<double>  *currs, double *currsBuffRe, double *currsBuffIm, int Na, int Nstop);

	virtual cplx HfElectronCurrent(double  _reB, double _imB, double _A); //warper для ВЧ тока
	cplx CurrentB(double  _reB, double _imB, double _A);    //ВЧ ток (интеграл)

	cplx findAstatDoubleScheme(double nextA, int Nit, double G1, double G2);  //Возвращает стационарную амплитуду во второй НЧ секции, амплитуда в первой присваивается A_stat

	void findAstatMultiModes(double* nextX, int N_it, double G1, double G2); //переписывает Х

	cplx FindBstatDoubleScheme(cplx b0,	int Nb_it,  double Gb, double Astat, cplx A2_stat);

	double DeltaEnergyDoubleScheme();

	cplx CurrentA2(double A1, double  _reA, double _imA);		//В двойной схеме, ток во второй секции?
	cplx CurrentB2(double  _reB, double _imB, double _A, cplx _A2stat); //ВЧ ток в двойной схеме

	void PrintParams(char *filename, char *comment);
	void PrintParamsDoubleScheme(char *filename, char*comment);
	void PrintParamsMultiModes(char *filename, char *comment);

	void setParamsFromXML(char *xmlfilename);

public:

	Multiplier ();
	Multiplier (QDomDocument* doc);
	
//	~Multiplier();

	double A0, B0;

	double getHFoutputPower(double _Current_ampers, double _period, int _Nperiods, double _Ld, double _Lb, double _k1, double _Norma, double _NormaB, double _voltage, double _inputPower_watts, double _delta, double _Qa, double _Qb, double _wall, double _sQ1, double _sQ3, char *filename, char *comment, double *monitor = 0, char *difrFlag = 0);

	double getHFoutputPowerDoubleScheme(double _Current_ampers, double _period, int _Nperiods, double _Ld1, double _Ld2, double _La2, double _Lb, double _k1, double _Norma, double Norma1B, double _voltage, double _inputPower_watts, double _delta, double _Qa1, double _Qa2, double _Qb, double _wall, char *filename, char *comment);

	double getHFoutputPowerMultiModes(double _Current_ampers, double _period, int _Nperiods, double _Ld, double _Lb, double _k1, double _Norma, double Norma2, double _voltage, double _inputPower_watts, double _delta, double _Qa1, double _Qa2, double _Qb, double _wall, char *filename, char *comment);

	double retriveBPower(cplx b, double omega, double Qb, double NormaB);

	bool initMultiplierSolver(int nz, double lsolver, double groupSpeedCoeff, char *_solvername);
	
 	void releaseMultiplierMemory();

	void setCurrent(double current);


};

#endif