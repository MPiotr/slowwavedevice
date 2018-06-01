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
#define NR 16 //���������� ����� � ��������� ������������� ������� ����������������� ������ �� �������. �� ������ ���� ������ 512/4 = 128, ����� �� ����� ��������� ���������� rAk, iAk ...
#define NMESHMAX  128  //������������ ���������� ����� � ������� ��������� � twt_1d
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
struct ParamsM  //structure for packing all the parameters to be sent to device
{
	double la, lb, ld, k1, la1, la2, ld1, ld2, delta, h, voltage, L, wall, g1, g3, lossKappa, h_cycl, omega_cycl = 0, v_trans_max = 0, beamThickX, beamThickY, *enMax, *enMin;

	double *rJ3, *iJ3, *avEN, *int_rJ3, *int_iJ3, *int_rJ3_1, *int_iJ3_1, *tAr, *tAi, *Wmax, *Wmin;

	void *Amps;

	double *Qk, *Wk, *Q0, *W0, *ar0, *ai0, *ar0_t, *ai0_t, *rAk, *iAk, G, *int_rQ1, *int_iQ1, *rAq1, *iAq1, *rAq1k, *iAq1k; //���������������� �����

	double *radii;			//���������� ���������������� �����

	int *ifnotdestroyed;

	int Nz, Ngrid, Na1, Na2, Na3; // Na - additional dimensions for spreading parameters
};
#endif

class Device
{
protected:
	const static int NV = 4;

	int N, Nt, Nq, Np, Ns, Nx, Nv;
	unsigned int N_Wz, N_Qz;
	int Nmax;
	int fileIndex;
	double dTime, dz, Delta;
	double A_stat;
	double k0;
	double Lmax;
    double  grpSpeedCoeff;
	double g1;
	double g3;								//��������� ������� ���������
//
	double *d_rJ3, *d_iJ3, *d_avEN, *d_int_rJ3, *d_int_iJ3, *d_int_rJ3_1, *d_int_iJ3_1;
	double *d_Qz, *d_Wz, *d_tAr, *d_tAi, *tAr, *tAi;
	double *d_ar0, *d_ai0, *d_ar0_t, *d_ai0_t, *d_rAk, *d_iAk, *d_radii;
	double *d_Qk, *d_Wk, *d_Q0, *d_W0;
	double *d_Wmax, *d_Wmin, *d_int_rQ1, *d_int_iQ1, *d_rAq1, *d_iAq1, *d_rAq1k, *d_iAq1k;
	cplx   *d_Amps;
	int    *d_ifnotdestroyed;
	ParamsM *d_par;
	
	vector< vec > iteratedParams;               //���������, ������� ���� ��������� (������ ���������� (�������) - ����� ����� (������ ������) ��� ����������� ���������)
	vector< string > iteratedNames;             //����� ����������, ������� ���� ���������
	virtual void iterate(int) { throw ("iterate(int) is not implemented"); };                        //������� ����������, ������������ �� �� ���� �������
	virtual void changeParam(string, double);   //�������� �������� �� �����, ������������ ���, ��� ������������ iterate();
	double paramG();
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
	bool initialized;						//�������� ���� �� ��������� ������������� ����������.
	bool debug;							    //������ ���������� ����������
	bool firstRun = true;					//��� ������� ������ �����

	ParamsM setPar();
	virtual void PrintResults() { ; };
	FILE *results;

	double Current_ampers;					//���
	double period,  Nperiods;		    	//������, ����� �������� ������ ���    
	double structureWidth, beamWidth, beamHeight;       //������ ��������� � ������ ����� (��� ��� ��������)
	double group_speed, synch_angle;		//��������� �������� (��� ���/��� ��������) � ���� �����������
	double k1;								//�������� �����
	double Norma;							//����� ������ ������ (��� ������)
	double inputPower_watts, delta_freq;	//������� ��������, � ����� �������
	double voltage;							//����������
	double Q;							    //��������� ����������� ������ ������ (��) 
	double wall;							//���������� �� ������ � ������� �����
	double spchQ1;							//����������� ����� ���������������� ������� (������ ���������)
	char *solverName;						//����� �������
	char *fieldFileName;					//���� � ���������� ����������.
	char *dispersionFileName;				//���� � ���������� (����������� ������� ��� �� ��������� ����� � ��������)

	double s_stat;
	double monitor_par;

	void setParamsFromXML(char *xmlfilename);

public:

	Device ();
	Device (QDomDocument* doc);
	
	double A0, B0;

	bool initSolver(int nz, double lsolver, double groupSpeedCoeff, char *_solvername);
	
	void releaseDeviceMemory();

	void setCurrent(double current);


};

#endif