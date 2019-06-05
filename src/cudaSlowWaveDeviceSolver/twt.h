#include "device.h"
#include <QtCore\qsharedmemory.h>

#ifndef __TWT_H
#define __TWT_H
class TWT :public Device
{
	bool initSolver(int nz, double lsolver);
protected:
//	vector<BWO_2D> parasites;

	cplx *A;
	double *ar, *ai;

	double  *d_fAr, *d_fAi, *d_mesh;

	virtual cplx solve(cplx *A, double *ar, double *ai, double inputAmp, double lossKappa, double delta,
		double *fieldStructureRe, double *fieldStructureIm, double G,
		double enPrint, bool printField = false,
		double *longStrRe = NULL, double *longStrIm = NULL, double *qStr = NULL, double *mesh = NULL) = 0; //решает уравнения для ЛБВ
	virtual void iterate(int);

	virtual void readTransversalStructure() = 0;
	void generateLongitudinalStructure(double h);
	void generateQStructure(double h);

	virtual void printResults(FILE *file, cplx *A);
	virtual void printCurrentParams(FILE *file);
	virtual void printParamsHeader(FILE *file);
	virtual void printResultHeader(FILE *file);
	virtual void printDataHeader(FILE *file);

	double paramG(double h);
	double power(double absAmplitude, double h);
	double efficiency(double absAmplitude, double h);
	double amplitude(double inpPower, double h);
	double lossKappa(double h);

	QSharedMemory *sharedMemory;
	bool shMemoryCreated = false;
	void printAbsAtoSharedMemory(int N);
public:
	TWT(QDomDocument *doc);
	TWT(QDomDocument *doc, TWT *instance);
	TWT(QDomDocument *doc, int a); //без инициализации CUDA солвера
	double solveTWT();
};

#endif