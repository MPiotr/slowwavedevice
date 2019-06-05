#include "twt.h"
#include <QtCore\qsharedmemory.h>

#ifndef TWT_2D__H
#define TWT_2D__H
using namespace std;
extern void parseFieldFile(FILE *fieldFile, int numX, int NbeamY, double beamWidth, double beamHeight, double width, double **structReal, double **structImag);
extern void parseTable(FILE *fieldFile, int Nmax, double dz, double **structReal, double **structImag);
extern void parseTable(FILE *fieldFile, Interpolation **structReal, Interpolation **structImag);

class BWO_2D;

class TWT_2D : public TWT
{
	bool initSolver(int, double);
	cplx solveTWT_2d(cplx *A, double *ar, double *ai, double inputAmp, double lossKappa, double delta,
		double *fieldStructureRe, double *fieldStructureIm, double G,
		double enPrint, bool printField = false,
		double *longStrRe = NULL, double *longStrIm = NULL, double *qStr = NULL); //������ ��������� ��� ���
protected:
//	vector<BWO_2D> parasites;

//	double *ar, *ai;

	virtual cplx solve(cplx *A, double *ar, double *ai, double inputAmp, double lossKappa, double delta,
					double *fieldStructureRe, double *fieldStructureIm, double G,
					double enPrint, bool printField = false,
					double *longStrRe = NULL, double *longStrIm = NULL, double *qStr = NULL, double *mesh = NULL);
//	void iterate(int);
	void readTransversalStructure();	
public:
	TWT_2D(QDomDocument *doc);
	TWT_2D(QDomDocument *doc, TWT_2D *instance);
	TWT_2D(QDomDocument *doc, int a); //������ �������� ��� ������������� CUDA �������
};
#endif