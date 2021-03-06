#ifndef __TWT_0D_H
#define __TWT_0D_H

#include "twt_1d.h"
class TWT_0D : public TWT_1D
{
	bool initSolver(int nz, double lsolver);
protected:
	void changeParam(string, double);
	void printCurrentParams(FILE * file);
	void printParamsHeader(FILE *file);

	double clinotronAngle; 
	Interpolation *clinotronShiftStrRe = NULL;
	Interpolation *clinotronShiftStrIm = NULL;
	double *clinotronStructure = NULL;

	void generateClinotronStructure(double h);

	int NumMesh;
	cplx solveTWT_0d(cplx *A, double *ar, double *ai, double inputAmp, double lossKappa, double delta,
		double *fieldStructureRe, double *fieldStructureIm, double *mesh, double G,
		double enPrint, bool printField = false,
		double *longStrRe = NULL, double *longStrIm = NULL, double *qStr = NULL); //������ ��������� ��� ���
	cplx solve(cplx *A, double *ar, double *ai, double inputAmp, double lossKappa, double delta,
		double *fieldStructureRe, double *fieldStructureIm, double G,
		double enPrint, bool printField = false,
		double *longStrRe = NULL, double *longStrIm = NULL, double *qStr = NULL, double *mesh = NULL);
	void readTransversalStructure();	
public:
	TWT_0D(QDomDocument *doc);
	TWT_0D(QDomDocument *doc, TWT_0D *instance);
	TWT_0D(QDomDocument *doc, int a); //��� ������������� CUDA �������

};
#endif