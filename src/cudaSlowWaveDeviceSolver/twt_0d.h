#ifndef __TWT_0D_H
#define __TWT_0D_H

#include "twt_1d.h"
class TWT_0D : public TWT_1D
{
protected:
	double bTesla, omega_cycl, h_cycl, v_trans_max;

	int NumMesh;
	cplx solveTWT_0d(cplx *A, double *ar, double *ai, double inputAmp, double lossKappa, double delta,
		double *fieldStructureRe, double *fieldStructureIm, double *mesh, double G,
		double enPrint, bool printField = false,
		double *longStrRe = NULL, double *longStrIm = NULL, double *qStr = NULL); //решает уравнения для ЛБВ
	cplx solve(cplx *A, double *ar, double *ai, double inputAmp, double lossKappa, double delta,
		double *fieldStructureRe, double *fieldStructureIm, double G,
		double enPrint, bool printField = false,
		double *longStrRe = NULL, double *longStrIm = NULL, double *qStr = NULL, double *mesh = NULL);
	void readTransversalStructure();
	bool initSolver(int nz, double lsolver);
public:
	TWT_0D(QDomDocument *doc);
	TWT_0D(QDomDocument *doc, TWT_0D *instance);
	TWT_0D(QDomDocument *doc, int a); //без инициализации CUDA солвера

};
#endif