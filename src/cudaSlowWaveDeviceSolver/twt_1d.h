#ifndef __TWT_2D_H
#define __TWT_2D_H
#include "twt_2d.h"
class TWT_1D : public TWT_2D
{
protected: 
	double bTesla, omega_cycl, h_cycl, v_trans_max;

	int NumMesh;
	cplx solveTWT_1d(cplx *A, double *ar, double *ai, double inputAmp, double lossKappa, double delta,
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
	TWT_1D(QDomDocument *doc);
	TWT_1D(QDomDocument *doc, TWT_1D *instance);
	TWT_1D(QDomDocument *doc, int a); //без инициализации CUDA солвера
	
};
#endif