#include "twt_0d.h"


TWT_0D::TWT_0D(QDomDocument *doc) :TWT_1D(doc)
{

}

TWT_0D::TWT_0D(QDomDocument *doc, TWT_0D *instance) : TWT_1D(doc, instance)
{
	;
}

TWT_0D::TWT_0D(QDomDocument *doc, int a) : TWT_1D(doc, a) //без инициализации CUDA солвера
{
	;
}

cplx TWT_0D::solve(cplx *A, double *ar, double *ai, double inputAmp, double lossKappa, double delta,
	double *fieldStructureRe, double *fieldStructureIm, double G,
	double enPrint, bool printField,
	double *longStrRe, double *longStrIm, double *qStr, double*mesh)
{
	return solveTWT_0d(A, ar, ai, inputAmp, lossKappa, delta, fieldStructureRe, fieldStructureIm, mesh, G, enPrint, printField, longStrRe, longStrIm, qStr);
}

void TWT_0D::readTransversalStructure()
{
	double* strAr = new double;
	double* strAi = new double;

	*strAi = 0;
	*strAr = 1./sqrt(2 * Norma / double(Nperiods));      //*2/Nperiods : -- Norma is the norma of the full section (input 'periodNorma' is multiplierd by Nperiods and divided by 2, this is backward transformation)		
	
	structReal = strAr;
	structIm = strAi;

	return;
}