#include "twt_0d.h"
#include <io.h>
#include "xml_routines.h"


TWT_0D::TWT_0D(QDomDocument *doc) :TWT_1D(doc)
{
	if (!setXMLEntry(doc, "clinotronAngle", &clinotronAngle, &iteratedParams, &iteratedNames)) clinotronAngle = 0;
	char clinotronStructureFile[200];
	if (setXMLEntry(doc, "clinotronStructureTable", (char*)clinotronStructureFile))
	{
		if (_access(clinotronStructureFile, 0) == 0)
		{
			FILE *strFile = fopen(clinotronStructureFile, "r");
			parseTable(strFile, &clinotronShiftStrRe, &clinotronShiftStrIm);
			fclose(strFile);
			clinotronStructure = new double[Nmax];
		}

	}
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
void TWT_0D::changeParam(string	 name, double value)
{
	Device::changeParam(name, value);
	if (name == "clinotronAngle") clinotronAngle = value;
}

void TWT_0D::printCurrentParams(FILE *file)
{
	TWT_1D::printCurrentParams(file);
	fprintf(file, "%g,", clinotronAngle);

}

void TWT_0D::printParamsHeader(FILE *file)
{
	TWT_1D::printParamsHeader(file);
	fprintf(file, "%g,", clinotronAngle);
}

void TWT_0D::generateClinotronStructure(double h)
{
	if (clinotronStructure == NULL) return;
	double La = Nperiods*period*h;
	double dz = Lmax / double(Nmax);
	int Nstop = ceil(La / dz);
	for (int i = 0; i < Nmax; i++)
	{
		double hz = i*dz;
		clinotronStructure[i] = clinotronShiftStrRe->at(hz / h);
	}
}
