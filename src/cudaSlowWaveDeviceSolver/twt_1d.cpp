#include <io.h>
#include "twt_1d.h"
#include "xml_routines.h"

TWT_1D::TWT_1D(QDomDocument *doc) :TWT_2D(doc)
{
	initSolver(Nmax, Lmax);
	readTransversalStructure();

	if (setXMLEntry(doc, "magneticField", &bTesla))	{
		omega_cycl = e*bTesla*1.e4 / (m*c);
		double gamma = 1 + voltage / 511.;
		double beta = sqrt(gamma*gamma - 1) / gamma;
		h_cycl = omega_cycl/(beta*c);

	}
	else {
		bTesla = 0, omega_cycl = 0, h_cycl = 0;
	}
	if (setXMLEntry(doc, "maxTransversalBeta", &v_trans_max))	{
		printf("maximal cyclortorn radius is %g mm\n", v_trans_max*c*10. / (omega_cycl));
	}
	else {
		v_trans_max = 0;
	}

}
TWT_1D::TWT_1D(QDomDocument *doc, TWT_1D *instance) : TWT_2D(doc, instance)
{
	;
}

TWT_1D::TWT_1D(QDomDocument *doc, int a):TWT_2D(doc, a) //без инициализации CUDA солвера
{
	;
}

void parseFieldFile1d(FILE *fieldFile, int *NumY, double **structReal, double **structImag, double **mesh)
{
	std::vector<double> struct1d_fieldRe;
	std::vector<double> struct1d_fieldIm;
	std::vector<double> struct1d_y;

	float y, fieldRe, fieldIm;
	int index = 0;
	while (fscanf(fieldFile, "%g,%g,%g\n", &y, &fieldRe, &fieldIm) == 3)
	{
		struct1d_fieldRe.push_back(fieldRe);
		struct1d_fieldIm.push_back(fieldIm);
		struct1d_y.push_back(y);
		index++;
	}
	int numY = struct1d_y.size();

	double *strReal = new double[numY];
	double *strImag = new double[numY];
	double *strMesh = new double[numY];
	
	int ny = 0;
	for (int y = 0; y < numY; y++)
	{
		strMesh[y] = struct1d_y[y];
		strReal[y] = struct1d_fieldRe[y];
		strImag[y] = struct1d_fieldIm[y];
	}
	*structReal = strReal;
	*structImag = strImag;
	*mesh = strMesh;
	*NumY = numY;
}

void generateHomogeneousStructure(int Ns, int *NumY, double constant, double **structReal, double **structImag, double **mesh)
{
	double *strReal = new double[Ns];
	double *strImag = new double[Ns];
	double *strMesh = new double[Ns];
	int numY = Ns;
	*NumY = numY;
	for (int i = 0; i < numY; i++)
	{
		strMesh[i] = i;
		strReal[i] = 1./constant;
		strImag[i] = 0;
	}
	*structReal = strReal;
	*structImag = strImag;
	*mesh = strMesh;
}
void TWT_1D::readTransversalStructure()
{
	A = new cplx[Nmax];
	ar = new double[Nmax];
	ai = new double[Nmax];
	memset(ar, 0, sizeof(double)*Nmax);
	memset(ai, 0, sizeof(double)*Nmax);
	memset(A, 0, sizeof(cplx)*Nmax);

	if (_access(fieldFileName, 0) == 0)
	{
		FILE *fieldFile = fopen(fieldFileName, "r");
		if (!fieldFile) {
			printf("\nError opening field file: %s \n", fieldFileName);
			initialized = false;
		}
		parseFieldFile1d(fieldFile, &NumMesh, &structReal, &structIm, &mesh);
		fclose(fieldFile);
	}
	else
	{
		generateHomogeneousStructure(Ns, &NumMesh, sqrt(2 * Norma / double(Nperiods)), &structReal, &structIm, &mesh);      //*2/Nperiods : -- Norma is the norma of the full section (input 'periodNorma' is multiplierd by Nperiods and divided by 2, this is backward transformation)		
	}

}
cplx TWT_1D::solve(cplx *A, double *ar, double *ai, double inputAmp, double lossKappa, double delta,
	double *fieldStructureRe, double *fieldStructureIm, double G,
	double enPrint, bool printField,
	double *longStrRe, double *longStrIm, double *qStr, double*mesh)
{
	return solveTWT_1d(A, ar, ai, inputAmp, lossKappa, delta, fieldStructureRe, fieldStructureIm, mesh, G, enPrint, printField, longStrRe, longStrIm, qStr);
}
