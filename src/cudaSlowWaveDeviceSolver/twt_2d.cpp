#include <io.h>
#include "twt_2d.h"
#include "xml_routines.h"

void parseFieldFile(FILE *fieldFile, int numX, int NbeamY, double beamWidth, double beamHeight, double width, double **structReal, double **structImag)
{
	std::vector<double> struct1d_fieldRe;
	std::vector<double> struct1d_fieldIm;
	std::vector<double> struct1d_y;

	float y, fieldRe, fieldIm;
	while (fscanf(fieldFile, "%g,%g,%g\n", &y, &fieldRe, &fieldIm) == 3)
	{
		struct1d_fieldRe.push_back(fieldRe);
		struct1d_fieldIm.push_back(fieldIm);
		struct1d_y.push_back(y);
	}
	int numY = struct1d_y.size();

	Interpolation reField(struct1d_y.data(), struct1d_fieldRe.data(), numY);
	Interpolation imField(struct1d_y.data(), struct1d_fieldIm.data(), numY);


	double *strReal = new double[numX*NbeamY];
	double *strImag = new double[numX*NbeamY];


	int ny = 0;
	for (int y = 0; y < NbeamY; y++)
	{
		double c_y = double(y - NbeamY / 2) / double(NbeamY)*beamHeight;
		for (int x = 0; x < numX; x++)
		{
			double X = (double(x) / double(numX) - 0.5)*beamWidth;
			strReal[y*numX + x] = reField.at(c_y)*cos(Pi / width*X);
			strImag[y*numX + x] = imField.at(c_y)*cos(Pi / width*X);
		}
	}
	*structReal = strReal;
	*structImag = strImag;
}
void parseTable(FILE *fieldFile, int Nmax, double dz, double **structReal, double **structImag)
{
	std::vector<double> struct1d_fieldRe;
	std::vector<double> struct1d_fieldIm;
	std::vector<double> struct1d_y;

	float y, fieldRe, fieldIm;
	while (fscanf(fieldFile, "%g,%g,%g\n", &y, &fieldRe, &fieldIm) == 3)
	{
		struct1d_fieldRe.push_back(fieldRe);
		struct1d_fieldIm.push_back(fieldIm);
		struct1d_y.push_back(y);
	}
	int size = struct1d_y.size();

	Interpolation reField(struct1d_y.data(), struct1d_fieldRe.data(), size);
	Interpolation imField(struct1d_y.data(), struct1d_fieldIm.data(), size);


	double *strReal = new double[Nmax];
	double *strImag = new double[Nmax];
	memset(strReal, 0, sizeof(double)*Nmax);
	memset(strImag, 0, sizeof(double)*Nmax);


	for (int i = 0; i < Nmax; i++)
	{
		double z = i*dz;
		strReal[i] = reField.at(z);
		strImag[i] = imField.at(z);
		
	}
	*structReal = strReal;
	*structImag = strImag;
}

void unitityStruct(int numX, int NbeamY, double norma, double **structReal, double **structImag)
{
	double *strReal = new double[numX*NbeamY];
	double *strImag = new double[numX*NbeamY];

	for (int y = 0; y < NbeamY; y++)
	{
		for (int x = 0; x < numX; x++)
		{
			strReal[y*numX + x] = 1./norma;
			strImag[y*numX + x] = 0;
		}
	}
	*structReal = strReal;
	*structImag = strImag;

}
TWT_2D::TWT_2D(QDomDocument *doc) : TWT(doc)
{
	initSolver(Nmax, Lmax);

// ...............reading parasites data ........................
/*	int NumParasites = doc->elementsByTagName("parasite").length();
	for (int i = 0; i < NumParasites; i++)
	{
		QDomNode parasite = doc->elementsByTagName("parasite").item(i);
		int chLength = parasite.childNodes().length();
		
		char dispFile[200];
		char fieldFile[200];
		if( !setXMLEntry(&parasite, "parasiteDispersion", (char*) dispFile))  continue;
		if (!setXMLEntry(&parasite, "parasiteField", fieldFile)) continue;

		syncwave->addMode(dispFile);
		
		parasites.push_back(BWO_2D(doc,this, i+1, fieldFile));
	}*/
// ........................................................

// .........  reading transversal structure .............
	readTransversalStructure();
	
}
TWT_2D::TWT_2D(QDomDocument *doc, int a) : TWT(doc)
{
	
}
TWT_2D::TWT_2D(QDomDocument *doc,TWT_2D *copy) :TWT_2D(doc, 0)
{

	if (inputPower_watts == 0) {
		printf("Warning: input power is zero for TWT solver\n");
		initialized = false;
	}
	char filename[200];
	sprintf(filename, "twt_result_%i.dat", fileIndex);
	results = fopen(filename, "w");
	printDataHeader(results);

	d_par = copy->d_par;

	d_rJ3 = copy->d_rJ3; d_iJ3 = copy->d_iJ3, d_avEN = copy->d_avEN; 
	d_int_rJ3 = copy->d_int_rJ3; d_int_iJ3 = copy->d_int_iJ3;
	d_int_rJ3_1 = copy->d_int_rJ3_1; d_int_iJ3_1 = copy->d_int_iJ3_1;
	
	d_Qz = copy->d_Qz; d_Wz = copy->d_Wz; d_tAr = copy->d_tAr; d_tAi = copy->d_tAi, tAr = copy->tAr; tAi = copy->tAi;

	d_ar0 = copy->d_ar0; d_ai0 = copy->d_ai0; d_ar0_t = copy->d_ar0_t; d_ai0_t = copy->d_ai0_t;
	d_rAk = copy->d_rAk; d_iAk = copy->d_iAk; d_radii = copy->d_radii; 
	d_fAr = copy->d_fAr; d_fAi = copy->d_fAi;

	d_Qk = copy->d_Qk; d_Wk = copy->d_Wk;   d_Q0 = copy->d_Q0;  d_W0 = copy->d_W0;

	d_Wmax = copy->d_Wmax; d_Wmin = copy->d_Wmin; d_int_rQ1 = copy->d_int_rQ1; d_int_iQ1 = copy->d_int_iQ1; d_rAq1 = copy->d_rAq1; d_iAq1 = copy->d_iAq1; d_rAq1k = copy->d_rAq1k; d_iAq1k = copy->d_iAq1k;
	
	d_Amps = copy->d_Amps;
	d_ifnotdestroyed = copy->d_ifnotdestroyed;


}
void TWT_2D::readTransversalStructure()
{

	if (_access(fieldFileName, 0) == 0)
	{
		FILE *fieldFile = fopen(fieldFileName, "r");
		if (!fieldFile) {
			printf("\nError opening field file: %s \n", fieldFileName);
			initialized = false;
		}
		parseFieldFile(fieldFile, Nq, Ns, beamWidth, beamHeight, structureWidth, &structReal, &structIm);
		fclose(fieldFile);
	}
	else
	{
		unitityStruct(Nq, Ns, sqrt(2 * Norma / double(Nperiods)), &structReal, &structIm);    //*2/Nperiods : -- Norma is the norma of the full section (multiplier legacy: input 'periodNorma' is multiplierd by Nperiods and divided by 2, this is backward transformation)
	}

}

cplx TWT_2D::solve(cplx *A, double *ar, double *ai, double inputAmp, double lossKappa, double delta,
	double *fieldStructureRe, double *fieldStructureIm, double G,
	double enPrint, bool printField,
	double *longStrRe, double *longStrIm, double *qStr, double *mesh)
{
	return solveTWT_2d(A, ar, ai, inputAmp, lossKappa, delta, fieldStructureRe, fieldStructureIm, G, enPrint, printField, longStrRe, longStrIm, qStr);
}

