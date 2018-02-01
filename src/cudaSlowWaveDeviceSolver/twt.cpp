#include "twt.h"
#include <io.h>
#include "xml_routines.h"

void parseTable(FILE *fieldFile, Interpolation **structReal, Interpolation **structImag)
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

	*structReal = new Interpolation(struct1d_y.data(), struct1d_fieldRe.data(), size);
	*structImag = new Interpolation(struct1d_y.data(), struct1d_fieldIm.data(), size);
}

TWT::TWT(QDomDocument *doc) : Multiplier(doc)
{
	initSolver(Nmax, Lmax);

	if (inputPower_watts == 0) {
		printf("Warning: input power is zero for TWT solver\n");
		initialized = false;
	}
	// ........ shared memory for result exchange
	sharedMemory = new QSharedMemory("slowwavedeviceresults");
	int size = sizeof(int) + Nmax * 2 * sizeof(cplx);    //int(flag)<<int(Nstop)<<Nstop*cplx(A)
	if (sharedMemory->isAttached())	sharedMemory->detach();
	try {
		shMemoryCreated = sharedMemory->create(size);
	}
	catch ( std::exception  exeption) {
		throw("Error creating shared memory segment for calculation results\n");
	}
	if (!shMemoryCreated)
		throw("Error creating shared memory segment for calculation results\n");
	sharedMemory->lock();
	memset(sharedMemory->data(), 0, 2 * sizeof(int));
	sharedMemory->unlock();

	A = new cplx[Nmax];
	ar = new double[Nmax];
	ai = new double[Nmax];
	memset(ar, 0, sizeof(double)*Nmax);
	memset(ai, 0, sizeof(double)*Nmax);
	memset(A, 0, sizeof(cplx)*Nmax);

	// ...............reading parasites data ........................
/*	int NumParasites = doc->elementsByTagName("parasite").length();
	for (int i = 0; i < NumParasites; i++)
	{
		QDomNode parasite = doc->elementsByTagName("parasite").item(i);
		int chLength = parasite.childNodes().length();

		char dispFile[200];
		char fieldFile[200];
		if (!setXMLEntry(&parasite, "parasiteDispersion", (char*)dispFile))  continue;
		if (!setXMLEntry(&parasite, "parasiteField", fieldFile)) continue;

		syncwave->addMode(dispFile);

		parasites.push_back(BWO_2D(doc, this, i + 1, fieldFile));
	}*/
	// ........................................................

	// .........  reading transversal structure .............
	
	// ........................................................

	// ...................reading longitudinal structure.......
	char logitudinalStructureFile[200];
	if (setXMLEntry(doc, "logitudinalStructureFile", (char*)logitudinalStructureFile))
	{
		if (_access(logitudinalStructureFile, 0) == 0)
		{
			FILE *longStrFile = fopen(logitudinalStructureFile, "r");
			parseTable(longStrFile, &longStructRealRe, &longStructRealIm);
			fclose(longStrFile);
			longitudinalStructureRe = new double[Nmax];
			longitudinalStructureIm = new double[Nmax];
		}

	}
	// ........................................................
	// ...................reading longitudinal structure of Q factor.......
	char QStructureFile[200];
	if (setXMLEntry(doc, "QStructureFile", (char*)QStructureFile))
	{
		if (_access(QStructureFile, 0) == 0)
		{
			FILE *qStrFile = fopen(QStructureFile, "r");
			parseTable(qStrFile, &qStructRe, &qStructIm);
			fclose(qStrFile);
			qStructure = new double[Nmax];
		}

	}
	// ........................................................


}
TWT::TWT(QDomDocument *doc, int a) : Multiplier(doc)
{

	if (inputPower_watts == 0) {
		printf("Warning: input power is zero for TWT solver\n");
		initialized = false;
	}
	char filename[200];
	sprintf(filename, "twt_result_%i.dat", fileIndex);
	results = fopen(filename, "w");
	printDataHeader(results);


	A = new cplx[Nmax];
	ar = new double[Nmax];
	ai = new double[Nmax];
	memset(ar, 0, sizeof(double)*Nmax);
	memset(ai, 0, sizeof(double)*Nmax);
	memset(A, 0, sizeof(cplx)*Nmax);

}
TWT::TWT(QDomDocument *doc, TWT *copy) :TWT(doc, 0)
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
double TWT::paramG(double h)
{
	return period / (pow(h, 2)*group_speed) *Current_ampers / 17045.; //FQield norm is in field structure!
}
double TWT::power(double absAmplitude, double h)
{
	double gamma = 1. + voltage / 511.;
	return 1000.*voltage*Current_ampers*efficiency(absAmplitude, h);//	0.5*pow(absAmplitude, 2)*1000.*voltage*Current_ampers / (paramG(h) * (gamma - 1.));
}
double TWT::efficiency(double absAmplitude, double h)
{
	double gamma = 1. + voltage / 511.;
	return 0.5*pow(absAmplitude, 2) / (paramG(h)*(gamma - 1.));
}
double TWT::amplitude(double inpPower, double h)
{
	double gamma = 1. + voltage / 511.;
	//	return 2.*sqrt(2.*inpPower*paramG(h) / 511000.);
	return sqrt(2.*paramG(h)*inpPower*(gamma - 1) / (Current_ampers*voltage * 1000));
}
double TWT::lossKappa(double h)
{
	return k1 / h / (2 * Q*group_speed);
}

void TWT::printResults(FILE *file, cplx *A)
{
	double Amax, Pmax, Lmx;
	double Aout, Pout;
	double eff;

	double dz = Lmax / double(Nmax);
	double h = 2 * Pi / period*(synch_angle / 360.);
	double G = paramG(h);
	double La = Nperiods*period*h;
	int Nstop = ceil(La / dz);
	if (Nstop > Nmax) Nstop = Nmax - 1;
	int nd = round(period / dz);

	cplx *maxElement = max_element(A, A + Nstop - 1, [](cplx a, cplx b) {return (abs(a) < abs(b)); });
	int maxIndex = maxElement - A;
	Amax = abs(*maxElement);
	Lmx = dz*maxIndex;
	Aout = abs(A[Nstop - 1]);
	Pmax = power(Amax, h);
	Pout = power(Aout, h);
	eff = efficiency(Aout, h);
	double deltaPh = arg(A[Nstop - 1] / A[0]);

	fprintf(file, "%15g,%15g,%15g,%15g,%15g,%15g,%15g", Aout, Pout, Amax, Pmax, Lmx, eff, deltaPh);

	if(shMemoryCreated) printAbsAtoSharedMemory(Nstop);

}
void TWT::printAbsAtoSharedMemory(int Nstop)
{
	double dz = Lmax / double(Nmax);
	do{
		if (!sharedMemory->lock()) continue;
		int nstop = Nstop;
		char *to = (char*)sharedMemory->data();
		char *from = (char*)A;
		int flag = 1;
		memcpy(to, &flag, sizeof(int));
		//	printf("flage = %i\n", flag);
		memcpy(to + sizeof(int), &nstop, sizeof(int));
		for (int i = 0; i < nstop; i++){
			*(double*)(to + 2 * sizeof(int) + i*sizeof(double)) = i*dz;
		}
		for (int i = 0; i < nstop; i++){
			*(double*)(to + 2 * sizeof(int) + (Nstop + i)*sizeof(double)) = abs(A[i]);
			//		printf("A = %g\n", abs(A[i]));
		}
		break;
	} while (true);
	sharedMemory->unlock();

}
void TWT::printCurrentParams(FILE *file)
{
	fprintf(file, "%g,%g,%g,%g,%g,%g,%g,%g,", Nperiods, synch_angle, Current_ampers, voltage, group_speed, Q, inputPower_watts, (syncwave!=NULL) ? syncwave->frequency() : 0);
}
void TWT::printParamsHeader(FILE *file)
{
	fprintf(file, "number of periods, synch angle[deg],Current[A],Voltage[kV],group beta,Q (one period),inpur power [W],input frequency [GHz],");
}
void TWT::printResultHeader(FILE *file)
{
	fprintf(file, "output amp, output power, max amp, max power, max length, efficiency, phase shift");
}
void TWT::printDataHeader(FILE *file)
{
	printParamsHeader(file);
	printResultHeader(file);
	fprintf(file, "\n");
}
void TWT::generateLongitudinalStructure(double h)
{
	if (longStructRealRe == NULL) return;
	Nperiods = longStructRealRe->xMax() / period;
	double La = Nperiods*period*h;
	double dz = Lmax / double(Nmax);
	int Nstop = ceil(La / dz);
	for (int i = 0; i < Nmax; i++)
	{
		double hz = i*dz;
		longitudinalStructureRe[i] = longStructRealRe->at(hz / h);
		longitudinalStructureIm[i] = longStructRealIm->at(hz / h);
	}

}
void TWT::generateQStructure(double h)
{
	if (qStructRe == NULL) return;
	if (longStructRealRe == NULL) Nperiods = longStructRealRe->xMax() / period;
	double La = Nperiods*period*h;
	double dz = Lmax / double(Nmax);
	int Nstop = ceil(La / dz);
	for (int i = 0; i < Nmax; i++)
	{
		double hz = i*dz;
		qStructure[i] = qStructRe->at(hz / h);
	}

}
double TWT::solveTWT()
{
	if (!initialized) {
		printf("Error:  Not initialized\n"); return -1;
	}


	if (syncwave != NULL)
	{
		syncwave->setFrequency(k1*299.8 / (2.*Pi), synch_angle);
		synch_angle = syncwave->wavenumber();
		group_speed = syncwave->groupSpeed();
	}

	double h = 2.*Pi / period*(synch_angle / 360.);
	double G = paramG(h);
	double inputAmp = amplitude(inputPower_watts, h);
	double lossKappa = k1 / h*0.5 / (Q*group_speed);


	if (longitudinalStructureRe) generateLongitudinalStructure(h);
	if (qStructure) generateQStructure(h);

	printDataHeader(results);
	if (iteratedParams.size() > 0)
	{
		iterate(0);
	}
	else
	{
		printf("power =%g, voltage = %g, synch_angle = %g, hA = %g, beta_gr = %g\n", inputPower_watts, voltage, synch_angle, h*inputAmp, group_speed);
		fflush(0);

//		double ar, ai;
		cplx am = solve(A, ar, ai, inputAmp, lossKappa, 0, structReal, structIm, G, 0, true, longitudinalStructureRe, longitudinalStructureIm, qStructure, mesh);
		printResults(results, A);
		printf("power = %g\n", power(abs(am), h));
/*		for (unsigned int i = 0; i < parasites.size(); i++)
		{
			double synch_angle = syncwave->wavenumber(voltage, i + 1);
			double fre = syncwave->frequency(voltage, i + 1);
			double speed = syncwave->groupSpeed(voltage, i + 1);
			printf("Considering parasite mode %i: \n", i);
			printf("\tsynchronous angle = %g, synchronous frequency = %g\n", synch_angle, fre);
			printf("\tstart current = %g\n", parasites[i].findStartCurrent(Current_ampers, -0.001, fre, synch_angle, speed, Q));
		}*/
	}
	return 1;

}
void TWT::iterate(int paramsInd)
{
	for (vec::iterator s = iteratedParams.at(paramsInd).begin(); s != iteratedParams.at(paramsInd).end(); ++s)
	{
		if (paramsInd != iteratedParams.size() - 1)
		{
			changeParam(iteratedNames[paramsInd], *s);
			iterate(paramsInd + 1);
		}
		else
		{
			changeParam(iteratedNames[paramsInd], *s);
			if (syncwave)
			{
				synch_angle = syncwave->wavenumber();
				group_speed = syncwave->groupSpeed();
				k1 = syncwave->frequency() * 2 * Pi / 299.8;
			}

			double h = 2.*Pi / period*(synch_angle / 360.);
			double G = paramG(h);
			double inputAmp = amplitude(inputPower_watts, h);
			double lossK = lossKappa(h);
			if (longitudinalStructureRe) generateLongitudinalStructure(h);
			if (qStructure) generateQStructure(h);

			printf("power =%g, voltage = %g, synch_angle = %g, input_amp = %g\n", inputPower_watts, voltage, synch_angle, inputAmp);
			fflush(0);

//			double ar, ai;
			solve(A, ar, ai, inputAmp, lossK, 0, structReal, structIm, G, 0, false, longitudinalStructureRe, longitudinalStructureIm, qStructure, mesh);
			printCurrentParams(results);
			printResults(results, A);
			fprintf(results, "\n");
		}
	}
}