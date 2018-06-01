#include <io.h>
#include "device.h"
#include "orotron.h"
#include "Interpolation.h"
#include "synchronism.h"
#include "cu_mult.h"
#include "xml_routines.h"





Device::Device()
{
	solverStatus = 0;
	fieldLoaded = 0;
	b = new cplx [1000];
	A_stat = 0;
	Np = NP; Nq = 6*NQ; Ns = 10*NS; Nv = 1; // умножительный расчёт (1, 4, 4, 4)
//	Np = NP; Nq = NQ; Ns = NS; Nv = 1; // умножительный расчёт (1, 4, 4, 4)
//	Np = NP; Nq = 4*NQ; Ns = 8*NS; Nv = 8;   // оротронный расчёт
//	Np = NP; Nq =  NQ; Ns = NS; Nv = 1;   // отладка -- один warp в grid
	k0 = 2.*Pi*86.667/299.8;
	s_stat = 0;
	monitor_par = 0;
	initialized = false;
	debug = true;

	
}
Device::Device(QDomDocument* doc)
{
	solverStatus = 0;
	fieldLoaded = 0;
	b = new cplx [1000];
	A_stat = 0;	
//	Np = NP; Nq = 2*2*4*NQ; Ns = 8*NS; Nv = 1; // умножительный расчёт (1, 4, 4, 4)
//	Np = NP; Nq = NQ; Ns = NS; Nv = 1; // умножительный расчёт (1, 4, 4, 4)
//	Np = NP; Nq = 4*NQ; Ns = 8*NS; Nv = 8;   // оротронный расчёт
	Np = NP; Nq = 4*NQ; Ns = 16*NS; Nv = 1;   // отладка -- один warp в grid
//	Np = NP; Nq = NQ; Ns = NS; Nv = 1;   // отладка -- один warp в grid
	k0 = 2.*Pi*86.667/299.8;
	s_stat = 0;
	monitor_par = 0;
	fieldFileName = new char[200]; fieldFileName[0] = '\0';
	dispersionFileName = new char[200]; dispersionFileName[0] = '\0';

	if(!setXMLEntry(doc, "MaxPointNumber", &Nmax)) {printf("Error: MaxPointNumber is not found\n"); return;}		 
	if(!setXMLEntry(doc, "MaxLength", &Lmax)) {printf("Error: MaxLength is not found\n"); return;}		 
	if(!setXMLEntry(doc, "groupSpeedCoeff", &grpSpeedCoeff)) {printf("Warning: groupSpeedCoeff is not found\n"); }		 

	char comment[500];
	char solverName[100];
	setXMLEntry(doc, "comment", comment);
	setXMLEntry(doc, "solverName", solverName);
	setXMLEntry(doc, "fieldFileName", fieldFileName);
	setXMLEntry(doc, "dispersionFileName", dispersionFileName);

	if (strcmp(solverName, "twt")!= 0) initSolver(Nmax, Lmax, grpSpeedCoeff, solverName);
	
	double norm1;
	QDomNode LFsection = doc->elementsByTagName("LFsection").item(0); 	
	if(LFsection.isNull()){ printf("Error: LFsection entry is not found\n"); return;}
	if(!setXMLEntry(&LFsection, "period", &period)) {printf("Error: LF period is not found\n"); return;}
	if (!setXMLEntry(doc, "periodNumbers", &Nperiods, &iteratedParams, &iteratedNames)) { printf("Error: Number of LF period is not found\n"); return; }      //период, число периодов первой сек (если вывод излучения адиабатический, то можно немного добавить длины, так как синус немного "выползает"
	if(!setXMLEntry(&LFsection, "periodNorma", &norm1)) {printf("Error: periodNorma of LF period is not found\n"); return;} ; //для прямоугольной гофрировки wall = 0.1: Norm = ; wall = 0.07 Norm = ?
	Norma = 0.5*double(Nperiods)*norm1;																						  //для синусовой гофрировки wall = 0.05: Norm = 33.8; wall = 0.07 Norm = 42;
	if(! setXMLEntry(&LFsection, "QFactor", &Q)) {printf("Error: QFactor of LF period is not found\n"); return;} ; //добротности первой (полная) секции
	if (!setXMLEntry(&LFsection, "structureWidth", &structureWidth)) { printf("Warning: structureWidth of LF period is not found\n"); }; //ширина структуры
	if (!setXMLEntry(&LFsection, "groupSpeed", &group_speed)) { printf("Warning: groupSpeed of LF period is not found\n"); }; //групповая скороть

	if (!setXMLEntry(doc, "current", &Current_ampers, &iteratedParams, &iteratedNames)) { printf("Error: beam current is not found\n"); return; }		//Ток
	k1 = 0;
	setXMLEntry(doc, "frequency", &k1, &iteratedParams, &iteratedNames);	k1 *= 2 * Pi / 299.8;    	//волновое число (частота возбуждения)
	if(!setXMLEntry(doc, "power", &inputPower_watts, &iteratedParams, &iteratedNames)) inputPower_watts = 0;					//входная мощность
	setXMLEntry(doc, "deltaFrequency", &delta_freq);				//	  и сдвиг частоты
	setXMLEntry(doc, "voltage", &voltage, &iteratedParams, &iteratedNames);							//	  напряжение			
 	setXMLEntry(doc, "beamStructureGap", &wall);					//расстояние от стенки и толщина пучка
	setXMLEntry(doc, "beamWidth", &beamWidth, &iteratedParams, &iteratedNames);					    //ширина пучка
	setXMLEntry(doc, "beamHeight", &beamHeight);					//высота пучка 
	setXMLEntry(doc, "angle", &synch_angle);					    //угол синхронизма
    if (!setXMLEntry(doc, "fileIndex", &fileIndex)) fileIndex = 0;	//
	double Time, initAmpMax;
	setXMLEntry(doc, "Time", &Time);								// Время счёта

	setXMLEntry(doc, "initialAmplitudeMaximum", &initAmpMax);		//Максимум амплитуды в начальный момент
	initialized = true;

	if (_access(dispersionFileName, 0) == 0)
		syncwave = new Synchronism(dispersionFileName, period, 299.8*k1 / (2.*Pi));	

	QDomNode solversection = doc->elementsByTagName("solver").item(0);
	if (solversection.isNull()){
		printf("Error: solver entry is not found\n"); 
		return;
	}
	else
	{
		QDomNode calcgridsection = solversection.namedItem("GPUcalculatingGrid");
		if (calcgridsection.isNull())
		{
			Np = NP; Nq = 4 * NQ; Ns = 16 * NS; Nv = 1;   // отладка -- один warp в grid
			printf("GPU calculating grid is not found, using following grid: (%i, %i, %i, %i)\n", Np, Nq, Nv, Ns);
		}
		else
		{
			if (calcgridsection.childNodes().size() < 4){
				Np = NP; Nq = 4 * NQ; Ns = 16 * NS; Nv = 1;   // отладка -- один warp в grid
				printf("Warning: GPU calculating grid size is less than 4, using following grid: (%i, %i, %i, %i)\n", Np, Nq, Nv, Ns);
			}
			else
			{
				int Nptmp, Nqtmp, Nstmp, Nvtmp;
				QDomNodeList children = calcgridsection.childNodes();
				Nptmp = children.at(0).toElement().text().toInt();
				Nqtmp = children.at(1).toElement().text().toInt();
				Nstmp = children.at(2).toElement().text().toInt();
				Nvtmp = children.at(3).toElement().text().toInt();

				Np = (Nptmp / NP)*NP;
				Nq = (Nqtmp / NQ)*NQ;
				Ns = (Nstmp / NS)*NS;
				Nv = Nvtmp;
				if (Np < NP) { printf("Warning: wrong grid dimension 1\n"); Np = NP; }
				if (Nq < NQ) { printf("Warning: wrong grid dimension 2\n"); Nq = 4*NQ; }
				if (Ns < NS) { printf("Warning: wrong grid dimension 3\n"); Ns = 16*NS; }
				if (Nv < 1)  { printf("Warning: wrong grid dimension 4\n"); Nv = 1; }
				printf("Using following GPU calculating grid: (%i, %i, %i, %i)\n", Np, Nq, Nv, Ns);
			}
		}

	}

	// ........ opening output file ............................
	char filename[200];
	int attr = 0;
	setXMLattribute(doc, "fileIndex", "append", &attr);
	sprintf(filename, "result_%i.dat", fileIndex);

	char write[2] = "w";
	char writeAppend[3] = "a+";
	if (attr == 0)
		results = fopen(filename, write);
	else
		results = fopen(filename, writeAppend);

	// ........................................................

}



double Multiplier::findAstat(double nextA, int N_it, double G, double *resS)
{
	double A0, A1;
	double F0, F1;

	double xi = 0.000001;
	double omega = k1*(10.*c);

	double delta1, delta2, delta;
	double svF1, svF2;

																	/////////
	double inputPower_dimensionless = inputPower_watts*1.E7/(omega/1.*(m*c*c/e)*(m*c*c/e)*(Current_ampers/17045.)*(1./(k1*10.))); // = S a
																////////
	
	printf("Qa = %g,\t G = %g,\t inputPowerDimensionless = %g \n", Q, G, inputPower_dimensionless);
	
//	FILE *debfile = fopen("d:\\Piotr\\w_orotron_Data\\debug.txt", "w");
//	printf("The Delta Energy  =  %g\n", ElectronsDeltaEnergy(1e-5));
	
	for(int i = 0; i < N_it; i++)
	{
	/*	A0 = double(i)/double(N_it)*0.005;
		fprintf(debfile, "%g\t%g\n", A0, DeltaEnergy(A0));*/
		
		A0 = nextA;
		F0 = A0*A0 + 2.*Q*G*( -inputPower_dimensionless + ElectronsDeltaEnergy(A0) ); 
	
		A1 = nextA*(1.+xi);
		F1 = A1*A1 + 2.*Q*G*( -inputPower_dimensionless + ElectronsDeltaEnergy(A1) );

		if(A0*A1 < 0) {xi *= 0.5; continue;}
		if(fabs(F1 - F0) < 1e-30) break;
		if((i>2) &&((fabs(delta) < 1e-20)||fabs(svF2 - svF1) < 1e-20)) break;

		delta = - F0/(F1- F0)*(A1 - A0);

		if(i == 0) { delta1 = delta; svF1 = F0;}
		if(i == 1) { delta2 = delta; svF2 = F0;}
		if(i > 1)  { delta1 = delta2; delta2 = delta; svF1 = svF2; svF2 = F0;}

		/*if(i > 1)*/ printf("\tnextA  = %g\t delta = %g\tRes. = %g\tF0 = %g\tdE = %g\tA0 = %g\n", nextA, fabs(delta), fabs(svF2 - svF1), F0, ElectronsDeltaEnergy(A0), A0);
		
		nextA = A0 + delta;

		
	
					
	}
	printf("A_{stat}  = %g\n", nextA);
	printf(" 2 q G = %g \n", 2.*Q*G);
	if(resS != 0) *resS = inputPower_dimensionless/nextA;
	//fclose(debfile);

	return nextA;

}
void Device::setCurrent(double curr)
{
	Current_ampers = curr;
}
double Device::paramG()
{
	return Current_ampers / 17045.*1. / (k1*Norma); //k1 в мм^{-1}, Norma в мм^3.
}
double Device::power(double kpd)
{
	return voltage*Current_ampers*kpd;

}
void Device::changeParam(string name, double par)
{
	if (name == "voltage") {
		voltage = par;
		if (syncwave != NULL){
			syncwave->setVoltage(par);
		}
	}
	if (name == "QFactor") Q = par;
	if (name == "power") inputPower_watts = par;
	if (name == "current") Current_ampers = par;
	if (name == "frequency"){
		if(syncwave != NULL) syncwave->setFrequency(par);
	}
	if (name == "periodNumbers") {
		Nperiods = par;
		if (Nperiods*period > Lmax) {
			Nperiods = 0.99*Lmax / period;
		}
	}
}