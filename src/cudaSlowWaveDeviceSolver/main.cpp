#include"device.h"
#include <QtCore\qsharedmemory.h>
#include <QtCore\qbuffer.h>
#include <QtCore\qdatastream.h>
#include<iostream>
#include"twt_2d.h"
#include"twt_1d.h"
#include"twt_0d.h"
#include"bwo.h"
#include "orotron.h"
#include "orotron_multimode.h"
#include "multiplier_spcharge_2d.h"
#include "Interpolation.h"
#include "synchronism.h"
#include<mpi.h>
#include<fftw3.h>
#include<direct.h>
#include<time.h>
#include<sys\timeb.h>
#include<Windows.h>



#include"cu_mult.h"
#include "xml_routines.h"


//const cplx I(0, 1);
/*int N = 100;//20000;
double dt = 400/double(N);
int Nalpha = 1;
int Ndelta = 100000;
int Nq = 32;
int Nx = 10;
double dDelta = 0.1/Ndelta;
double deltaAlpha = 0.2/double(Nalpha);
double dx = 4./double(Nx);
double xGr = 4;

bool debflag  = false;


double* localField1 = new double [N];
double* localField2 = new double [N];
double* localField3 = new double [N];
double* localField4 = new double [N];
double* localField5 = new double [N];*/



//const double PI = 3.14159265358979323846264338328;
//const double Pi = 3.14159265358979323846264338328;

#define QX Ind(q0, x0, Nx)

//Индекс 1 после функций или его отсутствие расчитывает случай бесконечно толстого пучка
//Индекс 2 - случай для пучка конечной толщины - конец пучка - конец прибора
//Индекс 3 - конечная толщина, задаваемая произвольно
inline int Ind(int i, int j, int N) {return (i*N + j);}
void multiplier()
{

	char rootFolder[600], filename[300], file2name[300], file3name[100];
	_mkdir("MultiplierTest"); _chdir("MultiplierTest");
	_mkdir("Data1"); _chdir("Data1"); _getcwd(rootFolder, 600);

	//	_mkdir("D:\\Piotr\\Multiplier_260GHz_Data");



// Параметры (ВСЕ) 
	double current = 0.3;				//Ток
	double period = 0.84; int Nper = 21; 	//период, число периодов первой сек
    double ld = 11., lb = 15.;				//длина дрейфовой секции, длина второй секции
	double ld1 = 4., ld2 = 7., la2 = 4;
	double kw = 2.*Pi*86.667/299.8;		//волновое число (частота возбеждения)
	double norm1 = 3.018;				//норма одного периода первой секции (расстояние от стенки до пучка 0,07: 3.018; расстояние от стенки до пучка 0,1: 4.443)
	double norm1B= 3.848;			//норма одного периода второй секции (расстояние от стенки до пучка 0,07: 3.848; расстояние от стенки до пучка 0,1: 13.6380622)
	double p_in = 0.8, delta_freq =  0.; //входная мощность, и сдвиг частоты
	double voltage = 15.9;			    //напряжение
	double Qa = 450, Qb = 1200;			//добротности первой (полная) и второй (омическая) секций
	double Qa1 = Qa, Qa2 = 500;
	double wall = 0.07;					//расстояние от стенки и толщина пучка
	double spchQ1 = 6.83e-5;			//амплитуда пространственного заряда на НЧ; для 1d 3.46e-5; для 2d 2.2e-5
	double spchQ3 = 6.83e-5;		  	   //амплитуда пространственного заряда на ВЧ;  для 1d 5e-5; для 2d  4.23e-5
	char solverName[] = "multiplier_spcharge_2d";		//имя солвера "multiplier", "multiplier_spcharge", "multiplier_spcharge_2d"
//	Multiplier_SpCharge_2D sol;  TODO закоментировал - нет пустого конструктора, поправить, позже
	char difrFlag[] = "Empiric"; //выбор зависимости дифр. добротности от длины: "Ohm","Min", "Empiric"
	char comment[] = "До 615 в FindBStatDetunded вместо полной стояла\n омическая пересчёт для полной добротности равной половине омической\n Тоже, что в 631, ширина полосы больше плюс смещение по длине";

//	sol.initMultiplierSolver(3000, 110, 0., solverName);

	int file_index = 5471;				//Индекс файла с результатами

	sprintf(file3name,"P_vs_Np_Lb_Ld");
	sprintf(filename,  "F:\\Piotr\\Multiplier_260GHz_Data\\%s_%i.txt", file3name, file_index);
	sprintf(file2name, "F:\\Piotr\\Multiplier_260GHz_Data\\%s_params_%i.txt", file3name, file_index);

	FILE *file =  fopen(filename, "w");

	int Nj = 1, Ni = 1;

	char file1[300];
	char file2[300];


    double inputPow; 
//	for(voltage = 15.6; voltage <= 16.; voltage += 0.2)
//	for(int Nper = 11; Nper < 21; Nper+=2)
	for(int j = 0; j < Nj; j++)
	for(int i = 0; i < Ni; i++)
	for(int u = 0; u < 1; u++)
	{///////
		
		lb = 14;// + 1. +  32.*double(i)/double(Ni);
		ld = 3.5;//2. +  30.*3./11.;//double(j)/double(Nj);
		delta_freq = 0.253;
//		delta_freq = 0.23 + (0.28-0.23)/15.*(double)u;// + 0.018/double(Ni)*i;;//0.195 + (0.245-0.195)/25.*(double)u;
		
	/*	double output_power = sol.getHFoutputPower( current,		 //Ток
													period, Nper,	 //период, число периодов первой секции
													ld, lb,			 //длина дрейфовой секции, длина второй секции
													kw,				 //волновое число
													norm1,			 //норма одного периода первой секции
													norm1B,			 //норма одного периода первой секции
													voltage,	  	 //напряжение
													p_in, delta_freq,//входная мощность, и сдвиг частоты
													Qa, Qb,			 //добротности первой и второй секций
													wall,			 //расстояние  от  стенки и толщина пучка
													spchQ1, spchQ3,	 //амплитуда пространственного заряда на НЧ и ВЧ частоте соответственно
													file2name,		 //имя файла для печати параметров
													comment,			//коментарий в файл параметров
													&inputPow, 
													difrFlag);		 */

	/*	double output_power = sol.getHFoutputPowerDoubleScheme( current,		 //Ток
													period, Nper,	 //период, число периодов первой секции
													ld1, ld2, la2, lb,//длины секций
													kw,				 //волновое число
													norm1,			 //норма одного периода первой секции
													norm1B,			 //норма одного периода первой секции
													voltage,	  	 //напряжение
													p_in, delta_freq,//входная мощность, и сдвиг частоты
													Qa1, Qa2, Qb,	 //добротности первой, второй и выходной секции
													wall,			 //расстояние от стенки и толщина пучка
													file2name);		 //имя файла для печати параметров*/

	/*	double output_power = sol.getHFoutputPowerMultiModes( current,		 //Ток
													period, Nper,	 //период, число периодов первой секции
													ld, lb,			 //длины секций
													kw,				 //волновое число
													norm1,			 //норма одного периода первой секции
													norm1B,			 //норма одного периода первой секции
													voltage,	  	 //напряжение
													p_in, delta_freq,//входная мощность, и сдвиг частоты
													Qa1, Qa2, Qb,	 //добротности первой, второй и выходной секции
													wall,			 //расстояние от стенки и толщина пучка
													file2name,		 //имя файла для печати параметров
													comment);*/
 

//		fprintf(file, "%g\t%g\t%g\n", delta_freq, output_power, inputPow);
//		fprintf(file, "%g\t%g\t%g\n", lb, ld, output_power);
// TODO		fprintf(file, "%g\t%g\t%g\t%g\t%g\n", lb, ld, delta_freq, output_power, inputPow);
		fflush(file);
//		fprintf(file, "%g\t%i\t%g\t%g\t%g\n", voltage, Nper, ld2, lb, output_power);

		sprintf(file1, "F:\\Piotr\\Multiplier_260GHz_Data\\currentData\\field1_%i_%i.txt", file_index, i);
		sprintf(file2, "F:\\Piotr\\Multiplier_260GHz_Data\\currentData\\field3_%i_%i.txt", file_index, i);
//TODO		sol.PrintCurrentsSpaceCharge2D(file1, file2);

	}
	fclose(file);



//TODO	sol.releaseMultiplierMemory();
	cudaDeviceReset();
	
}
/* TODO default contstuctor
void orotron()
{

	char rootFolder[600], filename[300], file2name[300], file3name[100];
//	_mkdir("MultiplierTest"); 
//	_chdir("MultiplierTest");
//	_mkdir("Data1"); 
//	_chdir("Data1"); 
	_getcwd(rootFolder, 600);

	Orotron sol; 

	sol.initMultiplierSolver(1000, 40, 0., "orotron");


//	_mkdir("D:\\Piotr\\Multiplier_260GHz_Data");



// Параметры (ВСЕ)
	double current = 1.0;				//Ток
	double period = 1; int Nper = 21;	//период, число периодов первой сек (если вывод излучения адиабатический, то можно немного добавить длины, так как синус немного "выползает"
    double ld = 1., lb = 1.;				//длина дрейфовой секции, длина второй секции
	double kw = 2.*Pi*96.58/299.8;		//волновое число (частота возбуждения)
	double norm1 = 33.8;		  		//для прямоугольной гофрировки wall = 0.1: Norm = ; wall = 0.07 Norm = ?
										//для синусовой гофрировки wall = 0.05: Norm = 33.8; wall = 0.07 Norm = 42;
	double norm1B= 13.6380622;			//
	double p_in = 0., delta_freq =  0.; //входная мощность, и сдвиг частоты
	double voltage = 30.5;				//напряжение
	double Qa = 450, Qb = 1200;   		//добротности первой (полная) и второй (омическая) секций
	double wall = 0.06;					//расстояние от стенки и толщина пучка
	char comment[] = "гофрировка - синус; стартовый ток";

//	devID = findCudaDevice(argc, (const char **)argv);
//  cudaGetDeviceProperties(&deviceProp, devID);


	int file_index = 1;				//Индекс файла с результатами

	sprintf(file3name,"oroStartCurr");
	sprintf(filename,  "F:\\Piotr\\w_orotron_Data\\%s_%i.txt", file3name, file_index);

	sprintf(file2name, "F:\\Piotr\\w_orotron_Data\\%s_params_%i.txt", file3name, file_index);

	FILE *file =  fopen(filename, "w");

	int Nj = 1, Ni = 75;

	double *a00 = new double [Ni*Nj];

	double Qdifr = 2000; //32.*(Nper*Nper*pow(period, 2))/9.;

	Qa = 1./(1./2000. + 1./Qdifr);
	printf("Qa = %g\n", Qa);

//	Qa = 450;


//	for(voltage = 15.6; voltage <= 16.; voltage += 0.2)
//	for(int Nper = 11; Nper < 21; Nper+=2)
	for(int j = 0; j < Nj; j++)
    for(int i = 0; i < Ni; i++)
	{
		voltage =  29.3+ 2.5*double(i)/double(Ni);
	//	voltage =  15.7+ 2*double(i)/double(Ni);
	//	current =  1.2;//0.3 + 2.5*double(j)/double(Nj);
//		if(j == 0) a0 = 0; else a0 = a00[Ni*(j-1)+i];
//		kw = 2.*Pi*86.667/299.8*(1. +  (double(i) - 15.)/15.*0.0015);
		//
		double kpd				= sol.orotronStartCurrent//sol.solveOrotron
													( current,		 //Ток
													period, Nper,	 //период, число периодов первой секции
													ld, lb,			 //длина дрейфовой секции, длина второй секции
													kw,				 //волновое число
													norm1,			 //норма одного периода первой секции
													norm1B,			 //норма одного периода первой секции
													voltage,	  	 //напряжение
													0.  , delta_freq,//входная мощность, и сдвиг частоты
													Qa, Qb,			 //добротности первой и второй секций
													wall,			 //расстояние  от  стенки и толщина пучка
													file2name,	  	 //имя файла для печати параметров
													comment);	 //коментарий и начальная точка отсчёта	 
		a00[Ni*j + i] = kpd;


		double Norma= norm1*0.5*double(Nper);
		double G = current/17045.*1./(kw*Norma);
		double gamma = 1.+ voltage/511.;
		double nu = 1./(gamma*(gamma*gamma -1.));
		double beta = sqrt(gamma*gamma -1)/gamma;
		double Delta =  1 - kw/(2.*Pi/period*beta);
		double h = 2*Pi/period;

//		fprintf(file, "%g\t%g\t%g\n", current, voltage, -/*-Qa/Qdifr*//*voltage*current*kpd);
		fprintf(file, "%g\t%g\t\n", voltage, kpd);
//		Безразмерные величины (нужно значение a_stat)
//		fprintf(file, "%g\t%g\t%g\t%g\n", sqrt(fabs(astat)*nu*h)*period*double(Nper), (fabs(astat) > 1e-5)? sqrt(h)*Delta/sqrt(fabs(astat)*nu): 0, pow(fabs((astat)), 1.5)*sqrt(h*nu)/(2.*Qa*G), astat);
//		fprintf(file, "%g\t%i\t%g\t%g\t%g\n", voltage, Nper, ld2, lb, output_power);
//		fprintf(file, "%g\t%g\n", 3.*299.8*kw/(2.*Pi), output_power);
	}
	fclose(file);



	sol.releaseMultiplierMemory();
	cudaDeviceReset();
	
}*/
void twt()
{
	_chdir("F:\\Piotr\\CalcData\\twt_data");
	char rootFolder[600];
	_mkdir("MultiplierTest"); //_chdir("MultiplierTest");
	_mkdir("Data1");// _chdir("Data1"); 
	_getcwd(rootFolder, 600);

	QDomDocument doc("inputParameters");
	QString err_msg; int err_line; int err_column;
	QFile file("input.xml");
	if (!file.open(QIODevice::ReadOnly))
		return;
	bool set_cont = doc.setContent(&file, &err_msg, &err_line, &err_column);
	if (!set_cont)
	{
		cout << qPrintable(err_msg) << "; line number" << err_line << ",  column number" << err_column << "\n";
		file.close();
		return;
	}
	file.close();

	TWT_2D sol(&doc);

	sol.solveTWT();


}
void bwo()
{
	_chdir("F:\\Piotr\\CalcData\\twt_data\\W_bwo_FoldWg_3"); 

	char rootFolder[600];
	_getcwd(rootFolder, 600);

	QDomDocument doc("inputParameters");
	QString err_msg; int err_line; int err_column;
	QFile file("input.xml");
	if (!file.open(QIODevice::ReadOnly))
	{
		printf("\nError reading input file\n");
		return;
	}
	bool set_cont = doc.setContent(&file, &err_msg, &err_line, &err_column);
	if (!set_cont)
	{
		cout << qPrintable(err_msg) << "; line number" << err_line << ",  column number" << err_column << "\n";
		file.close();
		return;
	}
	file.close();

	BWO<TWT_2D> sol(&doc);

/*	double delta = 0.3;
	double curr = sol.currentFromNormalizedLength(2.5, &delta);
	printf("values from norm = %g, %g\n", curr, delta);
	curr = 0.249; delta = 0.14;
	double stCurr = sol.findStartCurrent(4*curr, &delta);
	printf("start current %g;\n", stCurr);
	printf("delta = %g\n", delta);*/
//	sol.getNormalizedLength(stCurr, &delta);
    
//	sol.solveBWO(0.00016, -0.02);
	
	sol.solveBWO();


}
int startXML(int argc, char *argv[])
{

	char rootFolder[600];
	_getcwd(rootFolder, 600);
	QDomDocument doc("inputParameters");

//	if (true)
	if (argc > 1 && strcmp(argv[1], "control") == 0)
	{
		printf("Reading shared memory...\n");
		QSharedMemory sharedMemory("slowwavedevicecontroller");
		if (!sharedMemory.attach())	{
			printf("Error ataching shared memory\n"); return -3;
		}
		bool att = sharedMemory.isAttached();

		QBuffer buffer;
		QString errmsg; int x, y;
		QString docstring;
		
		
		QDataStream stream(&buffer);
		char *mem;

		
		sharedMemory.lock();
		mem = (char*)sharedMemory.constData();
		buffer.setData(mem, sharedMemory.size());
		buffer.open(QBuffer::ReadOnly);
		stream >> docstring;
		sharedMemory.unlock();
		if(!doc.setContent(docstring, &errmsg, &x, &y)) printf("Error: %s: %i,%i ", qPrintable(errmsg), x, y );
		if (setXMLEntry(&doc, "workingDirectory", rootFolder)) _chdir(rootFolder);
	}
	else
	{		
		QString err_msg; int err_line; int err_column;
		if (argc > 2 && strcmp(argv[1], "-i") == 0)
		{
			try   { 
				_chdir(argv[2]); 
			}
			catch (std::exception  e){
				printf("Failed to find specified directory\n");
				return 3;
			}
		}
		else 
			_chdir("F:\\Piotr\\CalcData\\verifiedCpy\\95GHz_bwo_axial_twt0d_complexClinotron");

		QFile file("input.xml");
		if (!file.open(QIODevice::ReadOnly))
		{
			printf("\nError reading input file\n");
			return 3;
		}
		bool set_cont = doc.setContent(&file, &err_msg, &err_line, &err_column);
		if (!set_cont)
		{
			cout << qPrintable(err_msg) << "; line number" << err_line << ",  column number" << err_column << "\n";
			file.close();
			return 3;
		}
		file.close();
	}

	

	char problemType[30], solverName[30];
	if (!setXMLEntry(&doc, "problemType", problemType)) { printf("Error finding problem type in input\nStop\n"); return 1; }
	if (!setXMLEntry(&doc, "solverName", solverName)) { solverName[0] = '\0';  sprintf(solverName, "twt"); }
	

	timeb start, finish;
	ftime(&start);
	bool solverfound = false;
	if (strcmp(problemType, "bwo") == 0)
	{
		if (strcmp(solverName, "twt") == 0)
		{
			BWO<TWT_2D> sol(&doc);
			sol.solveBWO();
		}
		else if (strcmp(solverName, "twt1d") == 0)
		{
			BWO<TWT_1D> sol(&doc);
			sol.solveBWO();
		}
		else if (strcmp(solverName, "twt0d") == 0)
		{
			BWO<TWT_0D> sol(&doc);
			sol.solveBWO();
		}
		else
		{
			printf("solver name is not recognized\n");
		}

		solverfound = true;
	}
	if (strcmp(problemType, "bwostart") == 0)
	{
		if (strcmp(solverName, "twt") == 0)
		{
			BWO<TWT_2D> sol(&doc);
			sol.findStartCurrent();
		}
		else if (strcmp(solverName, "twt1d") == 0)
		{
			BWO<TWT_1D> sol(&doc);
			sol.findStartCurrent();
		}
		else if (strcmp(solverName, "twt0d") == 0)
		{
			BWO<TWT_0D> sol(&doc);
			sol.findStartCurrent();
		}
		else
		{
			printf("solver name is not recognized\n");
		}
		solverfound = true;
	}
	if (strcmp(problemType, "twt") == 0)
	{
		TWT_2D sol(&doc);
		sol.solveTWT();
		solverfound = true;
	}
	if (strcmp(problemType, "twt1d") == 0)
	{
		TWT_1D sol(&doc);
		sol.solveTWT();
		solverfound = true;
	}
	if (strcmp(problemType, "twt0d") == 0)
	{
		TWT_0D sol(&doc);
		sol.solveTWT();
		solverfound = true;
	}
	if (strcmp(problemType, "orotron") == 0)
	{
		Orotron sol(&doc);
		sol.solveOrotron();
		solverfound = true;
	}
	if (strcmp(problemType, "orotronstart") == 0)
	{
		Orotron sol(&doc);
		sol.orotronStartCurrent();
		solverfound = true;
	}

	if (solverfound)
	{
		ftime(&finish);
		printf("Elapsed time %g s\n", 0.001*(float)(1000 * (finish.time - start.time) + finish.millitm - start.millitm));
		Sleep(1000);
		return 0;
	}
	else
	{
		printf("Error in recognizing problem type in input\nStop\n");
		return 2;
	}


}
void orotron_XML()
{

	char rootFolder[600], filename[300], file2name[300], *file3name;
	_getcwd(rootFolder, 600);

	 QDomDocument doc("inputParameters");
	 QString err_msg; int err_line; int err_column;
	 QFile file("F:\\Piotr\\CalcData\\mm_orotron_Data\\input.xml");
	 if (!file.open(QIODevice::ReadOnly))
		return ;
	 bool set_cont =  doc.setContent(&file, &err_msg, &err_line, &err_column);
	 if (!set_cont) 
	 {
		 cout<<qPrintable(err_msg)<<"; line number"<<err_line<<",  column number"<<err_column<<"\n";
		 file.close();
      return ;
	}
	file.close();

	Orotron sol(&doc);
	int Nmax = doc.elementsByTagName("MaxPointNumber").item(0).toElement().text().toInt();
	double Lmax = doc.elementsByTagName("MaxLength").item(0).toElement().text().toDouble();
	double grpSpeedCoeff = doc.elementsByTagName("groupSpeedCoeff").item(0).toElement().text().toDouble();
	QString tmp = doc.elementsByTagName("comment").item(0).toElement().text();
	QByteArray btmp  = tmp.toLocal8Bit();
	char* solcomment; solcomment  = (char*)btmp.data();
	QByteArray btmp1 = ((doc.elementsByTagName("solverName").item(0).toElement().text()).toLocal8Bit());
	char* solverName; solverName =  (char*)btmp1.data();	

	QDomNodeList ndlist1  = doc.elementsByTagName("problemName"); 
	QByteArray btmp2;
	if(ndlist1.length() != 0) 
	{
		btmp2 =  (ndlist1.item(0).toElement().text()).toLocal8Bit();
		file3name	    =  (char*)btmp2.data();
	}
	else
	{
		file3name = new char [100];
		sprintf(file3name, "unnamedProblem");
	}
	
	sol.initSolver(Nmax, Lmax, grpSpeedCoeff, solverName);


//	_mkdir("D:\\Piotr\\Multiplier_260GHz_Data");



// Параметры (ВСЕ)
	double current = doc.elementsByTagName("current").item(0).toElement().text().toDouble();				//Ток
	QDomNode LFsection = doc.elementsByTagName("LFsection").item(0);
	if(LFsection.isNull()){ printf("Error: LFsection entry is not found\n"); return;}
	double period = LFsection.namedItem("period").toElement().text().toDouble();
	int    Nper   = LFsection.namedItem("periodNumbers").toElement().text().toInt();	//период, число периодов первой сек (если вывод излучения адиабатический, то можно немного добавить длины, так как синус немного "выползает"
	double norm1  = LFsection.namedItem("periodNorma").toElement().text().toDouble();	//для прямоугольной гофрировки wall = 0.1: Norm = ; wall = 0.07 Norm = ?
																						//для синусовой гофрировки wall = 0.05: Norm = 33.8; wall = 0.07 Norm = 42;
	double Qa =     LFsection.namedItem("QFactor").toElement().text().toDouble();       //добротности первой (полная) секции
	

	QDomNode HFsection = doc.elementsByTagName("HFsection").item(0);
	if(HFsection.isNull()) {printf("Error: LFsection entry is not found\n"); return;}
	double norm1B=  HFsection.namedItem("periodNorma").toElement().text().toDouble(); 
	double lb  =    HFsection.namedItem("length").toElement().text().toDouble();        //длина длина второй секции
	double  Qb =    HFsection.namedItem("QFactor").toElement().text().toDouble();		// добротности второй (омическая) секции

    double ld =	 doc.elementsByTagName("Dreif").item(0).toElement().text().toDouble(); //длина дрейфовой секции, длина второй секции
	double kw =  2.*Pi/299.8*doc.elementsByTagName("frequency").item(0).toElement().text().toDouble();//волновое число (частота возбуждения)

	
	double p_in = doc.elementsByTagName("power").item(0).toElement().text().toDouble();									//входная мощность, и сдвиг частоты
	double delta_freq =  doc.elementsByTagName("deltaFrequency").item(0).toElement().text().toDouble();
	double voltage =     doc.elementsByTagName("voltage").item(0).toElement().text().toDouble();						//напряжение
	
	double wall = doc.elementsByTagName("beamStructureGap").item(0).toElement().text().toDouble();						//расстояние от стенки и толщина пучка
	int Namm = doc.elementsByTagName("numberOfLongitudinalHarmonics").item(0).toElement().text().toDouble();																											//Количиество учитываемых продольных мод
	double Time = doc.elementsByTagName("Time").item(0).toElement().text().toDouble();									// Время счёта
	double initAmpMax = doc.elementsByTagName("initialAmplitudeMaximum").item(0).toElement().text().toDouble();			//Максимум амплитуды в начальный момент
	
	char comment[] = "стартовый ток";

//	devID = findCudaDevice(argc, (const char **)argv);
//  cudaGetDeviceProperties(&deviceProp, devID);


	int file_index = 14;				//Индекс файла с результатами

	sprintf(file3name,"LFsectionStartCurr");
	sprintf(filename,  "F:\\Piotr\\Multiplier_260GHz_Data\\%s_%i.txt", file3name, file_index);

	sprintf(file2name, "F:\\Piotr\\Multiplier_260GHz_Data\\%s_params_%i.txt", file3name, file_index);

	FILE *fileR =  fopen(filename, "w");

	int Nj = 1, Ni = 75;

	double *a00 = new double [Ni*Nj];

	double Qdifr = 32.*(Nper*Nper*pow(period, 2))/9.;

	Qa = 1./(1./2000. + 1./Qdifr);
	printf("Qa = %g\n", Qa);

	Qa = 450;


//	for(voltage = 15.6; voltage <= 16.; voltage += 0.2)
//	for(int Nper = 11; Nper < 21; Nper+=2)
	for(int j = 0; j < Nj; j++)
    for(int i = 0; i < Ni; i++)
	{
	//	voltage =  29.5+ 2.5*double(i)/double(Ni);
		voltage =  15.7+ 2*double(i)/double(Ni);
	//	current =  1.2;//0.3 + 2.5*double(j)/double(Nj);
//		if(j == 0) a0 = 0; else a0 = a00[Ni*(j-1)+i];
//		kw = 2.*Pi*86.667/299.8*(1. +  (double(i) - 15.)/15.*0.0015);
		//возращает электронный КПД
		double kpd				= sol.orotronStartCurrent//sol.solveOrotron
													( current,		 //Ток
													period, Nper,	 //период, число периодов первой секции
													ld, lb,			 //длина дрейфовой секции, длина второй секции
													kw,				 //волновое число
													norm1,			 //норма одного периода первой секции
													norm1B,			 //норма одного периода первой секции
													voltage,	  	 //напряжение
													0.  , delta_freq,//входная мощность, и сдвиг частоты
													Qa, Qb,			 //добротности первой и второй секций
													wall,			 //расстояние  от  стенки и толщина пучка
													file2name,	  	 //имя файла для печати параметров*/
													comment);	 //коментарий и начальная точка отсчёта	 
		a00[Ni*j + i] = kpd;


		double Norma= norm1*0.5*double(Nper);
		double G = current/17045.*1./(kw*Norma);
		double gamma = 1.+ voltage/511.;
		double nu = 1./(gamma*(gamma*gamma -1.));
		double beta = sqrt(gamma*gamma -1)/gamma;
		double Delta =  1 - kw/(2.*Pi/period*beta);
		double h = 2*Pi/period;

//		fprintf(file, "%g\t%g\t%g\n", current, voltage, -/*-Qa/Qdifr*/voltage*current*kpd);
		fprintf(fileR, "%g\t%g\t\n", voltage, kpd);
//		Безразмерные величины (нужно значение a_stat)
//		fprintf(file, "%g\t%g\t%g\t%g\n", sqrt(fabs(astat)*nu*h)*period*double(Nper), (fabs(astat) > 1e-5)? sqrt(h)*Delta/sqrt(fabs(astat)*nu): 0, pow(fabs((astat)), 1.5)*sqrt(h*nu)/(2.*Qa*G), astat);
//		fprintf(file, "%g\t%i\t%g\t%g\t%g\n", voltage, Nper, ld2, lb, output_power);
//		fprintf(file, "%g\t%g\n", 3.*299.8*kw/(2.*Pi), output_power);
	}
	fclose(fileR);



	sol.releaseDeviceMemory();
	cudaDeviceReset();
	
}
/*void multimode_orotron() TODO: solve problem with default constuctor
{

	char rootFolder[600], filename[300], file2name[300], file3name[100];
	_mkdir("MultiplierTest"); _chdir("MultiplierTest");
	_mkdir("Data1"); _chdir("Data1"); _getcwd(rootFolder, 600);

	Orotron sol;

	sol.initMultiplierSolver(1000, 40, 0., "orotron");


//	_mkdir("D:\\Piotr\\Multiplier_260GHz_Data");



// Параметры (ВСЕ)
	double current = 0.01;				//Ток
	double period = 3.1; int Nper = 10.;	//период, число периодов первой сек (если вывод излучения адиабатический, то можно немного добавить длины, так как синус немного "выползает"
    double ld = 15., lb = 15.;				//длина дрейфовой секции, длина второй секции
	double kw = 2.*Pi*35.2/299.8;		//волновое число (частота возбуждения)
	double norm1 = 143.3;					//для прямоугольной гофрировки wall = 0.1: Norm = ; wall = 0.07 Norm = ?
										//для синусовой гофрировки wall = 0.05: Norm = 33.8; wall = 0.07 Norm = 42;
	double norm1B= 13.6380622;			//
	double p_in = 0., delta_freq =  0.; //входная мощность, и сдвиг частоты
	double voltage = 200.;				//напряжение
	double Qa = 2539, Qb = 1200;		//добротности первой (полная) и второй (омическая) секций
	double wall = 0.0001;					//расстояние от стенки и толщина пучка
	int Namm = 3;					    //Количиество учитываемых продольных мод
	double Time = 30.;					// Время счёта
	double initAmpMax = 1e-6;			//Максимум амплитуды в начальный момент
	char comment[] = "Первый прогон";



	

//	devID = findCudaDevice(argc, (const char **)argv);
//  cudaGetDeviceProperties(&deviceProp, devID);


	int file_index = 1;				//Индекс файла с результатами



	sprintf(file3name,"P_vs_current_U_Ld");
	sprintf(filename,  "F:\\Piotr\\CalcData\\mm_orotron_Data\\%s_%i.txt", file3name, file_index);

	sprintf(file2name, "F:\\Piotr\\CalcData\\mm_orotron_Data\\%s_params_%i.txt", file3name, file_index);

//	FILE *file =  fopen(filename, "w");

	int Nj = 1, Ni = 1;

	double *a00 = new double [Ni*Nj];

	double Qdifr = 32.*(Nper*Nper*pow(period, 2))/9.;

	Qa = 1./(1./2000. + 1./Qdifr);
	printf("Qa = %g\n", Qa);

	Qa = 900;


//	for(voltage = 15.6; voltage <= 16.; voltage += 0.2)
//	for(int Nper = 11; Nper < 21; Nper+=2)
	for(int j = 0; j < Nj; j++)
    for(int i = 0; i < Ni; i++)
	{
//		kw = 2.*Pi*86.667/299.8*(1. +  (double(i) - 15.)/15.*0.0015);
		double kpd				= sol.solveMultiModeOrotron( current,		 //Ток
													period, Nper,	 //период, число периодов первой секции
													ld, lb,			 //длина дрейфовой секции, длина второй секции
													kw,				 //волновое число
													norm1,			 //норма одного периода первой секции
													norm1B,			 //норма одного периода первой секции
													voltage,	  	 //напряжение
													0.  , delta_freq,//входная мощность, и сдвиг частоты
													Qa, Qb,			 //добротности первой и второй секций
													wall,			 //расстояние  от  стенки и толщина пучка
												    Time, Namm, initAmpMax,
													filename,
													file2name,	  	 //имя файла для печати параметров*//*
													comment);	 //коментарий и начальная точка отсчёта	 
		a00[Ni*j + i] = kpd;


		double Norma= norm1*0.5*double(Nper);
		double G = current/17045.*1./(kw*Norma);
		double gamma = 1.+ voltage/511.;
		double nu = 1./(gamma*(gamma*gamma -1.));
		double beta = sqrt(gamma*gamma -1)/gamma;
		double Delta =  1 - kw/(2.*Pi/period*beta);
		double h = 2*Pi/period;

	//	fprintf(file, "%g\t%g\t%g\n", current, voltage, -/*-Qa/Qdifr*//*voltage*current*kpd);
//		Безразмерные величины (нужно значение a_stat)
//		fprintf(file, "%g\t%g\t%g\t%g\n", sqrt(fabs(astat)*nu*h)*period*double(Nper), (fabs(astat) > 1e-5)? sqrt(h)*Delta/sqrt(fabs(astat)*nu): 0, pow(fabs((astat)), 1.5)*sqrt(h*nu)/(2.*Qa*G), astat);
//		fprintf(file, "%g\t%i\t%g\t%g\t%g\n", voltage, Nper, ld2, lb, output_power);
//		fprintf(file, "%g\t%g\n", 3.*299.8*kw/(2.*Pi), output_power);
	}
//	fclose(file);



	sol.releaseMultiplierMemory();
	cudaDeviceReset();
	
}*/
void multimode_orotronXML(char *inputfilename)
{

	char rootFolder[600], filename[300], file2name[300], *file3name;
	_mkdir("MultiplierTest"); _chdir("MultiplierTest");
	_mkdir("Data1"); _chdir("Data1"); _getcwd(rootFolder, 600);

	
	 QDomDocument doc("inputParameters");
	 QString err_msg; int err_line; int err_column;
	 QFile file("F:\\Piotr\\CalcData\\mm_orotron_Data\\input.xml");
	 if (!file.open(QIODevice::ReadOnly))
		return ;
	 bool set_cont =  doc.setContent(&file, &err_msg, &err_line, &err_column);
	 if (!set_cont) 
	 {
		 cout<<qPrintable(err_msg)<<"; line number"<<err_line<<",  column number"<<err_column<<"\n";
		 file.close();
      return ;
	}
	file.close();

 // print out the element names of all elements that are direct children
 // of the outermost element.

	OrotronMutlimode sol(&doc);
	int Nmax = doc.elementsByTagName("MaxPointNumber").item(0).toElement().text().toInt();
	double Lmax = doc.elementsByTagName("MaxLength").item(0).toElement().text().toDouble();
	double grpSpeedCoeff = doc.elementsByTagName("groupSpeedCoeff").item(0).toElement().text().toDouble();
	QString tmp = doc.elementsByTagName("comment").item(0).toElement().text();
	QByteArray btmp  = tmp.toLocal8Bit();
	char* solcomment; solcomment  = (char*)btmp.data();
	QByteArray btmp1 = ((doc.elementsByTagName("solverName").item(0).toElement().text()).toLocal8Bit());
	char* solverName; solverName =  (char*)btmp1.data();	

	QDomNodeList ndlist1  = doc.elementsByTagName("problemName"); 
	QByteArray btmp2;
	if(ndlist1.length() != 0) 
	{
		btmp2 =  (ndlist1.item(0).toElement().text()).toLocal8Bit();
		file3name	    =  (char*)btmp2.data();
	}
	else
	{
		file3name = new char [100];
		sprintf(file3name, "unnamedProblem");
	}
	


	sol.initSolver(Nmax, Lmax, grpSpeedCoeff, solverName);


//	_mkdir("D:\\Piotr\\Multiplier_260GHz_Data");



// Параметры (ВСЕ)
	double current = doc.elementsByTagName("current").item(0).toElement().text().toDouble();				//Ток
	QDomNode LFsection = doc.elementsByTagName("LFsection").item(0);
	if(LFsection.isNull()){ printf("Error: LFsection entry is not found\n"); return;}
	double period = LFsection.namedItem("period").toElement().text().toDouble();
	int    Nper   = LFsection.namedItem("periodNumbers").toElement().text().toInt();	//период, число периодов первой сек (если вывод излучения адиабатический, то можно немного добавить длины, так как синус немного "выползает"
	double norm1  = LFsection.namedItem("periodNorma").toElement().text().toDouble();	//для прямоугольной гофрировки wall = 0.1: Norm = ; wall = 0.07 Norm = ?
																						//для синусовой гофрировки wall = 0.05: Norm = 33.8; wall = 0.07 Norm = 42;
	double Qa =     LFsection.namedItem("QFactor").toElement().text().toDouble();       //добротности первой (полная) секции
	

	QDomNode HFsection = doc.elementsByTagName("HFsection").item(0);
	if(HFsection.isNull()) {printf("Error: LFsection entry is not found\n"); return;}
	double norm1B=  HFsection.namedItem("periodNorma").toElement().text().toDouble(); 
	double lb  =    HFsection.namedItem("length").toElement().text().toDouble();        //длина длина второй секции
	double  Qb =    HFsection.namedItem("QFactor").toElement().text().toDouble();		// добротности второй (омическая) секции

    double ld =	 doc.elementsByTagName("Dreif").item(0).toElement().text().toDouble(); //длина дрейфовой секции, длина второй секции
	double kw =  2.*Pi/299.8*doc.elementsByTagName("frequency").item(0).toElement().text().toDouble();//волновое число (частота возбуждения)

	
	double p_in = doc.elementsByTagName("power").item(0).toElement().text().toDouble();									//входная мощность, и сдвиг частоты
	double delta_freq =  doc.elementsByTagName("deltaFrequency").item(0).toElement().text().toDouble();
	double voltage =     doc.elementsByTagName("voltage").item(0).toElement().text().toDouble();						//напряжение
	
	double wall = doc.elementsByTagName("beamStructureGap").item(0).toElement().text().toDouble();						//расстояние от стенки и толщина пучка
	int Namm = doc.elementsByTagName("numberOfLongitudinalHarmonics").item(0).toElement().text().toDouble();																											//Количиество учитываемых продольных мод
	double Time = doc.elementsByTagName("Time").item(0).toElement().text().toDouble();									// Время счёта
	double initAmpMax = doc.elementsByTagName("initialAmplitudeMaximum").item(0).toElement().text().toDouble();			//Максимум амплитуды в начальный момент
	
	char comment[] = "Первый прогон";



	

//	devID = findCudaDevice(argc, (const char **)argv);
//  cudaGetDeviceProperties(&deviceProp, devID);


	int file_index = 1;				//Индекс файла с результатами

	QDomElement docElem = doc.documentElement();

	sprintf(filename,  "F:\\Piotr\\CalcData\\mm_orotron_Data\\%s_%i.txt", file3name, file_index);	
/*	FILE *tmpFile;
	while (fopen_s(&tmpFile, filename, "r") == 0)
	{
		fclose(tmpFile);
		sprintf(filename,  "F:\\Piotr\\CalcData\\mm_orotron_Data\\%s_%i.txt", file3name, ++file_index);
	} */


	sprintf(file2name, "F:\\Piotr\\CalcData\\mm_orotron_Data\\%s_params_%i.txt", file3name, file_index);

//	FILE *file =  fopen(filename, "w");
	QDomElement outputFile1 = doc.createElement("outputFile");
	QDomText t = doc.createTextNode(QString(filename));
	outputFile1.appendChild(t);
	
	docElem.appendChild(outputFile1);
	
	

	sprintf(file2name, "F:\\Piotr\\CalcData\\mm_orotron_Data\\%s_params_%i.xml", file3name, file_index);
	FILE *xmlFile = fopen(file2name, "w");
	fprintf(xmlFile, "%s", qPrintable(doc.toString()));
	fclose(xmlFile);
	sprintf(file2name, "F:\\Piotr\\CalcData\\mm_orotron_Data\\%s_params_%i.txt", file3name, file_index);


	int Nj = 1, Ni = 1;

	double *a00 = new double [Ni*Nj];

	double Qdifr = 32.*(Nper*Nper*pow(period, 2))/9.;

	Qa = 1./(1./2000. + 1./Qdifr);
	printf("Qa = %g\n", Qa);

	Qa = 900;


//	for(voltage = 15.6; voltage <= 16.; voltage += 0.2)
//	for(int Nper = 11; Nper < 21; Nper+=2)
	for(int j = 0; j < Nj; j++)
    for(int i = 0; i < Ni; i++)
	{
//		kw = 2.*Pi*86.667/299.8*(1. +  (double(i) - 15.)/15.*0.0015);
		double kpd				= sol.solveMultiModeOrotron( current,		 //Ток
													period, Nper,	 //период, число периодов первой секции
													ld, lb,			 //длина дрейфовой секции, длина второй секции
													kw,				 //волновое число
													norm1,			 //норма одного периода первой секции
													norm1B,			 //норма одного периода первой секции
													voltage,	  	 //напряжение
													0.  , delta_freq,//входная мощность, и сдвиг частоты
													Qa, Qb,			 //добротности первой и второй секций
													wall,			 //расстояние  от  стенки и толщина пучка
												    Time, Namm, initAmpMax,
													filename,        //имя файла для вывода временной зависимости*/
													file2name,	  	 //имя файла для печати параметров*/
													comment);	 //коментарий и начальная точка отсчёта	 
		a00[Ni*j + i] = kpd;


		double Norma= norm1*0.5*double(Nper);
		double G = current/17045.*1./(kw*Norma);
		double gamma = 1.+ voltage/511.;
		double nu = 1./(gamma*(gamma*gamma -1.));
		double beta = sqrt(gamma*gamma -1)/gamma;
		double Delta =  1 - kw/(2.*Pi/period*beta);
		double h = 2*Pi/period;

	//	fprintf(file, "%g\t%g\t%g\n", current, voltage, -/*-Qa/Qdifr*/voltage*current*kpd);
//		Безразмерные величины (нужно значение a_stat)
//		fprintf(file, "%g\t%g\t%g\t%g\n", sqrt(fabs(astat)*nu*h)*period*double(Nper), (fabs(astat) > 1e-5)? sqrt(h)*Delta/sqrt(fabs(astat)*nu): 0, pow(fabs((astat)), 1.5)*sqrt(h*nu)/(2.*Qa*G), astat);
//		fprintf(file, "%g\t%i\t%g\t%g\t%g\n", voltage, Nper, ld2, lb, output_power);
//		fprintf(file, "%g\t%g\n", 3.*299.8*kw/(2.*Pi), output_power);
	}
//	fclose(file);



	sol.releaseDeviceMemory();
	cudaDeviceReset();
	
}
/*void Klinotron1()
{
	double K = 1;

	char rootFolder[600], currentFolder[600], alphaFolder[600], filename[300];
	_mkdir("CudaTest"); _chdir("CudaTest");
	_mkdir("Data5"); _chdir("Data5"); _getcwd(rootFolder, 600);

	
	Wave wave(2*Pi,0,K, 1, 1, 0, 0, 0, 0);
	Beam beam(Beam::VtoU(C*K/wave.geth()), 0.3,  0.0001, 1,     0.3,  1.6*wave.geth());
	         //   U                       //I    //X0   //Xmax   //w     //h l
	beam.alpha = 1./6.;//0.4;
	CData data(6, 100, 0.06*0.5, 12000,   32,  32 , beam);
			 //L          //Nz   //d    //Nt  //Nq //Nx 
	Klinotron sol(data, beam, wave);
//	sol.changeAlphaChangingHeight(0.4);
//	sol.execute();

	
	double L, L0, DL, alpha, alpha0, Dalpha; int NL, Nalpha;
	L0 = 7;	    alpha0 = 0.1;
	DL = 10;	Dalpha = 0.2;
	NL = 7;	  Nalpha = 7;

		
	double *X = new double [NL]; 
	double *Y1 = new double [NL]; double *Y2 = new double [NL]; double *Y3 = new double [NL];


	for(int j = 0; j < Nalpha; j++)
	{
		alpha = alpha0  + Dalpha/(double)Nalpha*(double)j;
		sprintf(currentFolder, "alpa-%g", alpha); _mkdir(currentFolder); _chdir(currentFolder); _getcwd(alphaFolder, 600);

		for(int i = 0; i < NL; i++)
		{
			L = L0 + DL/(double)NL*(double)i;
			sprintf(currentFolder, "Data-%g-%g", alpha, L); _mkdir(currentFolder); _chdir(currentFolder);

			sol.changeL(L);
			sol.changeAlphaChangingHeight(alpha);
			time_t startTime, finishTime; 
			time(&startTime);
			sol.execute();
		//	sol.StartCurrentEvaluation();
			time(&finishTime); printf("\nElapsed time is %g\n",difftime(finishTime, startTime));

			X[i] = L;
			Y1[i] = sol.getKPD();
			Y2[i] = sol.getKPDav();

			_chdir(alphaFolder);
		}
		_chdir(rootFolder);
		sprintf(filename, "kpd_vs_L_a%g.txt", alpha);
		PrintTableX(filename, Y1, X, NL);
		sprintf(filename, "kpdAV_vs_L_a%g.txt", alpha);
		PrintTableX(filename, Y2, X, NL);
	}

}
void MPI_Klinotron(int  argc,  char **argv)
{
	double K = 1;

	char rootFolder[600], currentFolder[600], alphaFolder[600], filename[300];
	_mkdir("ModData3"); _chdir("ModData3");
	_mkdir("Data2"); _chdir("Data2"); _getcwd(rootFolder, 600);

	
	Wave wave(2*Pi,0,K, 1, 1, 0, 0, 0, 0);
	Beam beam(Beam::VtoU(C*K/wave.geth()), 0.3,  0.0001, 1,     0.3,  1.6*wave.geth());
	         //   U                       //I    //X0   //Xmax   //w     //h l
	beam.alpha = 1./6.;//0.4;
	CData data(6, 100, 0.06*0.5, 8000,   33,  33 , beam);
			 //L          //Nz   //dt    //Nt  //Nq //Nx 
	


	
	double L, L0, DL, alpha, alpha0, Dalpha; int NL, Nalpha;
	L0 = 10;     alpha0 = 0.25;
	DL = 30; 	 Dalpha = 0.1;
	NL = 18; 	 Nalpha = 1;

		
	double *X = new double [NL]; 
	double *Y1 = new double [NL]; double *Y2 = new double [NL]; double *Y3 = new double [NL];

	int rank, size;  double buff[3], buff2[20];
/*
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);



	Klinotron sol(data, beam, wave);


	for(int j = 0; j < Nalpha; j++)
	{
		alpha = alpha0  + Dalpha/(double)Nalpha*(double)j;
		sprintf(currentFolder, "alpa-%g", alpha); _mkdir(currentFolder); _chdir(currentFolder); _getcwd(alphaFolder, 600);
		printf("size = %i\n", size);

		for(int i = 0; i < NL/size; i++)
		{
			L = L0 + DL/(double)NL*(double)(i*size+rank);
			sprintf(currentFolder, "Data-%g-%g", alpha, L); _mkdir(currentFolder); _chdir(currentFolder);

			sol.changeL(L);
			sol.changeAlphaChangingHeight(alpha);
			sol.changeNx((int) (L*alpha)*5);
			sol.execute();

			buff[0] = L;
			buff[1] = sol.getKPD();
			buff[2] = sol.getKPDav();

			_chdir(alphaFolder);
			MPI_Gather(buff, 3, MPI_DOUBLE, buff2, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

			if(rank == 0)
			{
			
				for(int r = 0; r < size; r++) {X[(size)*i+r] = buff2[3*r+0]; Y1[(size)*i+r] = buff2[3*r+1]; Y2[(size)*i+r] = buff2[3*r+2];}
			}

		}

		if(rank == 0)
		{	_chdir(rootFolder);
			sprintf(filename, "kpd_vs_L_a%g.txt", alpha);
			PrintTableX(filename, Y1, X, NL);
			sprintf(filename, "kpdAV_vs_L_a%g.txt", alpha);
			PrintTableX(filename, Y2, X, NL);
		}
	}
	printf("\n rank %i finalized\n", rank);
	fflush(0);
	MPI_Finalize();
}*/
int main(int argc, char *argv[])
{
	int devNum = -1;
	cudaGetDeviceCount(&devNum);
	printf("Input parameter: %s\n", argv[1]);
	printf("There are %i devices found\n", devNum);


	cudaDeviceProp devs[3], the_dev;
	for(int q = 0; q < devNum; q++)
	{
		cudaGetDeviceProperties(devs + q, q);
		printf("Device %i is %s has capability %i.%i\n", q, devs[q].name, devs[q].major, devs[q].minor);
	}

	sprintf(the_dev.name, "Tesla C2075");
	the_dev.major = 2;

	cudaChooseDevice(&devNum, &the_dev);
	printf("DEVICE %s IS SELECTED\n",  devs[devNum].name);
	cudaSetDevice(devNum);

	cudaDeviceReset();

	
//	multiplier();
//	orotron();
//	twt();
//	bwo();
//	multimode_orotron();
//	multimode_orotronXML("dummy");
	startXML(argc, argv);
	cudaDeviceReset();
	fflush(0);
	return 0;
}
