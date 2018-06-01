#include "multiplier.h"
#include "xml_routines.h"



Multiplier::Multiplier(QDomDocument *doc) : Device(doc)
{
	B_stat = 0.;
	QDomNode HFsection = doc->elementsByTagName("HFsection").item(0);
	if (HFsection.isNull()) {
		printf("Warning: HFsection entry is not found\n");
	}
	else {
		double norm1B;
		setXMLEntry(&HFsection, "periodNorma", &norm1B);
		setXMLEntry(&HFsection, "length", &Lb);  //длина длина второй секции
		setXMLEntry(&HFsection, "QFactor", &Qb); // добротности второй (омическая) секции

	}
	setXMLEntry(doc, "Dreif", &Ld);									//длина дрейфовой секции, длина второй секции

}

void lsolve(double *X, double *Y, double *dX, int N) // Y dX = X
{
	double tmp;

	for (int s = 0; s < N; s++)
	{
		tmp = Y[N*s + s];

		for (int l = s; l < N; l++) Y[N*s + l] /= tmp;
		X[s] /= tmp;

		for (int l = s + 1; l < N; l++)
		{
			tmp = Y[N*l + s];
			for (int k = s; k < N; k++)
			{
				Y[N*l + k] -= Y[N*s + k] * tmp;
			}
			X[l] -= tmp*X[s];
		}


	}

	for (int s = N - 1; s > -1; s--)
	{
		dX[s] = X[s];
		for (int l = s + 1; l < N; l++)
		{
			dX[s] -= Y[N*s + l] * dX[l];

		}
		dX[s] /= Y[N*s + s];
	}
}

double Multiplier::getHFoutputPower(double _Current_ampers, double _period, int _Nperiods, double _Ld, double _Lb, double _k1, double _Norma, double Norma1B,
	double _voltage, double _inputPower_watts, double _delta, double _Qa, double _Qb, double _wall, double _sQ1, double _sQ3,
	char *filename, char *comment, double *monitor, char * difrFlag)
{

	Current_ampers = _Current_ampers;   Ld = _Ld; Lb = _Lb;
	period = _period;				    Norma = _Norma;
	Nperiods = _Nperiods;			  	voltage = _voltage;
	delta_freq = _delta;				Q = _Qa; Qb = _Qb;
	inputPower_watts = _inputPower_watts;
	spchQ1 = _sQ1;
	spchQ3 = _sQ3;

	Q_difr = diffractionQ(difrFlag);

	//	delta_freq = 0;

	k1 = k0 + delta_freq*2.*PI / 299.8;


	//	loadField("F:\\Piotr\\Multiplier_260GHz_Data\\fieldStructure_hot.txt");

	if (!fieldLoaded) Norma = _Norma*0.5*double(Nperiods);
	Norma1b = Norma1B;
	wall = _wall;

	if (A_stat == 0) PrintParams(filename, comment);

	int N_it = 50;

	double G = paramG(); //k1 в мм^{-1}, Norma в мм^3.
	double omega = k1*(10.*c);
	double h = 2.*Pi / period, d = period;
	g1 = sqrt(h*h - k1*k1); g3 = 3 * g1;

	//	if(!solverStatus) 	solverStatus = initMultiplierSolver(3000, 10.1, 0., Ns);

	double nextA = 0.000019;
	double nextB = 0;
	// Без разброса
	double hot_displacement = 0.25; 		//Nper = 19 => hot_displavement = 0.2135; voltage += 0.081
	double deltaVolt = 0.095;				//Nper = 13 => hot_displavement = 0.1; voltage += 0.038	
	//Nper = 23 => hot_displacement = 0.293; voltage += 0.112
	// С разбросом 
	//Nper = 23 => hot_displacement = 0.3; voltage +=  .114
	//	   = 21 =>					= 0.25		   +=  .095	
	if (s_stat == 0)
	{
		double tmp = k1;
		k1 = k0 + hot_displacement*2.*PI / 299.8;
		voltage += deltaVolt;
		s_stat = findAstatDetuned_full(nextA, N_it, G);
		voltage -= deltaVolt;
		k1 = tmp;
	}

	if (A_stat != 0) { nextA = A_stat; N_it = 8; } //A_stat не зависит от длин секций расположенных правее, но зависит от напряжения и числа периодов. 


	if (delta_freq == 0)
	{

		A_stat = findAstat(nextA, N_it, G);
	}
	else
	{
		voltage += deltaVolt;
		A_stat = findAstatDetuned(nextA, N_it, G, s_stat);
		voltage -= deltaVolt;
	}


	//	nextA = 0.000947;

	cplx b_k[5], tm_b;
	double coeffk[4] = { 0, 0.5, 0.5, 1. };
	b[0] = 0.;

	k0 += hot_displacement*2.*PI / 299.8;

	double NormaB = Norma1B*Lb / (2.*0.28);
	double Gb = Current_ampers / 17045.*1. / (3.*k1*NormaB);
	int Nb_it;
	double deltaT = 200;


	printf("\n");
	printf("Gb = %g \n", Gb);



	cplx B00;
	if (abs(B_stat) < 1e-50) { B00 = 1.57e-6*(1. + 0.001*I); Nb_it = 50; }
	else { B00 = B_stat; Nb_it = 50; }


	if (delta_freq == 0)
		B_stat = FindBstat(B00, Nb_it, Gb, A_stat);
	else
	{
		voltage += deltaVolt;
		B_stat = FindBstatDetuned(B00, Nb_it, Gb, A_stat);
		voltage -= deltaVolt;
	}


	k0 -= hot_displacement*2.*PI / 299.8;

	printf("B_{stat}  = %g + %g i\n", real(B_stat), imag(B_stat));

	double energy;
	cudaMemcpy((void*)&energy, d_avEN, sizeof(double), cudaMemcpyDeviceToHost);
	energy /= double(Np*Nq*Ns*Nv);
	double Qfull = Q_difr*Qb / (Q_difr + Qb);
	printf("Energy balance = %g\n", abs(B_stat)*abs(B_stat) / (2.*Qfull) + Gb*energy);

	double outputPower = retriveBPower(B_stat, omega, Qb, NormaB);
	printf("HF output power  = %g\t Residual = %g\n", outputPower, abs(b[Nb_it - 1]) - abs(b[Nb_it - 2]));


	*monitor = A_stat;//imag(B_stat);//A_stat;//abs(B_stat)*abs(B_stat)/(2.*Qfull*Gb);//monitor_par;//;

	return outputPower;
}

double Multiplier::retriveBPower(cplx B, double omega, double Qb, double NormaB)
{
	if (Q_difr < 300) Q_difr = 300;
	double Qfull = Q_difr*Qb / (Q_difr + Qb);
	double electronPower = pow(m*c*c / e, 2)*(3.*omega) / Qfull*1e-7*pow(abs(B)*10., 2)*NormaB*1e-3*0.5;
	double outputPower = electronPower*Qb / (Qb + Q_difr);
	return outputPower;
}

double Multiplier::diffractionQ(char *difrQflag)
{
	sprintf(difractionFlag, "%s", difrQflag);
	if (strcmp(difrQflag, "Ohm") == 0)  return Qb;
	if (strcmp(difrQflag, "Min") == 0)  return 18.7*(Lb*Lb);
	if (strcmp(difrQflag, "Empiric") == 0)  return 314.286 + 10.23834*Lb*Lb;

	sprintf(difractionFlag, "%s", "Ohm");

	return Qb;


}
cplx Multiplier::FindBstat(cplx b0, int Nb_it, double Gb, double Astat)
{

	double xi = 0.001;
	cplx F00, F01, F10, B00, B01, B10;
	double dX, dY, J11, J12, J21, J22, det;
	//	B00 = 1.57e-5*(1.+0.01*I); Nb_it = 15;

	/*	double temp_re[3000];
	double temp_im[3000];*/

	B00 = b0;


	double Qfull = Q_difr*Qb / (Q_difr + Qb);
	/*
	cplx cB = HfElectronCurrent(real(B00), imag(B00), Astat);
	printf("real(J3) = %g\timag(J3) = %g\n", real(cB), imag(cB));*/



	for (int i = 1; i < Nb_it; i++)
	{
		/*	tm_b = 4e-5*I;
		cplx tmp_B = CurrentB(real(tm_b), imag(tm_b), double(Nperiods)*d, Lb, Ld, h, k1, nextA, d, Nperiods,15.8);
		printf("\n host:\t%g\t%g\n", real(tmp_B), imag(tmp_B));

		cudaMemcpy(temp_re, d_rJ3, sizeof(double)*3000, cudaMemcpyDeviceToHost);

		FILE *file = fopen("d:\\Piotr\\debug_re.txt", "w");
		for(int pr = 0; pr < 3000; pr++) fprintf(file, "%g\t%g\n",(double)pr*(17.1)/3000., temp_re[pr]/double(NP*NQ*NS));
		fclose(file);

		cudaMemcpy(temp_im, d_iJ3, sizeof(double)*3000, cudaMemcpyDeviceToHost);

		file = fopen("d:\\Piotr\\debug_im.txt", "w");
		for(int pr = 0; pr < 3000; pr++) fprintf(file, "%g\t%g\n",(double)pr*(17.1)/3000., temp_im[pr]/double(NP*NQ*NS));
		fclose(file);*/

		/*
		cplx deb_tmp;

		for(int k = 1; k <= 4; k++)
		{
		tm_b = coeffk[k-1]*b_k[k-1] + b[i-1];
		deb_tmp = CurrentB(real(tm_b), imag(tm_b), double(Nperiods)*d, Lb, Ld, h, k1, nextA, d, Nperiods,voltage);
		b_k[k] = deltaT*(-tm_b/(2.*Qb) + Gb*deb_tmp);//CurrentB(real(tm_b), imag(tm_b), double(Nperiods)*d, Lb, Ld, h, k1, nextA, d, Nperiods,15.8));
		}
		b[i] = b[i-1] + (b_k[1] +2.*b_k[2] + 2.*b_k[3] + b_k[4])/6.;
		funcF[i] = abs(deb_tmp);*/


		printf("B00:\t%g\t%g\n", real(B00), imag(B00));

		F00 = -B00 / (2.*Qfull) + Gb*HfElectronCurrent(real(B00), imag(B00), Astat);

		B01 = B00 + xi*abs(B00);
		F01 = -B01 / (2.*Qfull) + Gb*HfElectronCurrent(real(B01), imag(B01), Astat);

		B10 = B00 + I*xi*abs(B00);
		F10 = -B10 / (2.*Qfull) + Gb*HfElectronCurrent(real(B10), imag(B10), Astat);


		J11 = real(F01 - F00) / (xi*abs(B00));  J12 = real(F10 - F00) / (xi*abs(B00));
		J21 = imag(F01 - F00) / (xi*abs(B00));  J22 = imag(F10 - F00) / (xi*abs(B00));

		det = J11*J22 - J12*J21;

		dX = -(J22*real(F00) - J12*imag(F00)) / det;
		dY = -(-J21*real(F00) + J11*imag(F00)) / det;

		cplx temp = real(B00) + dX + I*(imag(B00) + dY);

		if (abs(temp) < 1e-20) break;
		B00 = temp;
		printf("|B|  = %g\t%g\t%g\tdet = %g\n", abs(B00), real(B00), imag(B00), det);
	}
	return B00;
}
cplx Multiplier::FindBstatDetuned(cplx b0, int Nb_it, double Gb, double Astat)
{

	double xi = 0.001;
	cplx F00, F01, F10, B00, B01, B10;
	double dX, dY, J11, J12, J21, J22, det;
	//	B00 = 1.57e-5*(1.+0.01*I); Nb_it = 15;

	B00 = b0;
	double detuning = 0.5*(k1*k1 - k0*k0) / k1*k1;


	double Qfull = Q_difr*Qb / (Q_difr + Qb);


	for (int i = 1; i < Nb_it; i++)
	{

		printf("B00:\t%g\t%g\n", real(B00), imag(B00));

		F00 = -B00*(1. / (2.*Qfull) + I*detuning) + Gb*HfElectronCurrent(real(B00), imag(B00), Astat);

		B01 = B00 + xi*abs(B00);
		F01 = -B01*(1. / (2.*Qfull) + I*detuning) + Gb*HfElectronCurrent(real(B01), imag(B01), Astat);

		B10 = B00 + I*xi*abs(B00);
		F10 = -B10*(1. / (2.*Qfull) + I*detuning) + Gb*HfElectronCurrent(real(B10), imag(B10), Astat);


		J11 = real(F01 - F00) / (xi*abs(B00));  J12 = real(F10 - F00) / (xi*abs(B00));
		J21 = imag(F01 - F00) / (xi*abs(B00));  J22 = imag(F10 - F00) / (xi*abs(B00));

		det = J11*J22 - J12*J21;

		dX = -(J22*real(F00) - J12*imag(F00)) / det;
		dY = -(-J21*real(F00) + J11*imag(F00)) / det;

		if (abs(dX + I*dY) / abs(B00) < 1e-12) break;

		cplx temp = real(B00) + dX + I*(imag(B00) + dY);

		B00 = temp;

		printf("|B|  = %g\t%g\t%g\tdet = %g\n", abs(B00), real(B00), imag(B00), det);
	}
	return B00;
}
void Multiplier::PrintParams(char *filename, char *comment)
{
	FILE *file2 = fopen(filename, "wa");

	fprintf(file2, "\"Current\"\t%g\n\"period\"\t%g\n\"Nper\"\t%i\n", Current_ampers, period, Nperiods);
	fprintf(file2, "\"ld\"\t%g\n\"lb\"\t%g\n\"kw\"\t%g\n", Ld, Lb, k1);
	fprintf(file2, "\"norm1\"\t%g\n\"norm2\"\t%g\n\"p_in\"\t%g\n\"delta_freq\"\t%g\n", Norma, Norma1b, inputPower_watts, delta_freq);
	fprintf(file2, "\"voltage\"\t%g\n\"Qa\"\t%g\n\"Qb\"\t%g\n", voltage, Q, Qb);
	fprintf(file2, "\"beam thickness\"\t%g\n\"LF space charge amp\"\t%g\n\"HF space charge amp\"\t%g\n", wall, spchQ1, spchQ3);
	fprintf(file2, "\"particle number phase\"\t%i\n\"particle number radius\"\t%i\n\"particle number cyclotron phase\"\t%i\n\"particle number energy spread\"\t%i\n", Np, Nq, Ns, Nv);
	fprintf(file2, "\"%s\"\n", solverName);
	fprintf(file2, "\"%с\"\n", difractionFlag);
	fprintf(file2, "\"Np\"\t%i\n\"Nq\"\t%i\n\"Ns\"\t%i\n\"Nv\"\t%i\n", Np, Nq, Ns, Nv);
	if (strcmp(difractionFlag, "Empiric") == 0) fprintf(file2, "\"314.286 + 10.23834*Lb*Lb\"\n");
	fprintf(file2, "\"%s\"\n", comment);
	fclose(file2);
}
double Multiplier::findAstatDetuned(double nextA, int N_it, double G, double Sin)
{
	cplx A00, A01, A10;
	cplx F00, F01, F10;

	double dX, dY, J11, J12, J21, J22, det;

	double xi = 1e-5; // 1e-6

	double d = period;

	A00 = nextA;

	double detuning = 0.5*(k1*k1 - k0*k0) / (k1*k1);
	printf("findAstatDetuned: detuning = %g\n", delta_freq);

	double omega = k1*(10.*c);
	/////////
	double inputPower_dimensionless = inputPower_watts*1.E7 / (omega / 1.*(m*c*c / e)*(m*c*c / e)*(Current_ampers / 17045.)*(1. / (k1 * 10))); // = S a
	////////
	double inputCurrent = (Sin == 0) ? sqrt(inputPower_dimensionless / (2.*Q*G)) : Sin;

	/*	printf("Qa = %g,\t G = %g,\t inPower = %g\tinCurrent = %g \n", Q, G, inputPower_dimensionless, inputCurrent);
	double linA = 1e-7;
	printf("(A0*A0 + 2.*Q*G*DeltaEnergy(A0)  = %g\n", abs(linA*linA + 2.*Q*G*DeltaEnergy(linA))); // Tecт на стартовый ток F > 0 -> старта нет.
	cplx cA = CurrentA(real(A00), imag(A00));
	cplx cA2 = ElectronCurrentA(real(A00), imag(A00));
	printf("cA = %g,  %g\ncA2 = %g,  %g\n", real(cA), imag(cA), real(cA2), imag(cA2));

	printf("findAstatDetuned: %g = %g\n", real(cA*conj(A00)), DeltaEnergy(abs(A00)));*/

	//	N_it = 0;
	for (int i = 0; i < 10 + 0 * N_it; i++)
	{

		F00 = A00*(1. / (2.*Q) + I*detuning) + G*(-inputCurrent - ElectronCurrentA(real(A00), imag(A00)));

		A01 = A00 + xi*abs(A00);
		F01 = A01*(1. / (2.*Q) + I*detuning) + G*(-inputCurrent - ElectronCurrentA(real(A01), imag(A01)));

		A10 = A00 + I*xi*abs(A00);
		F10 = A10*(1. / (2.*Q) + I*detuning) + G*(-inputCurrent - ElectronCurrentA(real(A10), imag(A10)));

		J11 = real(F01 - F00) / (xi*abs(A00));  J12 = real(F10 - F00) / (xi*abs(A00));
		J21 = imag(F01 - F00) / (xi*abs(A00));  J22 = imag(F10 - F00) / (xi*abs(A00));

		det = J11*J22 - J12*J21;

		dX = -(J22*real(F00) - J12*imag(F00)) / det;
		dY = -(-J21*real(F00) + J11*imag(F00)) / det;

		if (abs(dX + I*dY) / abs(A00) < 1e-15) break;
		if (abs(F00) < 1e-30) break;

		cplx temp = real(A00) + dX + I*(imag(A00) + dY);

		printf("%g\t%g\t%g\t%g\t%g\n", abs(A00*A00) + 2.*Q*G*(-inputPower_dimensionless /*+ ElectronsDeltaEnergy(abs(A00))*/),
			abs(A00),
			G*inputCurrent,
			ElectronsDeltaEnergy(abs(A00)),
			abs(A00*A00));
		//		printf("A_{st}= %g\tdet = %g\t%g\t%g\n", abs(A00), det, abs(A00 - temp), abs(CurrentA(real(A00), imag(A00))));

		A00 = temp;



	}
	printf("A_{stat}  = %g\t Residual = %g\n", abs(A00), abs(A10 - A00));

	monitor_par = //pow(abs(A00),2) + 2.*Q*G*(-inputPower_dimensionless + ElectronsDeltaEnergy(abs(A00)));
		real(A00)*inputCurrent / inputPower_dimensionless;

	return abs(A00);

}

double Multiplier::findAstatDetuned_full(double nextA, int N_it, double G)
{
	cplx A00, A01, A10;
	cplx F00, F01, F10;

	double xi = 0.0001;
	double omega = k1*(10.*c);
	double d = period;

	A00 = nextA;

	double detuning = 0.5*(k1*k1 - k0*k0) / (k1*k1);
	/////////
	double inputPower_dimensionless = inputPower_watts*1.E7 / (omega / 1.*(m*c*c / e)*(m*c*c / e)*(Current_ampers / 17045.)*(1. / (k1 * 10))); // = S a
	double inputCurrent = sqrt(inputPower_dimensionless / (2.*Q*G));

	double A[3], J[9];

	double X[3], deltaX[3];
	X[0] = nextA*cos(0.05);
	X[1] = nextA*sin(0.05);
	X[2] = inputCurrent;
	cplx curr[3];
	////////
	for (int i = 0; i < 20 + 0 * N_it; i++)
	{

		curr[0] = ElectronCurrentA(X[0], X[1]);
		A[0] = -0.5*X[0] / (Q)+detuning*X[1] + G*real(curr[0]) + G*X[2];
		A[1] = -0.5*X[1] / (Q)-detuning*X[0] + G*imag(curr[0]);
		A[2] = X[0] * X[2] - inputPower_dimensionless;

		for (int a = 0; a < 2; a++)
		{
			double tmp = X[a];
			X[a] *= (1. + xi);
			curr[a + 1] = ElectronCurrentA(X[0], X[1]);
			X[a] = tmp;
		}

		J[3 * 0 + 0] = -0.5 / Q + G*real(curr[1] - curr[0]) / (X[0] * xi);  //D A'/da'
		J[3 * 0 + 1] = detuning + G*real(curr[2] - curr[0]) / (X[1] * xi);  //D A'/da''
		J[3 * 0 + 2] = G;												  //DA'/ds

		J[3 * 1 + 0] = -detuning + G*imag(curr[1] - curr[0]) / (X[0] * xi);  //D A''/da'
		J[3 * 1 + 1] = -0.5 / Q + G*imag(curr[2] - curr[0]) / (X[1] * xi);  //D A''/da''
		J[3 * 1 + 2] = 0;												  //DA''/ds

		J[3 * 2 + 0] = X[2];	 //D P''/da'
		J[3 * 2 + 1] = 0;		 //D P''/da''
		J[3 * 2 + 2] = X[0];	 //D P''/ds

		lsolve(A, J, deltaX, 3);

		//	printf("dX:%g, %g, %g, %g\n",  deltaX[0], deltaX[1], deltaX[2]);

		double absdX = 0; for (int e = 0; e < 3; e++) absdX += deltaX[e] * deltaX[e]; absdX = sqrt(absdX);
		double absF = 0; for (int e = 0; e < 3; e++) absF += A[e] * A[e]; absF = sqrt(absF);

		//	printf("X:%g, %g, %g\t%g\t%g\n", X[0], X[1], X[2], sqrt(X[0]*X[0] + X[1]*X[1]), X[0]*X[2]);
		printf("%g\t%g\t%g\t%g\t%g\n", inputCurrent,
			X[2],
			X[0] * X[2],
			abs(X[0] + I*X[1]),
			sqrt(2.*Q*G*inputPower_dimensionless));

		if (absF < 1e-30) break;
		if (absdX < 1e-20) break;

		for (int e = 0; e < 3; e++) X[e] -= deltaX[e];

	}

	return  X[2];//abs(X[0] + I*X[1]);
}

double Multiplier::ElectronsDeltaEnergy(double A)
{
	double res = 0;
	if (strcmp(solverName, "multiplier") == 0)  return DeltaEnergy(A);

	return res;

}
cplx Multiplier::ElectronCurrentA(double reA, double imA)
{
	cplx res = 0;
	if (strcmp(solverName, "multiplier") == 0)  return CurrentA(reA, imA);
	return res;
}
cplx Multiplier::HfElectronCurrent(double rB, double iB, double Astat)
{
	cplx res = 0;
	if (strcmp(solverName, "multiplier") == 0) return CurrentB(rB, iB, Astat);
	return res;

}
