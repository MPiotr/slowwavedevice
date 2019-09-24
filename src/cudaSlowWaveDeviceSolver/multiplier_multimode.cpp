 #define lapack_complex_float std::complex<float>
 #define lapack_complex_double std::complex<double>

#include "device.h"
#include "multiplier_multimode.h"
#include "xml_routines.h"
//#include <mkl.h>
#include "cu_mult.h"


void lsolve(double *X, double *Y, double *dX) // Y dX = X
{
	double tmp;
	int N = 4;

		for(int s = 0; s < N; s++)
		{
			tmp = Y[N*s + s];

			for(int l = s; l < N; l++) Y[N*s + l] /= tmp; 
			X[s] /= tmp;

			for(int l = s+1; l < N; l++)
			{
				tmp = Y[N*l + s];
				for(int k = s; k < N; k++)
				{
					Y[N*l + k] -= Y[N*s + k]*tmp;
				}
				X[l] -= tmp*X[s];
			}


		}

		for(int s = N-1; s > -1; s--)
		{
			dX[s] = X[s];
			for(int l = s+1; l < N; l++)
			{
				dX[s] -= Y[N*s + l]*dX[l];
				
			}
			dX[s] /= Y[N*s+s];
		}



}
MultiplierMultiModes::MultiplierMultiModes():Multiplier()
{
	Namm = 1;					//Количиество учитываемых продольных мод (multimode oro)
}

MultiplierMultiModes::MultiplierMultiModes(QDomDocument * doc) :Multiplier(doc)
{
	setXMLEntry(doc, "beamStructureGap", &Namm);					//Количиество учитываемых продольных мод (multimode oro)
}

double MultiplierMultiModes::getHFoutputPowerMultiModes(double _Current_ampers, double _period, int _Nperiods, double _Ld, double _Lb, double _k1, double _Norma, double Norma2, double _voltage, double _inputPower_watts, double _delta, double _Qa1, double _Qa2, double _Qb, double _wall, char *filename, char *comment)
{

	Current_ampers = _Current_ampers;   Ld1= _Ld; Lb = _Lb;
	period	= _period;				  	k1 = _k1; 
	Nperiods= _Nperiods;			  	voltage = _voltage;
	delta_freq = _delta;		    	Qa1 = _Qa1; Qa2 = _Qa2;  Qb = _Qb; 
	inputPower_watts = _inputPower_watts; wall = _wall;

	
	Norma1 = _Norma*0.5*double(Nperiods);

	Norma1b = Norma2;

	if(A_stat == 0) PrintParamsMultiModes(filename, comment);


	int N_it = 150;
	double G1 = Current_ampers/17045.*1./(k1*Norma1); //k1 в мм^{-1}, Norma в мм^3.
	double G2 = Current_ampers/17045.*1./(k1*Norma2); //k1 в мм^{-1}, Norma в мм^3.

	double omega = k1*(10.*c);
	double h = 2.*Pi/period, d = period;
	g1 = sqrt(h*h - k1*k1); g3 = 3*g1;

	cplx A2_stat;
	

	double  nextA = 0.000019;
	double nextX[5] = {2.7e-5, 1e-6, 1e-5, -4e-5, 1e-6};
	
	double nextB = 0;
//	if(A_stat != 0) {nextA = A_stat; N_it = 1;}

	k1 = k0;

	if(k0 == k1)
		findAstatMultiModes(nextX, N_it, G1, G2);
	else
	{	printf("\n Error: detuned mode is not yet implemented\n"); return 0;}

	double NormaB = Norma1b*Lb/(2.*0.28);
	double Gb = Current_ampers/17045.*1./(3.*k1*NormaB);
	int Nb_it = 1;// 7;

	printf("\n");
	printf("Lb = %g \n", Lb);
	printf("B_stat:\t%g\t%g\n", real(B_stat), imag(B_stat));

/*
	cplx B00;
	if(abs(B_stat) < 1e-50) {B00 = 1.57e-6*(1.+0.01*I); Nb_it = 15;} else {B00 = B_stat; Nb_it = 8;}

	if(k0 == k1)
		B_stat = FindBstatDoubleScheme(B00, Nb_it, Gb, A_stat, A2_stat);
	else
	{printf("\n Error: detuned mode is not yet written\n"); return 0;}



	printf("B_{stat}  = %g + %g i\n", real(B_stat), imag(B_stat));
	
	double energy;
	cudaMemcpy((void*) &energy, d_avEN, sizeof(double), cudaMemcpyDeviceToHost);

	energy /= double(Np*Nq*Ns*Nv);

	printf("Energy balance = %g\n", abs(B_stat)*abs(B_stat)/(2.*Qb) + Gb*energy);

	double outputPower = retriveBPower(B_stat, omega, Qb, NormaB);

	printf("HF output power  = %g\n", outputPower);*/

	return 0.;//outputPower;

}
void MultiplierMultiModes::findAstatMultiModes(double* nextX, int N_it, double G1, double G2)
{
	double xi = 0.001;
	double omega = k1*(10.*c);
	/////////
	double inputPower_dimensionless = inputPower_watts*1.E7 / (omega / 1.*(m*c*c / e)*(m*c*c / e)*(Current_ampers / 17045.)*(1. / (k1*10.))); // = S a
	////////

	printf("Qa = %g,\t G = %g,\t inputPowerDimensionless = %g \n", Q, G1, inputPower_dimensionless);

	//FILE *debfile = fopen("d:\\Piotr\\w_orotron_Data\\debug.txt", "w");

	double A[5];
	double J[25];

	cplx current[6];
	double X[5]; memcpy(X, nextX, 5 * sizeof(double));
	double dX[5];
	X[4] = inputPower_dimensionless / X[0];
	cplx J1 = 0, J2 = 0;
	cplx cur1[5]; cplx cur2[5];


	for (int i = 0; i < N_it; i++)
	{
		/*	A0 = double(i)/double(N_it)*0.005;
		fprintf(debfile, "%g\t%g\n", A0, DeltaEnergy(A0));*/

		A[0] = -0.5*X[0] / (Qa1)+G1*real(current[0]) + G1*inputPower_dimensionless / (X[0]);
		A[1] = -0.5*X[1] / (Qa1)+G1*imag(current[0]);
		A[2] = -0.5*X[2] / (Qa2)+G2*real(current[0]);
		A[3] = -0.5*X[3] / (Qa2)+G2*imag(current[0]);
		//		A[4] = X[0]*X[4] - inputPower_dimensionless;


		//		CurrentAMultiModes(X, &J1, &J2); поправить. Испортилось после добавления мультимод
		cur1[0] = J1; cur2[0] = J2;
		printf("J:%g, %g\t", abs(J1), abs(J2));
		for (int a = 0; a < 4; a++)
		{
			double tmp = X[a];
			X[a] *= (1. + xi);
			//			CurrentAMultiModes(X, &J1, &J2); поправить. Испортилось после добавления мультимод
			cur1[a + 1] = J1;
			cur2[a + 1] = J2;
			X[a] = tmp;
		}
		for (int a = 0; a < 4; a++)
		{
			J[4 * 0 + a] = -((a == 0) ? 0.5 / Qa1 : 0) + G1*((a<4) ? 1. / (X[a] * xi) : 0)*real(cur1[a + 1] - cur1[0]) - ((a == 0) ? 1 : 0)*G1*inputPower_dimensionless / (X[0] * X[0]);		  //a1'
			J[4 * 1 + a] = -((a == 1) ? 0.5 / Qa1 : 0) + G1*((a<4) ? 1. / (X[a] * xi) : 0)*imag(cur1[a + 1] - cur1[0]);																 //a1''
			J[4 * 2 + a] = -((a == 2) ? 0.5 / Qa2 : 0) + G2*((a<4) ? 1. / (X[a] * xi) : 0)*real(cur2[a + 1] - cur2[0]);																//a2'
			J[4 * 3 + a] = -((a == 3) ? 0.5 / Qa2 : 0) + G2*((a<4) ? 1. / (X[a] * xi) : 0)*imag(cur2[a + 1] - cur2[0]);															   //a2''
			//		J[5*4 + a] = ((a==4)? X[0]: 0) + ((a==0)? X[4]: 0);											  //s
		}

		lsolve(A, J, dX);
		//	printf("dX:%g, %g, %g, %g\n",  dX[0], dX[1], dX[2], dX[3]);

		double absdX = 0; for (int e = 0; e < 4; e++) absdX += dX[e] * dX[e]; absdX = sqrt(absdX);

		printf("X:%g, %g, %g, %g \n", X[0], X[1], X[2], X[3]);
		if (absdX < 1e-20) break;

		for (int e = 0; e < 4; e++) X[e] -= dX[e];

	}

	memcpy(nextX, X, 4 * sizeof(double));



	printf(" 2 q G = %g \n", 2.*Qa1*G1);
	//fclose(debfile);



}

void MultiplierMultiModes::PrintParamsMultiModes(char *filename, char *comment)
{
	FILE *file2 = fopen(filename, "w");


	fprintf(file2, "\"Current\"\t%g\n\"period\"\t%g\n\"Nper\"\t%i\n", Current_ampers, period, Nperiods);
	fprintf(file2, "\"ld\"\t%g\n\"kw\"\t%g\n", Ld, k1);
	fprintf(file2, "\"lb\"\t%g\n", Lb);
	fprintf(file2, "\"norm1\"\t%g\n\"norm2\"\t%g\n\"p_in\"\t%g\n\"delta_freq\"\t%g\n", Norma, Norma1b, inputPower_watts, delta_freq);
	fprintf(file2, "\"voltage\"\t%g\n\"Qa1\"\t%g\n\"Qa2\"\t%g\n\"Qb\"\t%g\n", voltage, Qa1, Qa2, Qb);
	fprintf(file2, "\"beam thickness\"\t%g\n", wall);
	fprintf(file2, "\"%s\"\n", comment);

	fclose(file2);
}


