#include "multiplier_three_cavity.h"
#include "cu_mult.h"


double MultiplierTreeCavity::getHFoutputPowerDoubleScheme(double _Current_ampers, double _period, int _Nperiods, double _Ld1, double _Ld2, double _La2, double _Lb, double _k1, double _Norma, double Norma1B, double _voltage, double _inputPower_watts, double _delta, double _Qa1, double _Qa2, double _Qb, double _wall, char *filename, char *comment)
{

	Current_ampers = _Current_ampers;   Ld1 = _Ld1; Ld2 = _Ld2; La2 = _La2; Lb = _Lb;
	period	= _period;				  	k1 = _k1; 
	Nperiods= _Nperiods;			  	voltage = _voltage;
	delta_freq = _delta;		    	Qa1 = _Qa1; Qa2 = _Qa2;  Qb = _Qb; 
	inputPower_watts = _inputPower_watts; wall = _wall;

	
	Norma1 = _Norma*0.5*double(Nperiods);
	Norma2 = _Norma*0.5*La2/period;
	Norma1b = Norma1B;

	if(A_stat == 0) PrintParamsDoubleScheme(filename, comment);


	int N_it = 30;
	double G1 = Current_ampers/17045.*1./(k1*Norma1); //k1 β μμ^{-1}, Norma β μμ^3.
	double G2 = Current_ampers/17045.*1./(k1*Norma2); //k1 β μμ^{-1}, Norma β μμ^3.

	double omega = k1*(10.*c);
	double h = 2.*Pi/period, d = period;
	g1 = sqrt(h*h - k1*k1); g3 = 3*g1;

	cplx A2_stat;
	

	double  nextA = 0.000019;
	double nextB = 0;
	if(A_stat != 0) {nextA = A_stat; N_it = 1;}

	k1 = k0;

	if(k0 == k1)
		A2_stat = findAstatDoubleScheme(nextA, N_it, G1, G2);
	else
	{	printf("\n Error: detuned mode is not yet written\n"); return 0;}

	double NormaB = Norma1b*Lb/(2.*0.28);
	double Gb = Current_ampers/17045.*1./(3.*k1*NormaB);
	int Nb_it = 1;// 7;

	printf("\n");
	printf("Lb = %g \n", Lb);
	printf("B_stat:\t%g\t%g\n", real(B_stat), imag(B_stat));


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

	printf("HF output power  = %g\n", outputPower);

	return outputPower;

}

cplx MultiplierTreeCavity::findAstatDoubleScheme(double nextA, int N_it, double G1, double G2)
{
	double omega = k1*(10.*c);
																/////////
	double inputPower_dimensionless = inputPower_watts*1.E7/(omega/1.*(m*c*c/e)*(m*c*c/e)*(Current_ampers/17045.)*(1./(k1*10))); // = S a
																////////
	double F0, F1, A0, A1;
	double xi = 0.0001;

	printf("Qa1 = %g,\t G1 = %g,\t inputPowerDimensionless = %g \n", Qa1, G1, inputPower_dimensionless);
	
	for(int i = 0; i < N_it; i++)
	{
		
		A0 = nextA;
		F0 = A0*A0 - 2.*Qa1*G1*(inputPower_dimensionless + DeltaEnergy(A0)); 
	
		A1 = nextA*(1.+xi);
		F1 = A1*A1 - 2.*Qa1*G1*(inputPower_dimensionless + DeltaEnergy(A1));

		if(A0*A1 < 0) {xi *= 0.5; continue;}
		if(fabs(F1 - F0) < 1e-14) break;
		
		nextA = A0- F0/(F1- F0)*(A1 - A0);

		printf("nextA  = %g\t Residual = %g\n", nextA, F1 - F0);
						
	}	
	printf("A_{stat}  = %g\n", nextA);
	A_stat = nextA;
	
	cplx A00, A01, A10;
	cplx F00, F01, F10;

	double dX, dY, J11, J12, J21, J22, det;

	xi = 0.00001;

	A00 = nextA;

	double detuning  = 0;

	for(int i = 0; i < 10; i++)
	{
		
		F00 = -A00*(1./(2.*Qa2) + I*detuning) + G2*(CurrentA2(nextA, real(A00), imag(A00)));

		A01 = A00+xi*abs(A00);
		F01 = -A01*(1./(2.*Qa2) + I*detuning) + G2*(CurrentA2(nextA, real(A01), imag(A01)));

		A10 = A00+I*xi*abs(A00);
		F10 = -A10*(1./(2.*Qa2) + I*detuning) + G2*(CurrentA2(nextA, real(A10), imag(A10)));


		J11 = real(F01 - F00)/(xi*abs(A00));  J12 = real(F10 - F00)/(xi*abs(A00)); 
		J21 = imag(F01 - F00)/(xi*abs(A00));  J22 = imag(F10 - F00)/(xi*abs(A00)); 

		det = J11*J22 - J12*J21;

		dX = -( J22*real(F00)-J12*imag(F00))/det;
		dY = -(-J21*real(F00)+J11*imag(F00))/det;

		if((fabs(dX) < 1e-20)&&(fabs(dY) < 1e-20)) break;

		cplx temp = real(A00)+dX + I*(imag(A00) + dY);

		A00 = temp;

		printf("A2_{stat}  = %g\tdet = %g\t%g\t%g\n", abs(A00), det, dX, dY);
		
	}
	printf("A2_{stat}  = %g\t Residual = %g\n", abs(A00),abs(A10 - A00));

	return A00;

}

double MultiplierTreeCavity::retriveBPower(cplx B, double omega, double Qb, double NormaB)
{
	if (Q_difr < 300) Q_difr = 300;
	double Qfull = Q_difr*Qb / (Q_difr + Qb);
	double electronPower = pow(m*c*c / e, 2)*(3.*omega) / Qfull*1e-7*pow(abs(B)*10., 2)*NormaB*1e-3*0.5;
	double outputPower = electronPower*Qb / (Qb + Q_difr);
	return outputPower;
}

cplx MultiplierTreeCavity::FindBstatDoubleScheme(cplx b0, int Nb_it, double Gb, double Astat, cplx A2_stat)
{

	double xi = 0.001;
	cplx F00, F01, F10, B00, B01, B10; 
	double dX, dY, J11, J12, J21, J22, det;

//	B00 = 1.57e-5*(1.+0.01*I); Nb_it = 15;
	B00 = b0;

/*	
	double temp_re[6000];
	double temp_im[6000];
	Ld = 15.;
	cplx deb1 = CurrentB2(real(B00), imag(B00), Astat, 0.);


	printf("\t CurrentB2(B00, Astat, 0) = %g, %g\n", real(deb1), imag(deb1));
	
	cplx deb2 = CurrentB(real(B00), imag(B00), Astat);
	printf("\t CurrentB(B00, Astat) = %g, %g\n", real(deb2), imag(deb2));
			cudaMemcpy(temp_re, d_rJ3, sizeof(double)*3100, cudaMemcpyDeviceToHost);

		FILE *file = fopen("D:\\Piotr\\Multiplier_260GHz_Data\\B1_re.txt", "w");
		for(int pr = 0; pr < 3100; pr++) fprintf(file, "%g\t%g\n",(double)pr*(93)/6000, temp_re[pr]/double(NP*NQ*NS));
		fclose(file);

		cudaMemcpy(temp_im, d_iJ3, sizeof(double)*3100, cudaMemcpyDeviceToHost);

		file = fopen("D:\\Piotr\\Multiplier_260GHz_Data\\B1_im.txt", "w");
		for(int pr = 0; pr < 3100; pr++) fprintf(file, "%g\t%g\n",(double)pr*(93)/6000, temp_im[pr]/double(NP*NQ*NS));
		fclose(file);
	
*/

 
	for(int i = 1; i < Nb_it; i++)
	{
		printf("B00:\t%g\t%g\n", real(B00), imag(B00));

		F00 = -B00/(2.*Qb) + Gb*CurrentB2(real(B00), imag(B00), Astat, A2_stat);

	//	printf("current B2 returned %g \n", abs(CurrentB2(real(B00), imag(B00), Astat, A2_stat)));

		B01 = B00+xi*abs(B00);
		F01 = -B01/(2.*Qb) + Gb*CurrentB2(real(B01), imag(B01), Astat, A2_stat);

		B10 = B00+I*xi*abs(B00);
		F10 = -B10/(2.*Qb) + Gb*CurrentB2(real(B10), imag(B10), Astat, A2_stat);


		J11 = real(F01 - F00)/(xi*abs(B00));  J12 = real(F10 - F00)/(xi*abs(B00)); 
		J21 = imag(F01 - F00)/(xi*abs(B00));  J22 = imag(F10 - F00)/(xi*abs(B00)); 

		det = J11*J22 - J12*J21;

		dX = -( J22*real(F00)-J12*imag(F00))/det;
		dY = -(-J21*real(F00)+J11*imag(F00))/det;

		cplx temp = real(B00)+dX + I*(imag(B00) + dY);

		B00 = temp;
	
		printf("|B|  = %g\t%g\t%g\tdet = %g\n", abs(B00), real(B00), imag(B00), det);
	}
	return B00;
}

void MultiplierTreeCavity::PrintParamsDoubleScheme(char *filename, char *comment)
{
	FILE *file2=  fopen(filename, "w"); 
	
	
	fprintf(file2, "\"Current\"\t%g\n\"period\"\t%g\n\"Nper\"\t%i\n", Current_ampers, period, Nperiods);
	fprintf(file2, "\"ld1\"\t%g\n\"la2\"\t%g\n\"kw\"\t%g\n", Ld1, La2, k1);
	fprintf(file2, "\"ld2\"\t%g\n\"lb\"\t%g\n", Ld2, Lb);
	fprintf(file2, "\"norm1\"\t%g\n\"norm2\"\t%g\n\"p_in\"\t%g\n\"delta_freq\"\t%g\n", Norma, Norma1b, inputPower_watts, delta_freq);
	fprintf(file2, "\"voltage\"\t%g\n\"Qa1\"\t%g\n\"Qa2\"\t%g\n\"Qb\"\t%g\n", voltage, Qa1, Qa2, Qb);
	fprintf(file2, "\"beam thickness\"\t%g\n", wall);
	fprintf(file2, "\"%s\"\n", comment);

	fclose(file2);
}