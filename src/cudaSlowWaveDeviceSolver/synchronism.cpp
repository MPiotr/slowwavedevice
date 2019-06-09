#include "synchronism.h"

#define M_PI  3.141592653589793

double Synchronism::frequency()
{
	return _frequency;
}
double Synchronism::frequency(double voltage)
{
	if(voltage != _voltage) setVoltage(voltage);
	return _frequency;
}
double Synchronism::frequency(double voltage, int n)
{
	if (voltage != _voltage) setVoltage(voltage);
	if (n >= frequencies.size()) return -1; 
	return frequencies.at(n);
}
double Synchronism::groupSpeed()
{
	return _groupSpeed;
}
double Synchronism::groupSpeed(double voltage)
{
	if (voltage != _voltage) setVoltage(voltage);
	return _groupSpeed;
}
double Synchronism::groupSpeed(double voltage, int n)
{
	if (voltage != _voltage) setVoltage(voltage);
	return _groupSpeed;
}
double Synchronism::wavenumber()
{
	return wavenumbers[0];
}
double Synchronism::wavenumber(double voltage)
{
	if (voltage != _voltage) setVoltage(voltage);
	return wavenumbers[0];
}
double Synchronism::wavenumber(double voltage, int n)
{
	if (n >= wavenumbers.size())  return -1;
	if (voltage != _voltage) setVoltage(voltage);
	return wavenumbers[n];
}


void Synchronism::setVoltage(double voltage)
{
	double gamma = 1 + voltage / 511.;
	double beta = sqrt(gamma*gamma - 1) / gamma *299.8/(2*180*_d);//дисперсионка -- зависимость частоты от градусов

	int startmode = 0;
	if (_inputFreq != 0) startmode = 1;

	for (int i = startmode; i < Nmodes; i++)
	{
		wavenumbers[i] = findRootFixedBeta(modes[i], beta, wavenumbers[i]);
		frequencies[i] = modes[i]->at(wavenumbers[i]);
	}
	if (startmode == 0){
		_frequency = frequencies[0];
		_groupSpeed = (2 * _d * 180) / 299.8* (modes[0])->der(wavenumbers[0]);
	}

}
void Synchronism::setFrequency(double fre, double h0)
{
	_inputFreq = fre;
	double groupspeed;
	for (int i = 0; i < 1; i++)
	{
		wavenumbers[i] = findRootFixedFre(modes[i], _inputFreq, (h0 == 0) ? wavenumbers[i]: h0);
		groupspeed = (2 * _d * 180) / 299.8* (modes[0])->der(wavenumbers[0]);
		if (groupspeed < 0) wavenumbers[i] = findRootFixedFre(modes[i], _inputFreq, wavenumbers[i] + 360);
		frequencies[i] = modes[i]->at(wavenumbers[i]);
	}
	printf("%g\n", wavenumbers[0]);

	_frequency = frequencies[0];
	_groupSpeed = (2 * _d * 180) / 299.8* (modes[0])->der(wavenumbers[0]);
}

double Synchronism::findRootFixedFre(Interpolation *mode, double inputFreq, double h0)
{
	double h = h0;
	for (int i = 0; i < 100; i++)
	{
		if (mode->der(h) != 0){
			double f = (mode->at(h) - _inputFreq);
			double df = (mode->der(h));
			if (fabs(f / df) / fabs(h) < 1e-7) break;
			h -=  f/df;
		}
		else
			return h;
	}
	return h;
}
double Synchronism::findRootFixedBeta(Interpolation *mode, double beta, double h0)
{
	double h = h0;
	for (int i = 0; i < 150; i++)
	{
		if (mode->der(h) != 0){
			double f = (mode->at(h) - beta*h);
			double df = (mode->der(h)-beta);
			if (fabs(f / df) / fabs(h) < 1e-10) break;
			h -= f / df;
		}
		else
			return h;
	}
	return h;
}
void Synchronism::unitTesting()
{
	Interpolation mode1("F:\\Piotr\\CalcData\\twt_data\\testingDisp.txt");
	mode1.unitTesting();
	modes.push_back(&mode1);
	frequencies.push_back(mode1.at(0));
	wavenumbers.push_back(0);
	_d = 1.066;
	Nmodes = 1;

	setVoltage(23);
	printf("voltage = 23, frequency = %g, h = %g, groupSpeed = %g\n", _frequency, wavenumbers[0], _groupSpeed);
	printf("voltage = 24, frequency = %g, h = %g, groupSpeed = %g\n", frequency(24), wavenumbers[0], _groupSpeed);
	printf("voltage = 25, frequency = %g, h = %g, groupSpeed = %g\n", frequency(25), wavenumbers[0], _groupSpeed);
	setVoltage(17);
	printf("voltage = 17, frequency = %g, h = %g, groupSpeed = %g\n", frequency(17.), wavenumbers[0], _groupSpeed);

	

}
void Synchronism::addMode(char *filename)
{
	Interpolation *mode1 = new Interpolation(filename);
	modes.push_back(mode1);
	frequencies.push_back(mode1->at(400));
	wavenumbers.push_back(400);
	Nmodes++;
}
Synchronism::Synchronism()
{
	_d = 1.;
	_inputFreq = 0;
}
Synchronism::Synchronism(char *filename)
{
	Interpolation *mode1 = new Interpolation(filename);
	modes.push_back(mode1);
	frequencies.push_back(mode1->at(400));
	wavenumbers.push_back(400);
	_d = 1.;
	_inputFreq = 0;
	Nmodes = 1;

}
Synchronism::Synchronism(char *filename, double per, double inputFre)
{
	Interpolation *mode1 = new Interpolation(filename);
	modes.push_back(mode1);
	_inputFreq = inputFre;
	frequencies.push_back(_inputFreq);
	wavenumbers.push_back(330);
	_d = per;
	Nmodes = 1;

}
