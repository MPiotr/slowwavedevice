#include "Interpolation.h"
#include <algorithm>

void Interpolation::load(double *X, double *Y, int N)
{
	for (int i = 0; i < N; i++)
	{
		y.push_back(std::pair<double, double>(X[i], Y[i]));
	}	
	std::sort(y.begin(), y.end());
	xmin = y.front().first;
	xmax = y.back().first;

}
void Interpolation::load(char *filename)
{
	FILE *file = fopen(filename, "r");
	if (file == NULL) { printf("Error:  File %s not found!\n"); throw("file not found"); }
	float x, fx, fy;
	while (fscanf(file, "%g,%g,%g\n", &x, &fx, &fy) >= 2)
	{
		y.push_back(std::pair<double, double>(x, fx));
	}
	std::sort(y.begin(), y.end());
	xmin = y.front().first;  
	xmax = y.back().first;
	fclose(file);

}

double Interpolation::at(double x)
{
	if (y.size() == 0) throw("Interpolation data table is not loaded");

	std::pair<double, double> key(x, 0);
	std::vector<std::pair<double, double>>::iterator nplus  = std::upper_bound(y.begin(), y.end(), key);
	std::vector<std::pair<double, double>>::iterator nminus;

	
	if (nplus == y.end()) 
	{		--nplus; nminus = prev(nplus);		}
	else if (nplus == y.begin())
	{
		++nplus; nminus = y.begin();
	}
	else
		nminus = prev(nplus);
	
	

	double xplus = nplus->first;
	double xminus = nminus->first;

	double coefMinus = (nplus->first - x) / (nplus->first - nminus->first);
	double coefPlus = (x - nminus->first) / (nplus->first - nminus->first);


	return coefPlus*nplus->second + coefMinus*nminus->second;

}

double Interpolation::der(double x)
{
	if (y.size() == 0) throw("Interpolation data table is not loaded");

	std::pair<double, double> key(x, 0);
	std::vector<std::pair<double, double>>::iterator nplus = std::upper_bound(y.begin(), y.end(), key);
	std::vector<std::pair<double, double>>::iterator nminus;

	if (nplus == y.end())
	{
		--nplus; nminus = prev(nplus);
	}
	else if (nplus == y.begin())
	{
		++nplus; nminus = y.begin();
	}
	else
		nminus = prev(nplus);

	return (nplus->second - nminus->second) / (nplus->first - nminus->first);
}

void Interpolation::unitTesting()
{
//	load("F:\\Piotr\\CalcData\\twt_data\\testingDisp.txt");

	FILE *file = fopen("F:\\Piotr\\CalcData\\twt_data\\interpolationTest.txt", "w");
	for (int i = 0; i < 100; i++)
	{
		double x =  800. / 100.*double(i);
		fprintf(file, "%g,%g,%g\n", x, at(x), der(x));
	}
	fclose(file);


}
double Interpolation::xMin()
{
	return xmin;
}
double Interpolation::xMax()
{
	return xmax;
}
Interpolation::Interpolation(char *name)
{
	load(name);
}
Interpolation::Interpolation(double  *X, double *Y, int  N)
{
	load(X, Y, N);
}
Interpolation::Interpolation()
{

}