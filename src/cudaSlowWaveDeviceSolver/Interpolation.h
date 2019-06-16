#ifndef INTERPOLATION_H
#define INTERPOLATION_H
#include <vector>


class Interpolation
{
	std::vector<std::pair<double, double> > y;

	void load(double *x, double *y, int n);
	void load(char *filename);
	double xmin;
	double xmax;

public:
	Interpolation();
	Interpolation(char *filename);
	Interpolation(double  *X, double *Y, int  N);
	double at(double x);
	double der(double x);

	double xMin();
	double xMax();

	void reload(char * filename) { 
		y.clear(); 
		load(filename); 
	}

	void reload(double *_x, double *_y, int n) {
		y.clear();
		load(_x, _y, n);
	}

	void unitTesting();
};
typedef Interpolation* func;
#endif