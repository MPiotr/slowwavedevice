#include "twt.h"

#ifndef __BWO_H
#define __BWO_H

template  <class TWTsolver>
class BWO : public TWTsolver
{
	double a0;
	double delta;
	double reflectionCoeff, reflectionPhase;
	int parasiteIndex;
	virtual void iterate(int paramsInd);
	void printResults(FILE *file);
	void printResultsStartCurrent(FILE *file);
	void printDataHeader(FILE *file);
	void printParamsHeader(FILE *file);
	void printResultHeader(FILE *file);
	void printStCurrResultHeader(FILE *file);

public:
	BWO(QDomDocument *doc);
	//	BWO_2D(QDomDocument *doc, int parasiteIndex, char *fieldFile);
	//  BWO_2D(QDomDocument *doc, TWT_2D *twt, int parasiteIndex, char *fieldFile);
	void solveBWO();
	void solveBWO(double A, double delta);
	double findStartCurrent();
	double findStartCurrent(double G0, double *delta0);
	double findStartCurrent(double I0, double Delta0, double fre, double synch_angle, double group_speed, double Q);
	double getNormalizedLength(double curr, double *delta);
	double currentFromNormalizedLength(double normLength, double *Delta0);
};
#endif