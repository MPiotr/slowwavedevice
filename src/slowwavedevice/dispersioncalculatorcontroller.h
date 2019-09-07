#pragma once

#include <QObject>
#include <qprocess.h>
#include "projectviewer.h"

class DispersionCalculatorController : public QObject
{
	Q_OBJECT

private:
	QProcess resonatorGenerator;
	QProcess dispersionCalculator;
	QProcess fileMerger;
	QProcess modeSelector;
	QProcess *activeProcess = nullptr;

	QTextBrowser *logBrowser;

	projectviewer *proj;
	void errorExit(QString &message);

private slots:
	
	void resonatorGeneratorFinished(int exitCode, QProcess::ExitStatus stat);
	void dispersionCalculatorFinished(int exitCode, QProcess::ExitStatus stat);
	void fileMergerFinished(int exitCode, QProcess::ExitStatus stat);
	void modeSelectorFinished(int exitCode, QProcess::ExitStatus stat);
	void readConsole();

public:
	DispersionCalculatorController(projectviewer* project, QTextBrowser *loggerBrowser,  QObject *parent);
	~DispersionCalculatorController();

	void calculate();
	void abort();

signals:
	void finished(const int &exitStatus);


};

