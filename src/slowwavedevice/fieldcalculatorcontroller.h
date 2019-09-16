#pragma once

#pragma once

#include <QObject>
#include <qprocess.h>
#include "projectviewer.h"

class FieldCalculatorController : public QObject
{
	Q_OBJECT
	
	QProcess fieldCalculator;
	QTextBrowser *logBrowser;
	projectviewer *proj;

	void errorExit(QString & message);

public:
	FieldCalculatorController(projectviewer* project, QTextBrowser *loggerBrowser, QObject *parent);
	void calculate();
	void abort();

signals:
	void finished(int exitStatus);

private slots:
	void processFinished(int exitCode, QProcess::ExitStatus stat);
	void readConsole();
	
};