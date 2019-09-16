#include "fieldcalculatorcontroller.h"

FieldCalculatorController::FieldCalculatorController(projectviewer * project, QTextBrowser * logger, QObject * parent)
	: proj(project), logBrowser(logger), QObject(parent)
{
	fieldCalculator.setProgram("FreeFem++");

	connect(&fieldCalculator, static_cast<void(QProcess::*)(int, QProcess::ExitStatus)> (&QProcess::finished),
		[=](int exitCode, QProcess::ExitStatus exitStatus) { processFinished(exitCode, exitStatus); }
	);
}

void FieldCalculatorController::calculate()
{
	fieldCalculator.setArguments({ "fieldCalculator.edp" });
	fieldCalculator.start();
}


void FieldCalculatorController::abort()
{
	fieldCalculator.kill();
}

void FieldCalculatorController::readConsole()
{
	logBrowser->append(fieldCalculator.readAllStandardOutput());
}

void FieldCalculatorController::processFinished(int exitCode, QProcess::ExitStatus stat) {
	if (exitCode != 0) { errorExit(QString("exit code is not zero"));	return; }
	emit finished(0);
}

void FieldCalculatorController::errorExit(QString & message)
{
	logBrowser->append("<b><font color = \"red\">Error: </font></b>" + message);
	emit finished(1);
}



