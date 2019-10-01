#include "fieldcalculatorcontroller.h"

FieldCalculatorController::FieldCalculatorController(projectviewer * project, QTextBrowser * logger, QObject * parent)
	: proj(project), logBrowser(logger), QObject(parent)
{
	fieldCalculator.setProgram("FreeFem++");

	QProcessEnvironment env = QProcessEnvironment::systemEnvironment();
	QString path = env.value("Path");
	env.insert("FF_INCLUDEPATH", "C:\\Users\\Piotr\\Documents\\GitHub\\slowwavedevice\\bin\\x64\\Release;");

	fieldCalculator.setProcessEnvironment(env);

	connect(&fieldCalculator, static_cast<void(QProcess::*)(int, QProcess::ExitStatus)> (&QProcess::finished),
		[=](int exitCode, QProcess::ExitStatus exitStatus) { processFinished(exitCode, exitStatus); }
	);
}

void FieldCalculatorController::calculate()
{
	fieldCalculator.setArguments({ "fieldCalculator.edp" });

	QFile output("field_param.txt");
	if (!output.open(QIODevice::WriteOnly)) { errorExit(QString("can't open 'param.txt' for writing"));	return; }
	if (!proj->savePeriodParamsForFieldCalculation(&output)) { errorExit(QString("saving 'param.txt'"));	                return; }
	output.close();

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

void FieldCalculatorController::processFinished(int exitCode, QProcess::ExitStatus stat) 
{
	logBrowser->append(fieldCalculator.readAllStandardOutput());
	if (exitCode != 0) { errorExit(QString("exit code is not zero"));	return; }
	emit finished(0);
}

void FieldCalculatorController::errorExit(QString & message)
{
	logBrowser->append("<b><font color = \"red\">Error: </font></b>" + message);
	emit finished(1);
}



