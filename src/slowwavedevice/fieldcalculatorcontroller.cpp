#include "fieldcalculatorcontroller.h"

FieldCalculatorController::FieldCalculatorController(projectviewer * project, QTextBrowser * logger, QObject * parent)
	: proj(project), logBrowser(logger), QObject(parent)
{
	fieldCalculator.setProgram("FreeFem++");

	QProcessEnvironment env = QProcessEnvironment::systemEnvironment();
	QString path = env.value("Path");
	env.insert("FF_INCLUDEPATH", QCoreApplication::applicationDirPath()+";");

	fieldCalculator.setProcessEnvironment(env);

	connect(&fieldCalculator, static_cast<void(QProcess::*)(int, QProcess::ExitStatus)> (&QProcess::finished),
		[=](int exitCode, QProcess::ExitStatus exitStatus) { processFinished(exitCode, exitStatus); }
	);
	connect(&fieldCalculator, SIGNAL(readyReadStandardOutput()), this, SLOT(readConsole()));
}

bool FieldCalculatorController::getFemScriptName(QString &scriptName) {
	QString shapeType = proj->getShapeType();
	if (shapeType.isEmpty()) { errorExit(QString("shapeType is not specified")); return false; }
	if (shapeType != "planar" && shapeType != "axial") { errorExit(QString("shapeType is not recognized")); return false; }

	if (shapeType == "planar") scriptName = "fieldCalculator.edp";
	if (shapeType == "axial") scriptName = "fieldCalculatorAxial.edp";
	return true;
}

void FieldCalculatorController::calculate()
{
	QString scriptName;
	if(!getFemScriptName(scriptName)) return;
	fieldCalculator.setArguments({ scriptName });

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

	QString filename = proj->getFieldFileName();
	if (filename.isEmpty())
		filename = "tfield.csv";
	//	std::experimental::filesystem::rename("dispersion0.csv", filename.toStdString());

	remove(filename.toLocal8Bit().data());
	int renres = rename("tFieldFem.csv", filename.toLocal8Bit().data());

	emit finished(0);
}

void FieldCalculatorController::errorExit(QString & message)
{
	logBrowser->append("<b><font color = \"red\">Error: </font></b>" + message);
	emit finished(1);
}



