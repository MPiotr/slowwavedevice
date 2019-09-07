#include <filesystem>
#include "dispersioncalculatorcontroller.h"

DispersionCalculatorController::DispersionCalculatorController(projectviewer* project, QTextBrowser* logger, QObject *parent)
	: proj(project), logBrowser(logger), QObject(parent)
{
	resonatorGenerator.setProgram("resonatorGenerator.exe");
	
	dispersionCalculator.setProgram("mpiexec");
//    dispersionCalculator.setEnvironment({ "PATH=.:$PATH" });
	

	fileMerger.setProgram("file_merger.exe");
	fileMerger.setArguments({ "dispFem", "1", "2" });

	modeSelector.setProgram("mode_selector.exe");

	connect(&resonatorGenerator, static_cast<void(QProcess::*)(int, QProcess::ExitStatus)> (&QProcess::finished),
		    [=](int exitCode, QProcess::ExitStatus exitStatus) { resonatorGeneratorFinished(exitCode, exitStatus); } 
	       );

	connect(&dispersionCalculator, 	static_cast<void(QProcess::*)(int, QProcess::ExitStatus)> (&QProcess::finished),
		    [=](int exitCode, QProcess::ExitStatus exitStatus) {dispersionCalculatorFinished(exitCode, exitStatus); } 
	       );

	connect(&fileMerger, static_cast<void(QProcess::*)(int, QProcess::ExitStatus)> (&QProcess::finished),
		    [=](int exitCode, QProcess::ExitStatus exitStatus) {fileMergerFinished(exitCode, exitStatus); } 
	        );

	connect(&modeSelector, static_cast<void(QProcess::*)(int, QProcess::ExitStatus)> (&QProcess::finished),
	    	[=](int exitCode, QProcess::ExitStatus exitStatus) {modeSelectorFinished(exitCode, exitStatus); }
	       );

	connect(&resonatorGenerator,   SIGNAL(readyReadStandardOutput()), this, SLOT(readConsole()));
	connect(&dispersionCalculator, SIGNAL(readyReadStandardOutput()), this, SLOT(readConsole()));
	connect(&fileMerger,           SIGNAL(readyReadStandardOutput()), this, SLOT(readConsole()));
	connect(&modeSelector,         SIGNAL(readyReadStandardOutput()), this, SLOT(readConsole()));
}

DispersionCalculatorController::~DispersionCalculatorController()
{
}

void DispersionCalculatorController::calculate()
{
	QString shapeFileName;
	if (!proj->getShapeFileName(&shapeFileName)) {
		logBrowser->append("<b><font color = \"red\">Error</font></b>: Shape file name is not specified");
		return;
	}
	
	proj->recalculatePeriodFromShape(logBrowser);

	resonatorGenerator.setArguments({ shapeFileName });
	resonatorGenerator.start();
	activeProcess = &resonatorGenerator;
	logBrowser->append(resonatorGenerator.readAllStandardOutput());	
}



void DispersionCalculatorController::resonatorGeneratorFinished(int exitCode, QProcess::ExitStatus stat) 
{
	if (exitCode != 0) { errorExit(QString("resonatorGenerator exit code is not zero"));	return; }

	QFile output("param.txt");
	if (!output.open(QIODevice::WriteOnly)) { errorExit(QString("can't open 'param.txt' for writing"));	return; }
	if (!proj->savePeriodParamsForDispersionCalculation(&output)) { errorExit(QString("saving 'param.txt'"));	                return; }
	output.close();

	int numCores = proj->getNumCores();
	dispersionCalculator.setArguments({ "-n", QString::number(numCores), "FreeFem++-mpi" ,"dispersionCalculator.edp" });
	dispersionCalculator.start();
	activeProcess = &dispersionCalculator;
}

void DispersionCalculatorController::dispersionCalculatorFinished(int exitCode, QProcess::ExitStatus stat)
{
	readConsole();
	if (exitCode != 0) { 		
		errorExit(QString("dispersionCalculator exit code is not zero")); 
		return; 
	}
	fileMerger.start();
	activeProcess = &fileMerger;
}

void DispersionCalculatorController::fileMergerFinished(int exitCode, QProcess::ExitStatus stat)
{
	readConsole();
	if (exitCode != 0) { errorExit(QString("fileMerger exit code is not zero"));	return; }
	char buf[30], buf1[30];
	double per = proj->getPeriod();
	double freq = proj->getRefFreq();
	sprintf(buf, "%g", per);
	sprintf(buf1, "%g", freq);
	modeSelector.setArguments({ "dispFem.csv", "1", "2", "3", "dispersion", "0.8", buf, buf1 });
	modeSelector.start();
	activeProcess = &modeSelector;
}

void DispersionCalculatorController::modeSelectorFinished(int exitCode, QProcess::ExitStatus stat)
{
	readConsole();
	if (exitCode != 0) { errorExit(QString("modeSelector exit code is not zero"));	return; }
	QString filename = proj->getDispersionFileName();
	if (filename.isEmpty())
		filename = "dispersion.csv";
	std::experimental::filesystem::rename("dispersion0.csv", filename.toStdString());
	activeProcess = nullptr;
	emit finished(0);

}

void DispersionCalculatorController::readConsole()
{
	logBrowser->append(activeProcess->readAllStandardOutput());
}
void DispersionCalculatorController::errorExit(QString & message)
{
	logBrowser->append("<b><font color = \"red\">Error: </font></b>" + message);
	activeProcess = nullptr;
	emit finished(1);
}

void DispersionCalculatorController::abort()
{
	if (activeProcess != nullptr)
		activeProcess->kill();
}
