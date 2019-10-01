#include <qfiledialog.h>
#include "slowwavedevice.h"
#include "customxmldelegate.h"

//TODO: add calculation of the 0D wave norma from FreeFem

slowwavedevice::slowwavedevice(QWidget *parent)
	: QMainWindow(parent), sharedMemory("slowwavedevicecontroller"), sharedMemoryResult("slowwavedeviceresults"), resultStream(&resultBuffer)
{

	ui.setupUi(this);
	about.setupUi(&aboutDialog);
	for (int i = 0; i < 5; ++i) {
		recentFileActs[i] = ui.menu->addAction(QString(), this, SLOT(openRecentFile()));
		recentFileActs[i]->setVisible(false);
	}

	projmodel = new projectviewer(); //stays before loadRecentFiles, since the latter changes cwd
	loadRecentFiles();
	disableButtons();
	ui.AbortAction->setEnabled(false);
	ui.ExitAction->setEnabled(true);
	ui.OpenAction->setEnabled(true);
	
	qplot = new QCustomPlot(this);
	couplingPlot = new QCustomPlot(this); 
	resultPlot = new QCustomPlot(this);

	couplingPlot->setVisible(false);
	resultPlot->setVisible(false);

	resultStorageX = new double[3000];
	resultStorageY = new double[3000];

	ui.gridLayout_2->addWidget(qplot, 0, 1, 1, 1);
	ui.gridLayout_2->addWidget(couplingPlot, 0, 1, 1, 1);
	ui.gridLayout_2->addWidget(resultPlot, 0, 1, 1, 1);


	ui.treeView->setModel(projmodel);
	ui.treeView->setItemDelegate(new customxmldelegate (projmodel, &(projmodel->toolTipsXML)));
	ui.treeView->setContextMenuPolicy(Qt::CustomContextMenu);
	ui.menuBar->show();
	ui.textBrowser->setFontFamily("Consolas");

	solverProcess.setProgram("\"cudaSlowWaveDeviceSolver.exe\"");
	QStringList args; args.append("control");
	solverProcess.setArguments(args);
	whoIsActive = none;

	dispCalcContoller = new DispersionCalculatorController(projmodel, ui.textBrowser, this);
	fieldCalcController = new FieldCalculatorController(projmodel, ui.textBrowser, this);

	connect(ui.OpenAction, SIGNAL(triggered()), this, SLOT(loadFile()));
	connect(ui.SaveAction, SIGNAL(triggered()), this, SLOT(saveFile()));
	connect(ui.SaveAsAction, SIGNAL(triggered()), this, SLOT(saveFileAs()));
	connect(ui.ExitAction, SIGNAL(triggered()), this, SLOT(exit()));
	connect(ui.AbortAction, SIGNAL(triggered()), this, SLOT(abortSolver()));
	connect(ui.action, SIGNAL(triggered()), this, SLOT(start()));
	connect(ui.CalculateDispersionAction, SIGNAL(triggered()), this, SLOT(calculateDispersion()));
	connect(ui.CalculateField, SIGNAL(triggered()), this, SLOT(calculateField()));
	connect(ui.AboutAction, SIGNAL(triggered()), this, SLOT(showAbout()));
	connect(&solverProcess, static_cast< void(QProcess::*)(int, QProcess::ExitStatus) > (&QProcess::finished),
		[=](int exitCode, QProcess::ExitStatus exitStatus){
		solverFinished(exitCode, exitStatus);
	}
	);
	connect(&solverProcess, SIGNAL(readyReadStandardOutput()), this, SLOT(readConsole()));
	connect(dispCalcContoller, SIGNAL(finished(int)), this, SLOT(dispersionCalcFinished(int)));
	connect(fieldCalcController, SIGNAL(finished(int)), this, SLOT(dispersionCalcFinished(int)));
	connect(projmodel, SIGNAL(voltageChanged(QString)), this, SLOT(replotVoltage(QString)));
	connect(projmodel, SIGNAL(setVisiblePlot(int)), this, SLOT(showPlot(int)));
	connect(ui.treeView , SIGNAL(clicked(QModelIndex )), projmodel, SLOT(itemClicked(QModelIndex )));
	connect(ui.treeView, SIGNAL(customContextMenuRequested(QPoint)), this, SLOT(requestModelContextMenu(QPoint)));
	connect(this, SIGNAL(projmodelContextMenu(QTreeView *, QModelIndex, QPoint)), projmodel, SLOT(showContextMenu(QTreeView *, QModelIndex, QPoint)));
	connect(&timer, SIGNAL(timeout()), this, SLOT(readResults()));

}

slowwavedevice::~slowwavedevice()
{
	saveRecentFiles();
}

void slowwavedevice::loadFile()
{
	QDir::setCurrent(lastDirectory);
	inputFile = QFileDialog::getOpenFileName(this, tr("Open input file"), lastDirectory, tr("input file (*.xml)"));

	openFileFun(inputFile);

}
void slowwavedevice::saveFileAs()
{
	QFileInfo fileinfo(inputFile);
	QDir::setCurrent(fileinfo.path());
	QString outputFile = QFileDialog::getSaveFileName(this, tr("Save file"), fileinfo.path(), tr("input file (*.xml)"));

	QFile output(outputFile);
	output.open(QIODevice::WriteOnly);
	QTextStream out(&output);

	QString string = projmodel->getDoc().toString();
	out << string;
	bool flush = output.flush();
	output.close();
	
	if (flush){
		ui.textBrowser->append("<div align = \"left\">" + outputFile + "  saved </div>");
		recentFiles.push_front(outputFile);
		if (recentFiles.size() > 5) recentFiles.pop_back();
	}
	else
		ui.textBrowser->append("<div align = \"left\"> error while saving </div>");

}
void slowwavedevice::saveFile()
{

	QFile output(inputFile);
	output.open(QIODevice::WriteOnly);
	QTextStream out(&output);

	QString string = projmodel->getDoc().toString();
	out << string;
	bool flush = output.flush();
	output.close();

	if (flush)
	    ui.textBrowser->append("<div align = \"left\">" + inputFile + "  saved </div>");
	else
		ui.textBrowser->append("<div align = \"left\"> error while saving </div>");
}

void slowwavedevice::start()
{
	QBuffer buffer;
	buffer.open(QBuffer::ReadWrite);
	QDataStream out(&buffer);
	disableButtons();
	QString tmp = projmodel->getXMLString();
	out << projmodel->getXMLString();  
	  
	int size = buffer.size();  
	if (sharedMemory.isAttached())	sharedMemory.detach();
	if(!sharedMemory.create(size)) ui.textBrowser->append("<font color = \"red\">Shared memory segment is not created\n<\font>");
	
	sharedMemory.lock();
	char *to = (char*)sharedMemory.data();
	const char *from = buffer.data().data();
	memcpy(to, from, qMin(sharedMemory.size(), size));
	sharedMemory.unlock();

	solverProcess.setProcessChannelMode(QProcess::MergedChannels);
	while (solverProcess.waitForReadyRead())
		ui.textBrowser->insertPlainText(solverProcess.readAll());

	showPlot(3);
	solverProcess.start();
	whoIsActive = solver;
	timer.start(500);
	
	ui.textBrowser->append(solverProcess.readAllStandardOutput());
	
}

void slowwavedevice::calculateDispersion() 
{
	disableButtons();
	dispCalcContoller->calculate();
	whoIsActive = dispersion;
}

void slowwavedevice::calculateField()
{
	disableButtons();
	fieldCalcController->calculate();
	whoIsActive = field;
}

void slowwavedevice::abortSolver()
{
	if (whoIsActive == solver) {
		solverProcess.kill();
		timer.stop();
	}
	if (whoIsActive == dispersion) {
		dispCalcContoller->abort();
	}
	if (whoIsActive == field) {
		fieldCalcController->abort();
	}
	whoIsActive = none;
	enableButtons(); 
	ui.textBrowser->append("<font color = \"red\">Aborted by user<\font>");
}

void slowwavedevice::solverFinished(int exitCode, QProcess::ExitStatus stat)
{
	timer.stop();
	if (sharedMemoryResult.isAttached())	sharedMemoryResult.detach();
	enableButtons();
	ui.textBrowser->append("<font color = \"green\">finished<\font>");
}

void slowwavedevice::dispersionCalcFinished(int exitCode)
{
	enableButtons();
	whoIsActive = none;
	if (exitCode == 0) {
		ui.textBrowser->append("<b color = \"green\">Calculation has finished, reloading project<\b>");
		openFileFun(inputFile);
	}
	else {
		ui.textBrowser->append("<b color = \"red\">Calculation has failed<\b>");
	}
}

void slowwavedevice::readConsole()
{
	QProcess *process = nullptr;
	if (whoIsActive == solver) 	process = &solverProcess;
	if (process != nullptr) {
		ui.textBrowser->append(process->readAllStandardOutput());
	}
}
void slowwavedevice::replotVoltage(QString value)
{
	QVector<double> x(3), y(3);
	double voltage = value.toDouble();
	double period = projmodel->period;
	double gamma = 1 + voltage / 511.;
	double coef = 299.8 / (2.*M_PI)*sqrt(gamma*gamma - 1) / gamma*M_PI / (180.*period);
	
	for (int i = 0; i <= 2; i++)
	{
		x[i] = 720. / 2. * i;
		y[i] = x[i] * coef;
	}
	
	qplot->graph(1)->setData(x, y);
	qplot->replot();	
}
void slowwavedevice::showPlot(int n)
{
	if (n == 0) {
		couplingPlot->setVisible(false);
		resultPlot->setVisible(false);
		qplot->setVisible(true);
		qplot->replot();
	}
	if (n == 1)
	{
		resultPlot->setVisible(false);
		qplot->setVisible(false);
		couplingPlot->setVisible(true);
		couplingPlot->replot();
	}
	if (n == 2)
	{
		resultPlot->setVisible(true);
		qplot->setVisible(false);
		couplingPlot->setVisible(false);
	}
}
void slowwavedevice::requestModelContextMenu(QPoint point)
{
	QModelIndex index = ui.treeView->indexAt(point);
	emit projmodelContextMenu(ui.treeView,index, ui.treeView->mapToGlobal(point));
}
void slowwavedevice::exit()
{
	QApplication::quit();
}
void slowwavedevice::loadRecentFiles()
{
	lastDirectory = ".";
	cookie.setFileName("info.dat");
	cookie.open(QIODevice::ReadOnly);
	QTextStream info(&cookie);
	if (!info.atEnd()) lastDirectory = info.readLine();
	while (!info.atEnd())
	{
		recentFiles.push(info.readLine());
	}
	cookie.close();
	QDir::setCurrent(lastDirectory);

	updateRecentFilesActions();

	setWindowTitle("Slow-wave device solver");
}
void slowwavedevice::updateRecentFilesActions()
{
	for (int i = 0; i < recentFiles.size(); i++)
	{
		recentFileActs[i]->setText(recentFiles.at(i));
		recentFileActs[i]->setData(recentFiles.at(i));
		recentFileActs[i]->setVisible(true);
	}
	for (int i = recentFiles.size(); i < 5; i++)
	{
		recentFileActs[i]->setVisible(false);
	}

}
void slowwavedevice::saveRecentFiles()
{
	cookie.open(QIODevice::WriteOnly);
	QTextStream info(&cookie);
	info << lastDirectory << "\n";
	while (!recentFiles.empty())
	{
		info << recentFiles.first() << "\n";
		recentFiles.pop_front();
	}
	cookie.close();
}
void slowwavedevice::openRecentFile()
{
	if (const QAction *action = qobject_cast<const QAction *>(sender()))
	{
		openFileFun(action->data().toString());
		inputFile = action->data().toString();
	}
}
void slowwavedevice::openFileFun(QString inputFile)
{
	QFileInfo fileinfo(inputFile);
	QDir::setCurrent(fileinfo.path());
	lastDirectory = fileinfo.path();
	QFile input(inputFile);

	if (!input.open(QIODevice::ReadOnly)){
		ui.textBrowser->append("<b color = \"red\">Error opening file:  </b>" + inputFile);
		return;
	}
	setWindowTitle("Slow-wave device solver: " + inputFile);
	projmodel->setprojContent(&input, couplingPlot);
	ui.treeView->expandAll();
	projmodel->setDispersionsPlot(qplot, ui.textBrowser, 50);
	projmodel->setPlot(couplingPlot, "longitudinalStructureFile", 50, ui.textBrowser);
	qplot->replot();
	input.close();

	if (!recentFiles.contains(inputFile))
	{
		recentFiles.push_front(inputFile);
		if (recentFiles.size() > 5) recentFiles.pop_back();
	}
	else 
	{
		recentFiles.removeOne(inputFile);
		recentFiles.push_front(inputFile);
	}
	updateRecentFilesActions();

	ui.textBrowser->append("<div align = \"left\"> open " + inputFile + "</div><div align = \"center\"><b>" + projmodel->problemName + "</b></div><div><div align = \"left\"></div>");
	enableButtons();
	projmodel->debugBrowser = ui.textBrowser;
}
void slowwavedevice::readResults()
{
	bool att = sharedMemoryResult.isAttached();
	if (!att )
	{
		if (!sharedMemoryResult.attach()) return;
	}
	char *mem;

	int flag;
	int N;
	double *ar, *ai;

	sharedMemoryResult.lock();
	mem = (char*)sharedMemoryResult.constData();
	memcpy(&flag, mem, sizeof(int));
	if (!flag)
	{
		sharedMemoryResult.unlock();
		timer.start(500);
		return;
	}
	else
	{
		double tmpA;
		memcpy(&N, mem + sizeof(int), sizeof(int));
		if (!resultPlot->graph(0))
			createResultPlot(N);
		memcpy(resultStorageX, mem + 2 * sizeof(int), sizeof(double)*N);
		memcpy(resultStorageY, mem + 2 * sizeof(int) + sizeof(double)*N, sizeof(double)*N);
		*((int*)mem) = 0;
		sharedMemoryResult.unlock();
		updateResultPlot(N);
		timer.start(500);
		return;
	}
}
void slowwavedevice::createResultPlot(int Npoints)
{
	resultPlot->clearGraphs();
	resultPlot->addGraph();
	
	QPen pen(Qt::darkBlue);
	pen.setWidth(2);
	resultPlot->graph(0)->setPen(pen);
	showPlot(3);
}
void slowwavedevice::updateResultPlot(int Npoints)
{
	if (resultPlot->graphCount() == 0) return;
	if (resultPlot->graph(0))
		resData = resultPlot->graph(0)->data();
	else 
		return; 
	if (resData->size() != 0) resultPlot->graph(0)->clearData();
	
	
	for (int i = 0; i < Npoints; i++)
	{
		resData->insert(resultStorageX[i], QCPData(resultStorageX[i], resultStorageY[i]));
	}
	resultPlot->graph(0)->rescaleAxes();
	resultPlot->replot();
}
void slowwavedevice::enableButtons()
{
	ui.action->setEnabled(true);
	ui.treeView->setEnabled(true);
	ui.SaveAction->setEnabled(true);
	ui.SaveAsAction->setEnabled(true);
	ui.ExitAction->setEnabled(true);
	ui.OpenAction->setEnabled(true);
	ui.CalculateDispersionAction->setEnabled(true);
	ui.AbortAction->setEnabled(false);
}
void slowwavedevice::disableButtons()
{
	ui.treeView->setEnabled(false);
	ui.action->setEnabled(false);
	ui.SaveAction->setEnabled(false);
	ui.SaveAsAction->setEnabled(false);
	ui.ExitAction->setEnabled(false);
	ui.OpenAction->setEnabled(false);
	ui.CalculateDispersionAction->setEnabled(false);
	ui.AbortAction->setEnabled(true);
}
void slowwavedevice::showAbout()
{
	aboutDialog.show();
}