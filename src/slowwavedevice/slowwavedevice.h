#ifndef SLOWWAVEDEVICE_H
#define SLOWWAVEDEVICE_H

#include <QtWidgets/QMainWindow>
#include <QAbstractItemModel>
#include <qsharedmemory.h>
#include "qcustomplot.h"
#include <qprocess.h>
#include <qbuffer.h>
#include "projectviewer.h"

#include "ui_slowwavedevice.h"
#include "ui_aboutDialog.h"

class slowwavedevice : public QMainWindow
{
	Q_OBJECT

public:
	slowwavedevice(QWidget *parent = 0);
	~slowwavedevice();

private:
	Ui::slowwavedeviceClass ui;
	Ui::Dialog about;

	QMenu *fileMenu;
	QTimer timer;
	
	QAction *openFile;
	QAction *recentFileActs[5];
	
	QString inputFile;

	QCustomPlot *qplot;
	QCustomPlot *couplingPlot;
	QCustomPlot *resultPlot;
	QCPDataMap *resData;
	double *resultStorageX;
	double *resultStorageY;

	projectviewer *projmodel;
	QSharedMemory sharedMemory;
	QSharedMemory sharedMemoryResult;

	QProcess solverProcess;
	QProcess dispersionProcess;
	QToolBar *toolbar;
	QDialog aboutDialog;
	
	QVector<double> *plotx, *ploty;
	QString lastDirectory;
	QStack<QString> recentFiles;
	QFile cookie;

	void loadRecentFiles();
	void saveRecentFiles();
	void updateRecentFilesActions();
	void openFileFun(QString);
	void enableButtons();
	void disableButtons();

	double dXResult();
	QBuffer resultBuffer;
	QDataStream resultStream;
	void createResultPlot(int Npoints);
	void updateResultPlot(int Npoints);

signals:
	void projmodelContextMenu(QTreeView *, QModelIndex, QPoint);
public	slots: 
	void replotVoltage(QString voltage);
	void showPlot(int n);
	void requestModelContextMenu(QPoint);
	void openRecentFile();
	void saveFile();
	void saveFileAs();
	void exit();
	void readResults();
private slots:
	void loadFile();
	void start();
	void solverFinished(int sig, QProcess::ExitStatus stat);
	void abortSolver();
	void readConsole();
	void showAbout();
	void calculateDispersion();
};

#endif // SLOWWAVEDEVICE_H
