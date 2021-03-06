#ifndef PROJECTVIEWER_H
#define PROJECTVIEWER_H

#include <QTreeView>
#include <QDomDocument>
#include <qstandarditemmodel.h>
#include <map>
#include <qfile.h>
#include <qabstractxmlreceiver.h>
#include "qcustomplot.h"
#include "synchronism.h"
#include "xml_routines.h"

class projectviewer : public QStandardItemModel
{
	Q_OBJECT

public:
	projectviewer();
	void setprojContent(QFile *inputfile, QCustomPlot* plotArea);
	~projectviewer();
	bool setData(const QModelIndex & index, const QVariant & value, int role);
	QDomDocument getDoc();
	QString getXMLString();
	QString problemName;
	void setDispersionsPlot(QCustomPlot *plot, QTextEdit* console, int Npoints);
	void setPlot(QCustomPlot *plot, char *entryname,  int Npoints, QTextBrowser *);
	void setTablePlot(QCustomPlot *plot, char* entryName, int Npoints, QTextBrowser *browser);
	double period;	
	QTextBrowser *debugBrowser;
	QDomDocument toolTipsXML;
	void recalculatePeriodFromShape(QTextBrowser *browser);

	bool savePeriodParamsForDispersionCalculation(QFile* file);
	bool savePeriodParamsForFieldCalculation(QFile * file);
	bool getShapeFileName(QString *name);
	int  getNumCores();
	double getPeriod();
	double getRefFreq();
	QString getDispersionFileName();
	QString getFieldFileName();
	QString getShapeType();
	bool getLtransversalOrM(double &Ltransversal);


	static bool openHelperXML(QDomDocument *out, QString *err_msg, int *err_line, int *err_column);
public slots:
	void itemChanged(QModelIndex);
	void itemClicked(QModelIndex);
	void showContextMenu(QTreeView *view, QModelIndex index, QPoint point);
	void toggleIterate(QModelIndex *);
signals: 
	void errorOpeningFile(QString line, int row, int column);
	void voltageChanged(QString newvalue);
	void setVisiblePlot(int plotNum);

private:
	QDomDocument doc; 

	QCustomPlot *plotArea = nullptr;
	Interpolation *plotInterpol = nullptr;
	
	void addNode(QDomNode *, QStandardItem  *);
	void addOnlyNode(QDomNode *, QStandardItem  *, int );

	vector< vec > iteratedParams;
	vector< string > iteratedNames;

	map<QString, QModelIndex> nameToIndex;

	QIcon checkedIcon;
	QIcon uncheckedIcon;

	
};

#endif // PROJECTVIEWER_H
