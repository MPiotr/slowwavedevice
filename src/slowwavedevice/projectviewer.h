#ifndef PROJECTVIEWER_H
#define PROJECTVIEWER_H

#include <QTreeView>
#include <QDomDocument>
#include <qstandarditemmodel.h>
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
	void setprojContent(QFile *inputfile);
	~projectviewer();
	bool setData(const QModelIndex & index, const QVariant & value, int role);
	QDomDocument getDoc();
	QString getXMLString();
	QString problemName;
	void setDispersionsPlot(QCustomPlot *plot, QTextEdit* console, int Npoints);
	void setPlot(QCustomPlot *plot, char *entryname,  int Npoints, QTextBrowser *);
	double period;	
	QTextBrowser *debugBrowser;
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
	QDomDocument toolTipsXML;
	void addNode(QDomNode *, QStandardItem  *);
	void addOnlyNode(QDomNode *, QStandardItem  *, int );

	vector< vec > iteratedParams;
	vector< string > iteratedNames;
	QIcon checkedIcon;
	QIcon uncheckedIcon;

	
};

#endif // PROJECTVIEWER_H
