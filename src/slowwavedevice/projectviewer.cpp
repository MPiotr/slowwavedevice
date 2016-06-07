#include "projectviewer.h"
#include "standartxmlitem.h"
#include "projviewmenuaction.h"

projectviewer::projectviewer()
	: QStandardItemModel()
{
	
	}

projectviewer::~projectviewer()
{

}

void projectviewer::setprojContent(QFile* input)
{
	QString err_msg; int err_line; int err_column;
	doc.clear();
	clear();
	bool set_cont = doc.setContent(input, &err_msg, &err_line, &err_column);
	if (!set_cont)
	{
		emit(errorOpeningFile(qPrintable(err_msg), err_line, err_column));
		return;
	}

	QStandardItem *parentItem = invisibleRootItem();
	addNode(&doc.firstChild().firstChild().firstChild(), parentItem);
	
// -------- add cwd entry if necessary
	char dummy[300];
	if (setXMLEntry(&doc, "workingDirectory", dummy))
	{
		QDomNode currentItem = doc.elementsByTagName("workingDirectory").item(0);
		QDomNode parrent = currentItem.parentNode();
		parrent.removeChild(currentItem);
	}
		QDomElement newelement = doc.createElement("workingDirectory");
		QDomText newnodetext = doc.createTextNode(QDir::current().absolutePath());
		newelement.appendChild(newnodetext);
		doc.firstChild().firstChild().appendChild(newelement);

// ----------

}
void projectviewer::addOnlyNode(QDomNode *node, QStandardItem *parent, int insertrow = -1)
{
	QString elementName = node->nodeName();

	StandardXMLItem *  thisItemName = new StandardXMLItem(elementName, *node);
	thisItemName->setEditable(false);

	QString elementValue = node->toElement().text();
	StandardXMLItem *  thisItemValue = new StandardXMLItem(elementValue);
	thisItemValue->setEditable(true);

	if (elementName != QString("#comment"))
	{
		if (insertrow < 0)
			parent->appendRow(thisItemName);
		else
			parent->insertRow(insertrow, thisItemName);

		int par = 0;
		if (setXMLattribute(node, "iterable", &par)){
			StandardXMLItem* iter = new StandardXMLItem("iterable");
			QList<QStandardItem*> list;
			list.append(iter);
			parent->appendColumn(list);
		};

		if (!node->firstChild().isNull()){
			if (!node->firstChild().firstChild().isNull())
			{
				addNode(&(node->firstChild()), thisItemName);
			}
			else
			{
				thisItemName->appendRow(thisItemValue);
			}
		}
	}
	if (elementName == "problemName") problemName = node->toElement().text();

}
void projectviewer::addNode(QDomNode *node, QStandardItem *parent)
{
	QString elementName = node->nodeName();
	
	StandardXMLItem *  thisItemName = new StandardXMLItem(elementName, *node);
	thisItemName->setEditable(false);

	QString elementValue = node->toElement().text();
	StandardXMLItem *  thisItemValue = new StandardXMLItem(elementValue);
	thisItemValue->setEditable(true);
	
	if (elementName != QString("#comment"))
	{

		parent->appendRow(thisItemName);

		int par = 0;
		if (setXMLattribute(node, "iterable", &par)){
			StandardXMLItem* iter = new StandardXMLItem("iterable");
			QList<QStandardItem*> list;
			list.append(iter);
			parent->appendColumn(list);
		};

		if (!node->firstChild().isNull()){
			if (!node->firstChild().firstChild().isNull())
			{
				addNode(&(node->firstChild()), thisItemName);
			}
			else
			{
				thisItemName->appendRow(thisItemValue);
			}
		}
	}
	if (elementName == "problemName") problemName = node->toElement().text();
	
	if (!node->nextSibling().isNull())	addNode(&(node->nextSibling()), parent);
	
	
}
void projectviewer::itemChanged(QModelIndex ind)
{

	
}
void projectviewer::itemClicked(QModelIndex  index)
{
	StandardXMLItem *currentItem = (StandardXMLItem *)itemFromIndex(index);
	QDomElement current = currentItem->node.toElement();
	StandardXMLItem *parentItem = (StandardXMLItem *)itemFromIndex(parent(index));
	QString pname = currentItem->node.nodeName();

	if (parentItem && !parentItem->node.isNull())
	{
		QDomNode node = parentItem->node;
		QDomNode  parentparent = node.parentNode();
		QString name = node.nodeName();

		if (name == "dispersionFileName" || pname == "dispersionFileName" ||
			name == "voltage" || pname == "voltage" ) emit setVisiblePlot(0);
		if (name == "logitudinalStructureFile" || pname == "logitudinalStructureFile") emit setVisiblePlot(1);
		if (name == "fieldFileName" || pname == "fieldFileName") emit setVisiblePlot(2);
		
	}

}
void projectviewer::showContextMenu(QTreeView *view, QModelIndex index, QPoint point)
{
	
	StandardXMLItem *currentItem = (StandardXMLItem *)itemFromIndex(index);
	QDomElement current = currentItem->node.toElement();
	QString iterable = current.attribute("iterable");

	if (iterable == "0" || iterable == "1")
	{
		QMenu contextMenu(view);
		QString text;
		if (iterable == "0") text = "iterate value";
		else text = "single value";
		projviewMenuAction	action(text, &index, &contextMenu);
		contextMenu.addAction(&action);
		connect(&action, SIGNAL(triggered(QModelIndex*)), this, SLOT(toggleIterate(QModelIndex*)));
		contextMenu.exec(point);
	}

}
bool projectviewer::setData(const QModelIndex & index, const QVariant & value, int role)
{
	StandardXMLItem *currentItem = (StandardXMLItem *)itemFromIndex(index);
	StandardXMLItem *parentItem = (StandardXMLItem *)itemFromIndex(parent(index));
	StandardXMLItem *parentparentItem = (StandardXMLItem *)itemFromIndex(parent(parent(index)));
	
	if (!parentItem->node.isNull())
	{
		QDomNode node = parentItem->node;
		QDomNode  parentparent = node.parentNode();
		QString name = node.nodeName();
		
		QDomElement newelement = doc.createElement(name);
		QDomText newnodetext = doc.createTextNode(value.toString());
		newelement.appendChild(newnodetext);
		parentparent.replaceChild(newelement, node);

		parentItem->node = newelement;
		QDomNode tmpnode = newelement.parentNode();

		if (name == "voltage") emit voltageChanged(value.toString());

	}
	QString tmp = doc.toString();
	return QStandardItemModel::setData(index, value);
}

QDomDocument projectviewer::getDoc()
{
	return doc;
}

QString projectviewer::getXMLString()
{
	return doc.toString();
}

void projectviewer::setDispersionsPlot(QCustomPlot *plot, int Npoints)
{
	char dispersionFileName[200];
	double norm1;
	QDomNode LFsection = doc.elementsByTagName("LFsection").item(0);
	if (LFsection.isNull()){ printf("Error: LFsection entry is not found\n"); return; }
	if (!setXMLEntry(&LFsection, "period", &period)) { printf("Error: LF period is not found\n"); return; }
	
	Synchronism *syncwave;
	Interpolation *interpolation;

	double k1, fre;
	plot->clearGraphs();
	setXMLEntry(&doc, "frequency", &fre, &iteratedParams, &iteratedNames);	k1 = fre*2 * M_PI / 299.8;    	//волновое число (частота возбуждения)
	QVector<double> x(Npoints + 1);
	QVector<double> y(Npoints + 1);
	if (setXMLEntry(&doc, "dispersionFileName", dispersionFileName))
	{
		syncwave = new Synchronism(dispersionFileName, period, 299.8*k1 / (2.*M_PI));
		interpolation = new Interpolation(dispersionFileName);
		for (int i = 0; i <= Npoints; i++)
		{
			x[i] = 720 / (double)Npoints * i;
			y[i] = interpolation->at(x[i]);
		}
		plot->addGraph();
		plot->graph(0)->setData(x, y);
		plot->graph(0)->rescaleAxes();
		QPen pen(Qt::darkBlue);
		pen.setWidth(3);
		plot->graph(0)->setPen(pen);
	}
	else
	{
		plot->addGraph();
		plot->graph(0)->data()->insert(0, QCPData(0, fre));
		plot->graph(0)->data()->insert(720, QCPData(720, fre));
	}

	double voltage;
	setXMLEntry(&doc, "voltage", &voltage, &iteratedParams, &iteratedNames);
	double gamma = 1 + voltage / 511.;
	double coef = 299.8 / (2.*M_PI)*sqrt(gamma*gamma - 1) / gamma*M_PI / (180.*period);
	for (int i = 0; i <= Npoints; i++)
	{
		x[i] = 720. / (double) Npoints * i;
		y[i] = x[i] * coef;
	}
	plot->addGraph();
	plot->graph(1)->setData(x, y);
	plot->graph(1)->setPen(QPen(Qt::red));

	// ...............reading parasites data ........................
	int NumParasites = doc.elementsByTagName("parasite").length();
	int numgraph = 2;
	for (int i = 0; i < NumParasites; i++)
	{
		QDomNode parasite = doc.elementsByTagName("parasite").item(i);
		
		int chLength = parasite.childNodes().length();

		char dispFile[200];
		char fieldFile[200];
		if (!setXMLEntry(&parasite, "parasiteDispersion", (char*)dispFile))  continue;
		if (!setXMLEntry(&parasite, "parasiteField", fieldFile)) continue;

		syncwave->addMode(dispFile);

		Interpolation  parasite_disp(dispFile);
		for (int i = 0; i <= Npoints; i++)
		{
			x[i] = 720 / (double)Npoints * i;
			y[i] = parasite_disp.at(x[i]);
		}
		plot->addGraph();
		plot->graph(numgraph)->setData(x, y);
		numgraph++;
	}
	plot->rescaleAxes();
	// ........................................................

    


}
void projectviewer::setPlot(QCustomPlot *plot, char* entryName, int Npoints, QTextBrowser *browser)
{
	char fileName[200];
	QDomNode LFsection = doc.elementsByTagName("LFsection").item(0);
	if (LFsection.isNull()){	browser->append("<b>LFsection entry is not found</b>"); return;}
	if (!setXMLEntry(&LFsection, "period", &period)) { browser->append("<b>LFsection period is not found</b>"); return; }
	plot->clearGraphs();
	
	Interpolation *interpolation;

	double k1;

	setXMLEntry(&doc, "frequency", &k1, &iteratedParams, &iteratedNames);	k1 *= 2 * M_PI / 299.8;    	//волновое число (частота возбуждения)
	QVector<double> x(Npoints + 2);
	QVector<double> y(Npoints + 2);
	if (setXMLEntry(&doc, entryName, fileName))
	{
		interpolation = new Interpolation(fileName);
		
		double max = interpolation->xMax();
		double min = interpolation->xMin();
		for (int i = 0; i <= Npoints; i++)
		{
			double cx = min + (max - min) / (double)Npoints * i;
			x.append(cx);
			y.append(interpolation->at(cx));
		}
		
		plot->addGraph();
		plot->graph(0)->setData(x, y);
		plot->graph(0)->rescaleAxes();
		QPen pen(Qt::darkRed);
		pen.setWidth(2);
		plot->graph(0)->setPen(pen);
		plot->rescaleAxes();
	}
	else return;
}
void projectviewer::toggleIterate(QModelIndex* index)
{
	StandardXMLItem *currentItem = (StandardXMLItem *)itemFromIndex(*index);
	QDomElement current = currentItem->node.toElement();
	StandardXMLItem *parentItem = (StandardXMLItem *)itemFromIndex(parent(*index));
	QString iterable = current.attribute("iterable");

	QString tmp1 = doc.toString();

	if (!parentItem->node.isNull())
	{
		QDomNode node = current;
		QDomNode  parent = node.parentNode();
		QString name = node.nodeName();

		if (iterable == "0")
		{
			QString fromText = current.text();
			QDomNode newFrom = doc.createElement("from");
			newFrom.appendChild(doc.createTextNode(fromText));
			QDomNode newTo = doc.createElement("to");
			newTo.appendChild(doc.createTextNode(fromText));
			QDomNode newStep = doc.createElement("step");
			newStep.appendChild(doc.createTextNode("1"));

			QDomElement newelement = doc.createElement(name);
			newelement.setAttribute("iterable", "1");
			newelement.appendChild(newFrom);
			newelement.appendChild(newTo);
			newelement.appendChild(newStep);
			parent.replaceChild(newelement, node);

			int row = currentItem->row();
			parentItem->removeRow(row);
			addOnlyNode(&newelement, parentItem, row);
		}
		else
		{

			QString fromtext = current.childNodes().at(0).toElement().text();
			QDomText newnodetext = doc.createTextNode(fromtext);

			QDomElement newelement = doc.createElement(name);
			newelement.setAttribute("iterable", "0");
			newelement.appendChild(newnodetext);
			parent.replaceChild(newelement, node);

			int row = currentItem->row();
			parentItem->removeRow(row);
			addOnlyNode(&newelement, parentItem, row);
		}


	}
	QString tmp = doc.toString();


}
