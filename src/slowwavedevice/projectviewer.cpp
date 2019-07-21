#include <fstream>
#include <vector>
#include <algorithm>
#include "projectviewer.h"
#include "standartxmlitem.h"
#include "projviewmenuaction.h"

projectviewer::projectviewer()
	: QStandardItemModel()
{
	QString cwd = QDir::current().absolutePath();
	QFile input("helpInfo.xml");
	if (!input.open(QIODevice::ReadOnly)){
		return;
	}
	QString err_msg; int err_line; int err_column;
	bool readhelp =  toolTipsXML.setContent(&input, &err_msg, &err_line, &err_column);
	if (!readhelp)
	{
		emit(errorOpeningFile(qPrintable(err_msg), err_line, err_column));
		return;
	}
	input.close();
}

projectviewer::~projectviewer()
{

}

bool    projectviewer::openHelperXML(QDomDocument* helper, QString *err_msg,  int *err_line, int *err_column){
	QFile input("helpInfo.xml");
	if (!input.open(QIODevice::ReadOnly)){
		*err_msg = "File helpInfo.xml is not available";
		*err_line = 0;
		*err_column = 0;
		return false;
	}

	bool readhelp = helper->setContent(&input, err_msg, err_line, err_column);
	input.close();
	return readhelp;
}

void projectviewer::setprojContent(QFile* input, QCustomPlot* _plotArea)
{
	QString err_msg; int err_line; int err_column;
	plotArea = _plotArea;
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
void projectviewer::addOnlyNode(QDomNode *node, QStandardItem *parent, int insertrow = -1) //Adds only one node without recurssion. Used for toggling iterate
{
	QString elementName = node->nodeName();
	QDomNode helperNode = toolTipsXML.elementsByTagName(elementName).at(0);

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
		nameToIndex[elementName] = thisItemValue->index();
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

		QDomNodeList tooltips = toolTipsXML.elementsByTagName(thisItemName->nodeName);

		if (!tooltips.isEmpty())
		{
			QString  nodevalue = tooltips.item(0).toElement().firstChild().nodeValue();
			thisItemName->setToolTip(nodevalue);
		}
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
				nameToIndex[elementName] = thisItemValue->index();
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
	QString name = currentItem->node.nodeName();

	if (parentItem && !parentItem->node.isNull())
	{
		QDomNode node = parentItem->node;
		QDomNode  parentparent = node.parentNode();
		QString parname = node.nodeName();

		if (parname == "dispersionFileName" || name == "dispersionFileName" ||
			parname == "voltage" || name == "voltage") emit setVisiblePlot(0);
		if (parname == "longitudinalStructureFile"	|| parname == "fieldFileName" || parname == "QStructureFile") 
		{
			setPlot(plotArea, parname.toLocal8Bit().data(), 50, nullptr);
			emit setVisiblePlot(1);
		}
		if (name == "longitudinalStructureFile"||  name == "fieldFileName"|| name == "QStructureFile")
		{
			setPlot(plotArea, name.toLocal8Bit().data(), 50, nullptr);
			emit setVisiblePlot(1);
		}
		if (name == "periodShape") {
			setTablePlot(plotArea, name.toLocal8Bit().data(), 50, nullptr);
			emit setVisiblePlot(1);
		}
		if (parname == "periodShape") {
			setTablePlot(plotArea, parname.toLocal8Bit().data(), 50, nullptr);
			emit setVisiblePlot(1);
		}
		
	}

}
void projectviewer::showContextMenu(QTreeView *view, QModelIndex index, QPoint point)
{
	
	StandardXMLItem *currentItem = (StandardXMLItem *)itemFromIndex(index);
	QString elementName = currentItem->nodeName;
	QDomNode helperNode = toolTipsXML.elementsByTagName(elementName).at(0);

	QDomElement current = helperNode.toElement(); //  currentItem->node.toElement();
	QString iterable = current.attribute("iterable");

	QDomNode itemNode = currentItem->node;
	QString iterableItem = itemNode.toElement().attribute("iterable");

	if (iterable == "0" || iterable == "1")
	{
		QMenu contextMenu(view);
		QString text;
		if (iterableItem == "0") text = "iterate value";
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

void projectviewer::setDispersionsPlot(QCustomPlot *plot, QTextEdit* console, int Npoints)
{
	char dispersionFileName[200] = "";
	double norm1;
	QDomNode LFsection = doc.elementsByTagName("LFsection").item(0);
	if (LFsection.isNull()){
		console->append("<b><font color=\"red\">Error: LFsection entry is not found</font></b>"); return;
	}
	if (!setXMLEntry(&LFsection, "period", &period)) {
		console->append("<b><font color=\"red\">Error: LF period is not found</font></b>");  return;
	}
	
	Synchronism *syncwave;
	Interpolation *interpolation;

	double k1, fre;
	plot->clearPlottables();
	setXMLEntry(&doc, "frequency", &fre, &iteratedParams, &iteratedNames);	k1 = fre*2 * M_PI / 299.8;    	//волновое число (частота возбуждения)
	QVector<double> x(Npoints + 1);
	QVector<double> y(Npoints + 1);
	if (setXMLEntry(&doc, "dispersionFileName", dispersionFileName) && QFile(QString(dispersionFileName)).exists())
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
		if (!setXMLEntry(&doc, "dispersionFileName", dispersionFileName)) console->append("<b><font color=\"red\">Error</font></b> reading file <b>"+QString(dispersionFileName)+"</b>");
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
	//plot->rescaleAxes();
	// ........................................................

    


}
void projectviewer::setPlot(QCustomPlot *plot, char* entryName, int Npoints, QTextBrowser *browser)
{
	char fileName[200];
	QDomNode LFsection = doc.elementsByTagName("LFsection").item(0);
	if (LFsection.isNull() 
		&& browser != nullptr){	browser->append("<b>LFsection entry is not found</b>"); return;}
	if (!setXMLEntry(&LFsection, "period", &period)
		&& browser != nullptr){ browser->append("<b>LFsection period is not found</b>"); return; }
	plot->clearPlottables();
	
	Interpolation *interpolation = plotInterpol;

	double k1;
	 
	setXMLEntry(&doc, "frequency", &k1, &iteratedParams, &iteratedNames);	k1 *= 2 * M_PI / 299.8;    	//волновое число (частота возбуждения)
	QVector<double> x;
	QVector<double> y;
	if (setXMLEntry(&doc, entryName, fileName))
	{
		if (!QFile(QString(fileName)).exists()) {
			browser->append("<b><font color = \"red\">Error</font></b> reading file <b>" + QString(fileName) + "</b>"); 
			return;
		}
		if (interpolation == nullptr)
			interpolation = new Interpolation(fileName);
		else
			interpolation->reload(fileName);
	
		
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
void projectviewer::setTablePlot(QCustomPlot *plot, char* entryName, int Npoints, QTextBrowser *browser)
{
	char fileName[200];
	plot->clearPlottables();
	QVector<double> x;
	QVector<double> y;
	if (setXMLEntry(&doc, entryName, fileName))
	{	
		FILE *file = fopen(fileName, "r");			
		if (file == nullptr) {
			browser->append("<b><font color = \"red\">Error</font></b> opening file <b>" + QString(fileName) + "</b>");
			return;
		}
		//std::ifstream input(file);
		float _x, _y;
		while (fscanf(file, "%g,%g\n", &_x, &_y) == 2)
		{
			//input >> _x >> _y;
			x.push_back(_x);
			y.push_back(_y);
		}
		fclose(file);
		

		plot->addGraph();
		QCPCurve* curve = new QCPCurve(plot->graph(0)->keyAxis(), plot->graph(0)->valueAxis());
		QPen pen(Qt::darkGreen);
		pen.setWidth(2);
		curve->setData(x, y);
		curve->setPen(pen);
		curve->rescaleAxes();
		plot->addPlottable(curve);
		//plot->graph(0)->rescaleAxes();
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

		if (iterable == "0" || iterable == NULL)
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

void projectviewer::recalculatePeriodFromShape(QTextBrowser *browser)
{
	char fileName[200];
	if (!setXMLEntry(&doc, "periodShape", fileName)) {
		browser->append("<b><font color = \"red\">Error</font></b> XML entry <b>periodShape</b> not found");
		return;
	}
	FILE *file = fopen(fileName, "r");
	if (file == nullptr) {
		browser->append("<b><font color = \"red\">Error</font></b> opening file <b>" + QString(fileName) + "</b>");
		return;
	}
	//std::ifstream input(file);
	std::vector<double> x;
	std::vector<double> y;
	float _x, _y;
	while (fscanf(file, "%g,%g\n", &_x, &_y) == 2)
	{
		//input >> _x >> _y;
		x.push_back(_x);
		y.push_back(_y);
	}
	fclose(file);

	double xmax = *std::max_element(begin(x), end(x));
	double xmin = *std::min_element(begin(x), end(x));
	
	double period = xmax - xmin;

	auto el = doc.elementsByTagName("period").at(0);
	el.replaceChild(doc.createTextNode(QString::number(period)), el.firstChild());

	auto periodIndex = nameToIndex["period"];
	auto periodItem = itemFromIndex(periodIndex);

	if (periodItem == nullptr) {
		browser->append("<b><font color = \"red\">Error</font></b> can't find 'period' item");
		return;
	}
	periodItem->setText(QString::number(period));	
	
}

