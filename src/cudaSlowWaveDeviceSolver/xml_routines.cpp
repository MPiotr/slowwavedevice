#include "xml_routines.h"

bool setXMLEntry(QDomNode* node, char *name, double *par)
{
	QDomNode currentItem = node->namedItem(name);
	if (currentItem.isNull()) return false;
	*par = currentItem.toElement().text().toDouble();
	return true;
}
bool setXMLEntry(QDomNode* node, char *name, int *par)
{
	QDomNode currentItem = node->namedItem(name);
	if (currentItem.isNull()) return false;
	*par = currentItem.toElement().text().toInt();
	return true;
}
bool setXMLEntry(QDomNode* node, char *name, char *par)
{
	QDomNode currentItem = node->namedItem(name);
	if (currentItem.isNull()) return false;
	
	QString text = currentItem.toElement().text();
	QByteArray btmp = text.toLocal8Bit();

	memcpy(par, (void*)btmp.data(), sizeof(char)*text.length());
	par[text.length()] = '\0';
	return true;
}
bool setXMLEntry(QDomDocument* doc, char *name, char *par)
{
	QDomNode currentItem = doc->elementsByTagName(name).item(0);
	if (currentItem.isNull()) return false;
	QString text = currentItem.toElement().text();
	QByteArray btmp = text.toLocal8Bit();

	memcpy(par, (void*)btmp.data(), sizeof(char)*text.length());
	par[text.length()] = '\0';
	return true;
}
bool setXMLEntry(QDomDocument* doc, char *name, double *par, vector<vec> *iterPar, vector<string> *iterNames)
{
	QDomNode currentItem = doc->elementsByTagName(name).item(0);
	if (currentItem.isNull()) return false;
	if (currentItem.childNodes().length() > 1)
	{
		double from, to, step;
		setXMLEntry(&currentItem, "from", &from);
		setXMLEntry(&currentItem, "to", &to);
		if (!setXMLEntry(&currentItem, "step", &step)) step = 0;
		*par = from;
		if (iterPar != NULL)
		{
			iterNames->push_back(name);

			if (step == 0) step = to - from;
			if (step*(to - from) < 0) step *= -1;
			vec values;

			for (double value = from; (step > 0) ? value <= to : value >= to; value += step)
				values.push_back(value);

			iterPar->push_back(values);
		}
	}
	else
	{
		*par = currentItem.toElement().text().toDouble();
	}
	return true;
}
bool setXMLEntry(QDomDocument* doc, char *name, int *par)
{
	QDomNode currentItem = doc->elementsByTagName(name).item(0);
	if (currentItem.isNull()) return false;
	*par = currentItem.toElement().text().toInt();
	return true;
}
bool setXMLattribute(QDomDocument* doc, char *entryname, char* atributename, int *par)
{
	QDomNode currentItem = doc->elementsByTagName(entryname).item(0);
	if (currentItem.isNull()) return false;
	QDomNode atributeItem = currentItem.attributes().namedItem(atributename);
	if (atributeItem.isNull()) return false;
	*par = atributeItem.nodeValue().toInt();
	return true;
}
bool setXMLattribute(QDomNode*  node, char* atributename, int *par)
{
	if (node->isNull()) return false;
	QDomNode atributeItem = node->attributes().namedItem(atributename);
	if (atributeItem.isNull()) return false;
	*par = atributeItem.nodeValue().toInt();
	return true;
}