#include "standartxmlitem.h"

StandardXMLItem::StandardXMLItem()
	: QStandardItem()
{

}
StandardXMLItem::StandardXMLItem(const QString &text)
	: QStandardItem(text)
{
	
}
StandardXMLItem::StandardXMLItem(const QString &text, QDomNode _node)
	: QStandardItem(text)
{
	node = _node;
	nodeName = node.nodeName();
	if (node.toElement().attribute("type") == "filename") type = "filename";
}

StandardXMLItem::~StandardXMLItem()
{

}
