#ifndef STANDARTXMLITEM_H
#define STANDARTXMLITEM_H

#include <QObject>
#include <qstandarditemmodel.h>
#include <QDomDocument>


class StandardXMLItem : public QStandardItem
{
//	Q_OBJECT

public:
	StandardXMLItem();
	StandardXMLItem(const QString &text);
	StandardXMLItem(const QString &text, QDomNode node);
	~StandardXMLItem();
	
	QDomNode node;
	QString nodeName;
	QString type;

private:
		
};
Q_DECLARE_METATYPE(StandardXMLItem)
#endif // STANDARTXMLITEM_H
