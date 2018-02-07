#ifndef CUSTOMXMLDELEGATE_H
#define CUSTOMXMLDELEGATE_H

#include <QItemDelegate>
#include "projectviewer.h"

class customxmldelegate : public QItemDelegate
{
	Q_OBJECT

public:
	customxmldelegate(projectviewer * model, QDomDocument *helper, QObject *parent = 0) : QItemDelegate(parent), theModel(model), helperXML(helper) {};

	QWidget *createEditor(QWidget *parent, const QStyleOptionViewItem &option, const QModelIndex &index) const;
//	void setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const;
//	void paint(QPainter *painter, const QStyleOptionViewItem &option, const QModelIndex &index) const;

private:
	projectviewer *theModel;	
	QDomDocument *helperXML;
	
};

#endif // CUSTOMXMLDELEGATE_H
