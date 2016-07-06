#ifndef PROJVIEWMENUACTION_H
#define PROJVIEWMENUACTION_H

#include <QAction>

class projviewMenuAction : public QAction
{
	Q_OBJECT

public:
	projviewMenuAction(QString text, QModelIndex *index, QObject *parent);
	~projviewMenuAction();

signals:
	void triggered(QModelIndex *);
public slots:
	void itemTriggered();
private:
	QModelIndex *index;
	
};

#endif // PROJVIEWMENUACTION_H
