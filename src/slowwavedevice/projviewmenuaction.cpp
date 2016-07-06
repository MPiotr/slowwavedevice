#include "projviewmenuaction.h"

projviewMenuAction::projviewMenuAction(QString text, QModelIndex *ind, QObject *parent)
: QAction(text, parent)
{
	index = ind;
	connect(this, SIGNAL(triggered()), this, SLOT(itemTriggered()));
}

projviewMenuAction::~projviewMenuAction()
{

}
void projviewMenuAction::itemTriggered()
{
	emit triggered(index);
}