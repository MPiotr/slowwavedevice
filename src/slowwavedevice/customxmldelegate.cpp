#include "customxmldelegate.h"
#include "standartxmlitem.h"



QWidget* customxmldelegate::createEditor(QWidget *parent, const QStyleOptionViewItem &option, const QModelIndex &index) const
{
	StandardXMLItem *currentItem = (StandardXMLItem *) theModel->itemFromIndex(index.parent());
	QString name = currentItem->node.nodeName();
	QString value = currentItem->node.toElement().text();
	QDomElement helperElement = helperXML->elementsByTagName(name).at(0).toElement();
	QString  attribute = helperElement.attribute("type");
	if (helperElement.attribute("type") == QString("fixedOptions"))
	{
		QDomNodeList childs = helperElement.childNodes();
		QComboBox *editor = new QComboBox(parent);
		editor->setEditable(false);
		for (int i = 0; i < childs.count(); i++)
		{
			QDomNode child = childs.at(i);
			if (child.nodeName() == "option")
			{
				QString option = child.toElement().attribute("value");
				editor->addItem(option);
			}
		}
		int selected = editor->findText(value);
		editor->setCurrentIndex(selected);
		return editor;
	}
	
	
	return	QItemDelegate::createEditor(parent, option, index);



//	StandardXMLItem *currentItem = (StandardXMLItem *) QStandardItemModel.itemFromIndex(index);
//    return	QItemDelegate::createEditor(parent, option, index);
}



