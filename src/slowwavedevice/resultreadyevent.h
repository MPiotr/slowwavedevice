#ifndef RESULTREADYEVENT_H
#define RESULTREADYEVENT_H

#include <QEvent>

class resultReadyEvent : public QEvent
{
public:
	static const QEvent::Type myType = static_cast<QEvent::Type>(2000);
	resultReadyEvent(): QEvent((QEvent::Type)2000){}
	~resultReadyEvent(){};
	
	
};

#endif // RESULTREADYEVENT_H
