#include "slowwavedevice.h"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	slowwavedevice w;
	w.show();
	return a.exec();
}
