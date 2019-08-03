#include "datarow.h"

bool operator< (const DataRow & x, const DataRow &y)
{
	if (x.key() == y.key())
		return x.val() < y.val();
	else
		return x.key() < y.key();
}

ostream &  operator << (ostream &str, const DataRow & row)
{
	str << row.serialize();
	return str;
}