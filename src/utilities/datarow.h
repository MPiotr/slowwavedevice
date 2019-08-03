#pragma once
#include<vector>
#include<algorithm>
#include <iostream>
#include <sstream>

/*
   Simple wrap around vector<double> where two indices, key and value, have special meaning.
*/

using namespace std;

class DataRow {
public:
	DataRow(vector<double> & _data, int key_column, int value_column) : key_col(key_column), val_col(value_column), data(move(_data)) {
		if (data.size() <= max(key_col, val_col)) {
			cerr << "Max of (key_col=" << key_col << ", val_col=" << val_col << ") is greater than size of the data = " << data.size() << endl;
			throw exception("Wrong input");
		}
	}

	DataRow(int key_column, int value_column) : key_col(key_column), val_col(value_column) { ; }
	DataRow(const DataRow & other) :key_col(other.key_col), val_col(other.val_col), data(other.data) { ; }



	inline void push_back(double v) { data.push_back(v); }
	inline size_t size()      const { return data.size(); }
	inline int key_column()   const { return key_col; }
	inline int value_column() const { return val_col; }
	inline double key()       const { return data[key_col]; }
	inline double val()       const { return data[val_col]; }

	inline string serialize() const {
		ostringstream str;
		for (int i = 0; i < data.size(); i++) {
			double value = data[i];
			str << value;
			if (i != data.size() - 1) str << ",";
		}
		return str.str();
	}

private:
	int key_col, val_col;
	vector<double> data;
};

bool operator< (const DataRow & x, const DataRow &y);
ostream &  operator << (ostream &str, const DataRow & row);