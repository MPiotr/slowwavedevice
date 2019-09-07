#pragma once
#include  "datarow.h"

using namespace std;

vector<DataRow> readData(const string & filename, int key_column, int value_column);
void writeData(const string & filename, const vector<DataRow> & data, int write_key_column = -1, int write_value_column = -1, double key_transform = 1, double value_transform = 1);
void writeData(const string & filename, const vector<vector<DataRow >> & data, int write_key_column = -1, int write_value_column = -1, double key_transform = 1, double value_transform = 1);

void rescaleColumn(vector<DataRow> &data, int column, double factor);
void rescaleColumn(vector<vector<DataRow> > &data, int column, double factor);