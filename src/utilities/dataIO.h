#pragma once
#include  "datarow.h"

using namespace std;

vector<DataRow> readData(const string & filename, int key_column, int value_column);
void writeData(const string & filename, const vector<DataRow> & data);
void writeData(const string & filename, const vector<vector<DataRow >> & data);