#include "dataIO.h"
#include <fstream>


vector<DataRow> readData(const string & filename, int key_column, int value_column)
{
	vector<DataRow> data;
	ifstream input(filename);
	if (!input.good()) {
		cerr << "Can't open " << filename << "; " << __FILE__ << ":" << __LINE__ << endl;
		throw  exception::exception((string("Can't open" + filename + "; ")).data());
	}
	string line;
	while (getline(input, line)) {
		istringstream sstream(line);
		double x;
		char delim;
		DataRow row(key_column, value_column);
		while (!sstream.eof()) {
			sstream >> x >> delim;
			row.push_back(x);
		}
		data.push_back(row);
	}

	sort(data.begin(), data.end());

	return data;
}

void writeData(const string & filename, const vector<DataRow> & data, int save_key, int save_value, double key_transform, double value_transform)
{
	ofstream file(filename);
	for (auto row : data) {
		if (save_key != -1 && save_value != -1 && save_key < row.size() && save_value < row.size())
			file << key_transform*row[save_key] << "," << value_transform*row[save_value] << "\n";
		else 
			file << row << "\n";
	}
}

void writeData(const string & filename, const vector<vector<DataRow >> & data, int save_key, int save_value, double key_transform, double value_transform)
{
	int counter = 0;
	for (const auto & mode : data) {
		ostringstream mode_file_name;
		mode_file_name << filename << counter << ".csv";
		cout << "writing " << mode_file_name.str() << "\n";
		writeData(mode_file_name.str(), mode, save_key, save_value, key_transform, value_transform);
		counter++;
	}
}

void rescaleColumn(vector<DataRow> &data, int column, double factor)
{
	for (DataRow & row : data) {
		row[column] *= factor;
	}
}

void rescaleColumn(vector<vector<DataRow> > &data, int column, double factor)
{ 
	for (auto & mode : data) {
		rescaleColumn(mode, column, factor);
	}
}