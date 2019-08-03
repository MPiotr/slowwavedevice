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

void writeData(const string & filename, const vector<DataRow> & data)
{
	ofstream file(filename);
	for (auto row : data) {
		file << row << "\n";
	}
}

void writeData(const string & filename, const vector<vector<DataRow >> & data)
{
	int counter = 0;
	for (const auto & mode : data) {
		ostringstream mode_file_name;
		mode_file_name << filename << counter << ".csv";
		writeData(mode_file_name.str(), mode);
		counter++;
	}
}