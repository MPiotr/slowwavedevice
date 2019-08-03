// filemerger.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <algorithm>
#include <filesystem>
#include <iterator>
#include <utility>
#include "dataIO.h"

vector<DataRow> merge(const vector<string> & inputFiles, int key_col, int val_col) {
	
	vector<DataRow> result;
	for (const string & file : inputFiles) 
	{
		auto data = readData(file, key_col, val_col);
		if (result.empty()) result = move(data);
		else result.insert(end(result), make_move_iterator(begin(data)), make_move_iterator(end(data)));
		//move(begin(data), end(data), begin(result));		
	}

	sort(begin(result), end(result));

	return result;
}

vector <string> findFilesFromPrefix(const string &prefix) {
	vector<string> result;
	string name = prefix + "0.csv";
	int counter = 0;
	while (std::experimental::filesystem::exists(name)) {
		result.push_back(name);		
		counter++;
		name = prefix + to_string(counter) + ".csv";
	}
	return result;
}

int main(int argc, char** argv)
{
	if (argc < 4) {
		cout << "Usage: filemerger.exe [input_prefix for csv-file] [key column] [value column]" << endl;
		return 1;
	}

	string prefix = argv[1];
	int key_col = stoi(argv[2]);
	int val_col = stoi(argv[3]);

	vector<string> inputFiles = findFilesFromPrefix(prefix);
	auto result = merge(inputFiles, key_col, val_col);
	writeData(prefix + ".csv", result);

}
