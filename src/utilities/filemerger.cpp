// filemerger.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <algorithm>
#include <filesystem>
#include <iterator>
#include <utility>
#include <io.h>
#include "dataIO.h"

vector<DataRow> merge(const vector<string> & inputFiles, int key_col, int val_col) {
	
	vector<DataRow> result;
	cout << "Merging files: ";
	for (auto name : inputFiles)    cout << name << ", ";
	cout << "\n\n";


	for (const string & file : inputFiles) 
	{
		cout << file << "...";
		auto data = readData(file, key_col, val_col);
		cout << " read ...";
		if (result.empty()) result = move(data);
		else result.insert(end(result), make_move_iterator(begin(data)), make_move_iterator(end(data)));
		//move(begin(data), end(data), begin(result));	
		cout << "OK\n";
	}

	sort(begin(result), end(result));

	return result;
}

vector <string> findFilesFromPrefix(const string &prefix) {
	vector<string> result;
	string name = prefix + "0.csv";
	int counter = 0;	
//	while (std::experimental::filesystem::exists(name)) {
	while (_access(name.data(), 0) == 0) {
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
