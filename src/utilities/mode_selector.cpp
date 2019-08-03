#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <exception>
#include <algorithm>
#include "datarow.h"
#include "dataIO.h"

using namespace std;

vector< DataRow >::iterator selectNext(vector< DataRow >::iterator & nose, vector< DataRow>::iterator & tail, 
	                                   double & y, double  &vy, const double alpha) 
{
	double predicted_y = y + vy;
	auto iter = min_element(nose, tail, 
		       [predicted_y] (DataRow a, DataRow b)
	           {
	       	    return fabs(a.val() - predicted_y) < fabs(b.val() - predicted_y);
	           });

	double nexty = (*iter).val();
	vy += alpha * ((nexty - y) - vy);
	y = nexty;
	
	return iter;
}

void selectFirst(vector<DataRow >::iterator & nose, vector<DataRow >::iterator & tail,
	             const  vector<DataRow>::iterator & dataend,
	             vector<vector<DataRow > >& modes,
	                    vector<double> & y, vector<double> & vy,
	                    int num_modes) 
{

	while ((*nose).key() == (*tail).key() && tail != dataend )
		tail++;

	for (int i = 0; i < num_modes; i++)
	{
		if (i < tail - nose) {
			y[i] = (*(nose + i)).val();
			modes[i].push_back(*(nose+i));
		}
		else {
			y[i] = 0;	
			vector<double> zeroVec = vector<double>(nose->size(), 0);
			DataRow zeroRow(zeroVec, 0, 0);
			modes[i].push_back(zeroRow);
		}
		vy[i] = 0;		
	}

}

void moveIterators(vector<DataRow>::iterator & nose,  vector<DataRow >::iterator & tail,
	               const vector<DataRow>&data)
{
	nose = tail;
	if (tail == data.end()) return;
	double newkey = (*tail).key();
	while (tail != data.end() && (*tail).size() != 0 && (*tail).key() == newkey)
		tail++;	
}

vector<vector<DataRow> > propagate(vector<DataRow> &data, int num_modes, double alpha) {
	vector<vector<DataRow> > modes(num_modes);
	vector<DataRow>::iterator nose, tail;

	nose = data.begin();
	tail = data.begin();
	vector<double> vy(num_modes);
	vector<double>  y(num_modes);

	selectFirst(nose, tail, data.end(), modes, y, vy, num_modes);

	while(true) {
		vector<double> next_data(num_modes);
		moveIterators(nose, tail, data);
		if (nose == tail) break;
		for (int j = 0; j < num_modes; j++)
		{			
			auto next_iter = selectNext(nose, tail, y[j], vy[j], alpha);		
			modes[j].push_back(*next_iter);
		}		
	} 

	return modes;
}

int main(int argc, char** argv)
{
	cout << "arc = " << argc << endl;
	if (argc < 6) {
		cout << "Usage: modeselector [input - csv-file] [key_column - zero-based] [value_column - zero-based] [num_modes] [output prefix] <optionally: alpha [def 0.8]>" << endl;
		return 1;
	}
	string input_file = argv[1];
	string prefix = argv[5];
	int key_column = stoi(argv[2]);
	int value_column = stoi(argv[3]);
	int num_modes = stoi(argv[4]);
	double alpha = 0.8;
	if (argc >= 6)
		alpha = stod(argv[6]);
		


//	int key_column = 3, value_column = 4, num_modes = 4;
	auto data = readData(input_file, key_column, value_column);
	auto res = propagate(data, num_modes, alpha);
	writeData(prefix, res);
}

