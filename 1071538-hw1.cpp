#include<iostream>
#include<string>
#include<fstream>
#include<vector>
#include<cmath>
#include<iomanip>
#include<cstdlib>
using namespace std;

int main(int argc, char*argv[])
{
	ifstream file2("1.DNA");
	ifstream file1("1.PROTEIN");
	ifstream file3("1.MATRIX");
	/*if (!file) {
		cout << "無法開啟檔案1!!\n";
		system("pause");
		exit(0);
	}*/
	/*ifstream file1;
	ifstream file2;
	ifstream file3;
	file1.open(argv[1]);
	file2.open(argv[2]);
	file3.open(argv[3]);*/
	string world;
	char w[20];

	vector<double>buf;
	vector<vector<double> > dnaPosition;
	vector<vector<double> > proPosition;
	vector<vector<double> > matrix;
	vector<double>proPosition_x;
	vector<double>proPosition_y;
	vector<double>proPosition_z;
	vector<int>proPosition_aa;
	vector<double> Protein_DNA;
	vector<int>  binding_sites;
	while (getline(file1, world)) {
		for (int i = 22; i < 54; i++) {
			if (world[i] != '\0') {
				int m = 0;
				while (world[i] != '\0') {
					w[m] = world[i];
					m++;
					if (i == 26 || i == 38 || i == 46 || i == 54)
						break;
					i++;
				}
				buf.push_back(atof(w));
				for (int i = 0; i < 20; i++)
					w[i] = ' ';
			}
		}
		proPosition.push_back(buf);
		buf.clear();
	}
	while (getline(file2, world)) {
		for (int i = 22; i < 54; i++) {
			if (world[i] != '\0') {			
				int m = 0;
				while (world[i] != '\0' ) {
					w[m] = world[i];
					m++;
					if (i == 26 || i == 38 || i == 46 || i == 54 )
						break;
					i++;
				}
				buf.push_back(atof(w));
				for (int i = 0; i < 20; i++)
					w[i] = ' ';
			}					
		}
		dnaPosition.push_back(buf);
		buf.clear();
	}



	int e = 0;
	while (getline(file3, world)) {	
			if(e<6&&2<e){
				for (int i = 2; i < world.length(); i++) {
					if (world[i] != '\0'&& world[i]!='\t' && world[i] != ' ') {
						int m = 0;
						while (world[i] != '\0' && world[i] != '\t' && world[i] != '\n' && world[i] != ' ') {
							w[m] = world[i];
							m++;
							i++;
						}
						buf.push_back(atof(w));
						for (int i = 0; i < 20; i++)
							w[i] = ' ';
					}
				}
				matrix.push_back(buf);
				buf.clear();
		}
			e++;			
	}
	
	for (int i = 0; i < proPosition.size(); i++) {
		proPosition_aa.push_back(proPosition[i][0]);
		proPosition_x.push_back(matrix[0][0] + (matrix[0][1] * proPosition[i][1]) + (matrix[0][2] * proPosition[i][2]) + (matrix[0][3] * proPosition[i][3]));
		proPosition_y.push_back( matrix[1][0] + (matrix[1][1] * proPosition[i][1]) + (matrix[1][2] * proPosition[i][2]) + (matrix[1][3] * proPosition[i][3]));
		proPosition_z.push_back( matrix[2][0] + (matrix[2][1] * proPosition[i][1]) + (matrix[2][2] * proPosition[i][2]) + (matrix[2][3] * proPosition[i][3]));
		for (int j = 0; j < dnaPosition.size(); j++) {
					if ((sqrt(pow((proPosition_x[i] - dnaPosition[j][1]), 2) + pow((proPosition_y[i] - dnaPosition[j][2]), 2) + pow((proPosition_z[i] - dnaPosition[j][3]), 2))) < 4.5) {
						if (binding_sites.size() == 0) {
							binding_sites.push_back(proPosition_aa[i]);
						}
						else if (binding_sites[binding_sites.size() - 1] != proPosition_aa[i]) {
							binding_sites.push_back(proPosition_aa[i]);
						}
						break;
					}
				}
	}

	cout << binding_sites.size() << endl;
	for (int i = 0; i < binding_sites.size(); i++) {
		if (i != binding_sites.size() - 1)
			cout << binding_sites[i] << " ";
		else
			cout << binding_sites[i];
	}
	cout << endl;
	file1.close();
	file2.close();
	file3.close();
	//system("pause");
}