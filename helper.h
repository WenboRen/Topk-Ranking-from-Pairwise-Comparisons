#pragma once

#include<random>
#include<unordered_map>
#include<fstream>
#include<sstream>
#include<string>
#include<vector>
using namespace std;


bool compare(const double p) {
	return rand() < p * RAND_MAX;
}

void randomReorder(vector<int> &S) {
	int n = S.size();
	for (int i = 0; i < n; ++i) {
		swap(S[i], S[i + rand() % (n - i)]);
	}
}

vector<int> randomSelect(const vector<int> &S, int m) {
	int n = S.size();
	vector<int> ans;
	for (int i = 0; i < S.size() && m > 0; ++i, --n) {
		if (rand() % n < m) {
			ans.push_back(S[i]);
			--m;
		}
	}
	return ans;
}

bool isPACBestK(const vector<int> &R, const vector<int> &S, const vector<vector<double>> &P, const double epsilon) {
	int k = R.size(), n = S.size();
	unordered_map<int, bool> is_in_R;
	for (int i : R) {
		is_in_R[i] = true;
	}

	for (int i = 0; i < k; ++i) {
		for (int j = 0; j < n; ++j) {
			if (is_in_R[S[j]]) {
				continue;
			}
			if (P[R[i]][S[j]] < 0.5 - epsilon) {
				return false;
			}
		}
	}
	return true;
}

double doubleRand(const double lo, const double hi) {
	return lo + (double)rand() * (hi - lo) / (double)RAND_MAX;
}

vector<vector<double>> PrefLibPWGLoader(const string &path) {
	ifstream ifile(path);
	string line;
	vector<vector<int>> wins;
	int buff[3];
	int num_items;

	if (ifile.is_open()) {
		if (!getline(ifile, line)) {
			ifile.close();
			cout << "Error: Data format error";
		}
		num_items = stoi(line);

		for (int i = 0; i < num_items + 1; ++i) {
			if (!getline(ifile, line)) {
				ifile.close();
				cout << "Error: Data format error";
			}
		}

		wins = vector<vector<int>>(num_items, vector<int>(num_items));

		while (getline(ifile, line)) {
			int start = 0, end;
			for (int i = 0; i < 3; ++i) {
				end = line.find(',', start);
				if (end == string::npos) {
					end = line.size();
				}
				buff[i] = stoi(line.substr(start, end));
				start = end + 1;
			}
			wins[buff[1] - 1][buff[2] - 1] = buff[0];
		}

		ifile.close();
	}
	else {
		cout << "Error: Unable to open the file " + path + "!" << endl;
		return{};
	}

	int num = 0;
	for (int i = 0; i < num_items; ++i) {
		for (int j = 0; j < num_items; ++j) {
			num += wins[i][j];
		}
	}
	cout << "num = " << num << endl;

	vector<vector<double>> P(num_items, vector<double>(num_items));
	for (int i = 0; i < num_items; ++i) {
		for (int j = i + 1; j < num_items; ++j) {
			P[i][j] = (double)(wins[i][j]) / (double)(wins[i][j] + wins[j][i]);
			P[j][i] = 1 - P[i][j];
		}
		P[i][i] = 0.5;
	}

	return P;
}