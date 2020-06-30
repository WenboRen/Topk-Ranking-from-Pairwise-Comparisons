#include<stdio.h>
#include<tchar.h>
#include<vector>
#include<random>
#include<iostream>
#include<string>
#include<time.h>
using namespace std;

#include"RankingAlgorithms.h"
#include"helper.h"

int main()
{
	double epsilon = 0.001;
	double delta = 0.01;
	int repetition = 100;
	double Delta = 0.1;

	//string Alg = "EQS";
	//string Alg = "TKS";
	//string Alg = "SEEBS";
	//string Alg = "SEEKSv2";
	//string Alg = "SEEKS";
	//string Alg = "knockout";
	//string Alg = "seqEliminate";
	//string Alg = "optMaximize";
	//string Alg = "ActiveRanking";
	//string Alg = "BeatTheMean";
	string Alg = "MallowMPI";

	//string instance = "Homo";
	//string instance = "BTL";
	//string instance = "Random";
	string instance = "PrefLib";
	//string path = "F:/research/ranking/preflib/pwg/ED-00001-00000001.pwg";
	string path = "F:/research/ranking/preflib/pwg/ED-00015-00000047.pwg";

	string check_mode = "PAC";
	//string check_mode = "exact";
	bool worst = false;

	vector<int> Ns, Ks, true_order;
	int m = 0, n_max = 1000;
	vector<vector<double>> P;

	if (instance == "PrefLib") {
		P = PrefLibPWGLoader(path);

		vector<pair<double, int>> borda_scores(P.size());
		for (int i = 0; i < P.size(); ++i) {
			borda_scores[i].second = i;
			borda_scores[i].first = 0.0;
			for (int j = 0; j < P.size(); ++j) {
				borda_scores[i].first += P[i][j] / (double)(P.size());
			}
		}

		Ns = { (int)P.size() };
		Ks = { 1 };
		m = 1;
		n_max = Ns.back();

		true_order = vector<int>(n_max);
		sort(borda_scores.begin(), borda_scores.end(),
			[](const pair<double, int> &p1, const pair<double, int> &p2) {
			return p1.first > p2.first;
		});

		cout << "borda scores:\n";
		for (int i = 0; i < n_max; ++i) {
			true_order[i] = borda_scores[i].second;
			cout << borda_scores[i].second << ": " << borda_scores[i].first << '\n';
		}
		cout << endl;

		double sst = 1.0, sti = 1.0;
		for (int i = 0; i < n_max; ++i) {
			for (int j = i + 1; j < n_max; ++j) {
				for (int k = j + 1; k < n_max; ++k) {
					sst = max(sst, P[true_order[i]][true_order[j]] / P[true_order[i]][true_order[k]]);
					sst = max(sst, P[true_order[j]][true_order[k]] / P[true_order[i]][true_order[k]]);

					sti = max(sti, (P[true_order[i]][true_order[j]] + P[true_order[j]][true_order[k]] - 1.0) / (P[true_order[i]][true_order[k]] - 0.5));
				}
			}
		}
		cout << "k = " << Ks.back() << endl;
		cout << "SST fitting = " << sst << ", STI fitting = " << sti << endl;
	}
	else {
		Ks = { 1, 2,3,4,5, 6,7, 8,10,15,20,25,30,40,50,60,70,80,90,100,120,140,160,180,200,250,300,350, 400,450, 500 };
		Ns = vector<int>(Ks.size(), 1000);
		//Ns = { 2, 5, 10,15,20,25,30,40,50,60,70,80,90,100,120,140,160,180,200,250,300,400,500,600,700,800,900,1000 };
		//Ks = vector<int>(Ns.size(), 1);
		m = Ns.size();

		n_max = 1000;
		P = vector<vector<double>>(n_max, vector<double>(n_max, 0.5));

		true_order = vector<int>(n_max);
		for (int i = 0; i < n_max; ++i) {
			true_order[i] = i;
		}
		randomReorder(true_order);
		if (instance == "Homo") {
			for (int i = 0; i < n_max; ++i) {
				for (int j = i + 1; j < n_max; ++j) {
					P[true_order[i]][true_order[j]] = 0.5 + Delta;
					P[true_order[j]][true_order[i]] = 0.5 - Delta;
				}
			}
		}
		else if (instance == "BTL") {
			vector<double> thetas;
			for (int i = 0; i < n_max; ++i) {
				thetas.push_back(pow(1.1, i));
			}
			reverse(thetas.begin(), thetas.end());
			for (int i = 0; i < n_max; ++i) {
				for (int j = 0; j < n_max; ++j) {
					P[true_order[i]][true_order[j]] = thetas[i] / (thetas[i] + thetas[j]);
				}
			}
		}
		else if (instance == "Random") {
			for (int i = 0; i < n_max; ++i) {
				for (int j = i + 1; j < n_max; ++j) {
					double gap = doubleRand(0.5, 1.5) * Delta;
					P[true_order[i]][true_order[j]] = 0.5 + Delta;
					P[true_order[j]][true_order[i]] = 0.5 - Delta;
				}
			}
		}

		cout << "n:\t";
		for (int n : Ns) {
			cout << n << '\t';
		}
		cout << endl;

		cout << "k:\t";
		for (int k : Ks) {
			cout << k << '\t';
		}
		cout << endl;
	}

	cout << "Alg = " << Alg << " Instance = " << instance << endl;;
	if (instance == "PrefLib") {
		cout << "path = " << path << endl;
	}

	cout << "delta = " << delta << ", epsilon = " << epsilon << endl;


	for (int i = 0; i < m; ++i) {
		int n = Ns[i], k = Ks[i];
		int num_correct = 0;
		double cumulative_complexity = 0.0;

		for (int t = 0; t < repetition; ++t) {
			vector<int> input;
			input.insert(input.end(), true_order.begin(), true_order.begin() + n);
			randomReorder(input);
			vector<int> R;
			if (Alg == "EQS") {
				R = epsilonQuickSelect(input, P, k, epsilon, delta, cumulative_complexity, worst);
			}
			else if (Alg == "TKS") {
				R = tournamentKSelection(input, P, k, epsilon, delta, cumulative_complexity, worst);
			}
			else if (Alg == "SEEBS") {
				R.push_back(sequentialEliminationExactBestSelection(input, P, delta, cumulative_complexity, worst));
			}
			else if (Alg == "SEEKS") {
				R = sequentialEliminationExactKSelection(input, P, k, delta, cumulative_complexity, worst);
			}
			else if (Alg == "SEEKSv2") {
				R = sequentialEliminationExactKSelectionV2(input, P, k, delta, cumulative_complexity, worst);
			}
			else if (Alg == "knockout") {
				R.push_back(knockout(input, P, epsilon, delta, cumulative_complexity, worst));
			}
			else if (Alg == "seqEliminate") {
				R.push_back(seqEliminate(input, P, epsilon, delta, cumulative_complexity, worst));
			}
			else if (Alg == "optMaximize") {
				R.push_back(optMaximize(input, P, epsilon, delta, cumulative_complexity, worst));
			}
			else if (Alg == "ActiveRanking") {
				R = activeRankingForTopK(input, P, k, delta, cumulative_complexity, worst);
			}
			else if (Alg == "BeatTheMean") {
				R = { beatTheMeanBandit(input, P, delta, epsilon, cumulative_complexity, worst) };
			}
			else if (Alg == "MallowMPI") {
				R = { MallowsMPI(input, P, delta, cumulative_complexity, worst) };
			}

			if ((check_mode == "PAC" && isPACBestK(R, input, P, epsilon))
				|| (check_mode == "exact" && isPACBestK(R, input, P, 0.0))) {
				++num_correct;
			}
			//cout << cumulative_complexity / (1 + t) << endl;
		}

		//cout << "Alg = " << Alg << endl;
		//cout << "n = " << n << '\t' << ", k = " << k << endl;
		//cout << "epsilon = " << epsilon << ", delta = " << delta << endl;
		if (instance == "PrefLib")
			cout << "#correct = " << num_correct << endl;
		cout << cumulative_complexity / repetition << '\t';
	}

	while (true) {
		system("pause");
	}
	return 0;
}
