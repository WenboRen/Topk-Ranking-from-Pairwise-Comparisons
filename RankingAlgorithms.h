#pragma once

#include<vector>
#include<algorithm>
#include<functional>
using namespace std;

#include"helper.h"

#define PI 3.1415926

//mode = false, best selection
//mode = true, worst selection

void distributeItem(const int i, const double piv, const double epsilon, const double su, const double sd, const double delta, vector<int> &Sup, vector<int> &Smiddle, vector<int> &Sdown, double &complexity, const bool mode) {
	double tmax = ceil((2.0 / (epsilon * epsilon)) * log(4.0 / delta));
	double w = 0.0, t = 0.0;

	while (t < tmax) {
		t += 1.0;
		if (mode ^ compare(piv)) {
			++w;
		}
		complexity += 1.0;
		double bt = sqrt((0.5 / t) * log((PI * PI * t * t) / (3.0 * delta)));
		if (w / t - bt > 0.5 + su) {
			Sup.push_back(i);
			return;
		}
		else if (w / t + bt < 0.5 - sd) {
			Sdown.push_back(i);
			return;
		}
	}

	if (w / t > 0.5 + 0.5 * epsilon + su) {
		Sup.push_back(i);
	}
	else if (w / t < 0.5 - 0.5 * epsilon - sd) {
		Sdown.push_back(i);
	}
	else {
		Smiddle.push_back(i);
	}
}

vector<int> epsilonQuickSelect(const vector<int> &S, const vector<vector<double>> &P, const int k, const double epsilon, const double delta, double &complexity, const bool mode) {
	int n = S.size(), v = S[rand() % n];
	vector<int> Sup, Smiddle, Sdown;
	Smiddle.push_back(v);
	double delta_1 = delta / (n * (n - 1));

	for (int i = 0; i < n; ++i) {
		if (S[i] == v) {
			continue;
		}
		distributeItem(S[i], P[S[i]][v], epsilon / 2.0, 0.0, 0.0, delta_1, Sup, Smiddle, Sdown, complexity, mode);
	}

	if (Sup.size() > k) {
		return epsilonQuickSelect(Sup, P, k, epsilon, (n - 1) * delta / n, complexity, mode);
	}
	if (Sup.size() + Smiddle.size() >= k) {
		int kp = k - Sup.size();
		for (int i = 0; i < kp; ++i) {
			Sup.push_back(Smiddle[i]);
		}
		return Sup;
	}
	int kp = k - Sup.size() - Smiddle.size();
	vector<int> Sp = epsilonQuickSelect(Sdown, P, kp, epsilon, (n - 1) * delta / n, complexity, mode);
	Sup.insert(Sup.end(), Smiddle.begin(), Smiddle.end());
	Sup.insert(Sup.end(), Sp.begin(), Sp.end());
	return Sup;
}

vector<int> tournamentKSelection(const vector<int> &S, const vector<vector<double>> &P, const int k, const double epsilon, const double delta, double &complexity, bool mode) {
	int t = 0, n = S.size(), nRemain = n, m = 2 * k;
	double epsilon_t = 0.25 * epsilon, delta_t;
	vector<int> R = S;
	vector<vector<int>> splited;
	while (R.size() > k) {
		++t;
		epsilon_t *= 0.8;
		delta_t = 6.0 * delta / (PI * PI * t * t);
		splited.clear();
		for (int i : R) {
			if (splited.empty() || splited.back().size() == m) {
				splited.push_back({});
			}
			splited.back().push_back(i);
		}
		R.clear();
		for (vector<int> &part : splited) {
			if (part.size() <= k) {
				R.insert(R.end(), part.begin(), part.end());
			}
			else {
				vector<int> A = epsilonQuickSelect(part, P, k, epsilon_t, delta_t / k, complexity, mode);
				R.insert(R.end(), A.begin(), A.end());
			}
		}
	}

	return R;
}

int sequentialEliminationExactBestSelection(const vector<int> &S, const vector<vector<double>> &P, const double delta, double &complexity, const bool mode) {
	double epsilon_t = 1.0, delta_t;
	vector<int> R = S, Sup, Smiddle, Sdown;
	int t = 1;
	while (R.size() > 1) {
		epsilon_t *= 0.5;
		delta_t = 6.0 * delta / (PI * PI * t * t);
		int v = tournamentKSelection(R, P, 1, epsilon_t / 3.0, 2.0 * delta_t / 3.0, complexity, mode).back();
		Sup.clear(); Smiddle.clear(); Sdown.clear();
		Smiddle.push_back(v);
		for (int i : R) {
			if (i == v) {
				continue;
			}
			distributeItem(i, P[i][v], epsilon_t / 3.0, 0.0, epsilon_t / 3.0, delta_t / 3.0, Sup, Smiddle, Sdown, complexity, mode);
		}
		R = Sup;
		R.insert(R.end(), Smiddle.begin(), Smiddle.end());
		++t;
	}
	return R.back();
}

vector<int> sequentialEliminationExactKSelection(const vector<int> &S, const vector<vector<double>> &P, const int k, const double delta, double &complexity, const bool mode) {
	double epsilon_t = 1.0, delta_t;
	vector<int> R_t = S, S_t, Sup, Smiddle, Sdown;
	int t = 1, k_t = k;

	while (S_t.size() < k_t && S_t.size() + R_t.size() > k_t) {
		epsilon_t *= 0.5;
		delta_t = 6.0 * delta / (PI * PI * t * t);
		vector<int> A_t = tournamentKSelection(R_t, P, k_t, epsilon_t / 3.0, delta_t / 3.0, complexity, mode);
		int v_t = tournamentKSelection(A_t, P, 1, epsilon_t / 3.0, delta_t / 3.0, complexity, !mode).back();
		Sup.clear(); Smiddle.clear(); Sdown.clear();
		Smiddle.push_back(v_t);

		for (int i : R_t) {
			if (i == v_t) {
				continue;
			}
			distributeItem(i, P[i][v_t], epsilon_t / 3.0, epsilon_t / 3.0, epsilon_t / 3.0, delta_t / (3.0 * (R_t.size() - 1)), Sup, Smiddle, Sdown, complexity, mode);
		}
		S_t.insert(S_t.end(), Sup.begin(), Sup.end());
		R_t = Smiddle;
		k_t -= Sup.size();
		++t;
	}

	S_t.insert(S_t.end(), R_t.begin(), R_t.begin() + (k - S_t.size()));
	return S_t;
}

vector<int> sequentialEliminationExactKSelectionV2(const vector<int> &S, const vector<vector<double>> &P, const int k, const double delta, double &complexity, const bool mode) {
	double epsilon_t = 1.0, delta_t;
	vector<int> R_t = S, S_t, Sup, Smiddle, Sdown;
	int t = 1, k_t = k;

	while (S_t.size() < k_t && S_t.size() + R_t.size() > k_t) {
		epsilon_t *= 0.5;
		delta_t = 6.0 * delta / (PI * PI * t * t);
		vector<int> A_t = epsilonQuickSelect(R_t, P, k_t, epsilon_t / 3.0, delta_t / 3.0, complexity, mode);
		int v_t = tournamentKSelection(A_t, P, 1, epsilon_t / 3.0, delta_t / 3.0, complexity, !mode).back();
		Sup.clear(); Smiddle.clear(); Sdown.clear();
		Smiddle.push_back(v_t);

		for (int i : R_t) {
			if (i == v_t) {
				continue;
			}
			distributeItem(i, P[i][v_t], epsilon_t / 3.0, epsilon_t / 3.0, epsilon_t / 3.0, delta_t / (3.0 * (R_t.size() - 1)), Sup, Smiddle, Sdown, complexity, mode);
		}
		S_t.insert(S_t.end(), Sup.begin(), Sup.end());
		R_t = Smiddle;
		k_t -= Sup.size();
		++t;
	}

	S_t.insert(S_t.end(), R_t.begin(), R_t.begin() + (k - S_t.size()));
	return S_t;
}

bool compareSub1(const double pij, const double epsilon, const double delta, double &complexity, const bool mode) {
	double p_i = 0.5, c = 0.5, m = (0.5 / (epsilon * epsilon)) * log(2.0 / delta), r = 0.0, w_i = 0.0;
	while (fabs(p_i - 0.5) < c - epsilon && r <= m) {
		if (compare(pij) ^ mode) {
			w_i += 1.0;
		}
		complexity += 1.0;
		r += 1.0;
		p_i = w_i / r;
		c = sqrt((0.5 / r) * log(4.0 * r * r / delta));
	}
	if (p_i <= 0.5) {
		return false;
	}
	return true;
}

vector<int> knockoutRound(const vector<int> &S, const vector<vector<double>> &P, const double epsilon, const double delta, double &complexity, const bool mode) {
	int pos = 0, n = S.size();
	vector<int> top;
	if (n % 2) {
		top.push_back(S[0]);
		pos = 1;
	}
	while (pos < n) {
		if (compareSub1(P[S[pos]][S[pos + 1]], epsilon, delta, complexity, mode)) {
			top.push_back(S[pos]);
		}
		else {
			top.push_back(S[pos + 1]);
		}
		pos += 2;
	}
	return top;
}

int knockout(const vector<int> &S, const vector<vector<double>> &P, const double epsilon, const double delta, double &complexity, const bool mode) {
	double gamma = 1.0, t = 1.0;
	vector<int> R = S;
	double c = pow(2.0, 1.0 / 3.0) - 1;
	while (R.size() > 1) {
		R = knockoutRound(R, P, c * epsilon / (gamma * pow(2.0, t / 3.0)), delta / pow(2.0, t), complexity, mode);
		t += 1.0;
	}
	return R[0];
}

int compareSub2(const double pij, const double epsilon_l, const double epsilon_u, const double delta, double &complexity, const bool mode) {
	double epsilon_m = (epsilon_l + epsilon_u) / 2.0, p = 0.0, c = 0.5, t = 0.0, w = 0.0;
	double tmax = (2.0 / (epsilon_u - epsilon_l) / (epsilon_u - epsilon_l)) * log(2.0 / delta);
	while (fabs(p - epsilon_m) <= c && t <= tmax) {
		if (compare(pij) ^ mode) {
			w += 1.0;
		}
		complexity += 1.0;
		t += 1.0;
		p = w / t - 0.5;
		c = sqrt((0.5 / t) * log(4.0 * t * t / delta));
	}
	if (p <= epsilon_m) {
		return 1;
	}
	return 2;
}

int seqEliminate(const vector<int> &S, const vector<vector<double>> &P, const double epsilon, const double delta, double &complexity, const bool mode) {
	int n = S.size(), r = 0, c = 1;
	while (c < n) {
		if (compareSub2(P[S[c]][S[r]], 0.0, epsilon, delta / n, complexity, mode) == 2) {
			r = c;
		}
		++c;
	}
	return S[r];
}

int pickAnchor(const vector<int> &S, const vector<vector<double>> &P, const double np, const double epsilon, const double delta, double &complexity, const bool mode) {
	int m = min<int>(ceil((S.size() / np) * log(2.0 / delta)), S.size());
	vector<int> Q = randomSelect(S, m);
	return seqEliminate(Q, P, epsilon, delta / 2.0, complexity, mode);
}

vector<int> prune(const vector<int> &S, const vector<vector<double>> &P, const int a, const double np, const double epsilon_l, const double epsilon_u, const double delta, double &complexity, const bool mode) {
	double t = 1.0, threshold = pow(log(S.size()), 2);
	vector<int> R = S, nR;
	while (R.size() > 2 * np && t < threshold) {
		nR.clear();
		for (int e : R) {
			if (compareSub2(P[e][a], epsilon_l, epsilon_u, delta * pow(0.5, 1 + t), complexity, mode) != 1) {
				nR.push_back(e);
			}
		}
		t += 1.0;
		R.swap(nR);
	}
	return R;
}

int optMaximize(const vector<int> &S, const vector<vector<double>> &P, const double epsilon, const double delta, double &complexity, const bool mode) {
	int n = S.size();
	if (delta < 1.0 / n) {
		return seqEliminate(S, P, epsilon, delta, complexity, mode);
	}

	int a = pickAnchor(S, P, sqrt(6.0 * n * log((double)n)), epsilon / 3.0, delta / 4.0, complexity, mode);
	vector<int> Sp = prune(S, P, a, sqrt(6.0 * n * log((double)n)), epsilon / 3.0, 2.0 * epsilon / 3.0, delta / 4.0, complexity, mode);

	for (int e : Sp) {
		if (compareSub2(P[e][a], 2.0 * epsilon / 3.0, epsilon, 0.25 * delta / n, complexity, mode) == 2) {
			return seqEliminate(Sp, P, epsilon, delta / 4.0, complexity, mode);
		}
	}
	return a;
}


vector<int> activeRankingForTopK(const vector<int> &S, vector<vector<double>> &P, const int k, const double delta, double &complexity, const bool mode) {
	int n = S.size();
	vector<int> R, T, remain(n, 0);
	vector<double> estimated(n, 0.0);
	double t = 0.0;
	vector<bool> assured(n, false);
	int numAssured = 0;
	vector<pair<double, int>> fs;
	while (R.size() < k) {
		t += 1.0;
		for (int i = 0; i < n; i++) {
			if (assured[i]) {
				continue;
			}
			if (compare(P[S[i]][S[rand() % n]])) {
				estimated[i] = (1.0 + estimated[i] * (t - 1.0)) / t;
			}
			else {
				estimated[i] = estimated[i] * (t - 1.0) / t;
			}
			complexity += 1.0;
		}
		fs.clear();
		for (int i = 0; i < n; i++) {
			if (!assured[i]) {
				fs.push_back(make_pair(estimated[i], i));
			}
		}
		sort(fs.begin(), fs.end());
		double alpha_t = sqrt(log(125.0 * n * log(1.12 * t) / delta) / t);
		double anchorLow = fs[n - 1 - numAssured - (k - R.size())].first;
		double anchorHigh = fs[n - numAssured - (k - R.size())].first;
		for (int i = 0; i < n; i++) {
			if (assured[i]) {
				continue;
			}
			if (estimated[i] > anchorLow + 4.0 * alpha_t) {
				R.push_back(S[i]);
				assured[i] = true;
				numAssured++;
				if (R.size() >= k) {
					return R;
				}
			}
			else if (estimated[i] < anchorHigh - 4.0 * alpha_t) {
				T.push_back(S[i]);
				assured[i] = true;
				numAssured++;
				if (T.size() >= n - k) {
					for (int i = 0; i < n; i++) {
						if (!assured[i]) {
							R.push_back(S[i]);
						}
					}
					return R;
				}
			}
		}
	}
	return R;
}


//Cite from "Beat the mean bandit (2011). Yisong Yue et al."
int beatTheMeanBandit(const vector<int> &S, const vector<vector<double>> &P, const double delta, const double epsilon, double &complexity, const bool mode) {
	int n = S.size();

	double N_0 = 100;
	const double N = ceil((36.0 / (epsilon * epsilon)) * log(pow(n, 3) * N_0 / delta));
	const double log_delta = log(pow(n, 3) * N / delta);

	double gamma2 = 1.0;

	vector<vector<int>> compares(n, vector<int>(n, 0));
	vector<int> nums(n, 0);
	vector<vector<int>> wins(n, vector<int>(n, 0));
	vector<double> scores(n, 0.0);
	vector<double> c(n);

	vector<int> R;
	for (int i = 0; i < S.size(); ++i) {
		R.push_back(i);
	}

	while (R.size() > 1) {
		double min_ucb = 10000, max_lcb = -10000;
		for (int b : R) {
			int bp = b;
			while (bp == b) {
				bp = R[rand() % R.size()];
			}
			++nums[b];
			++compares[b][bp];
			if (compare(P[S[b]][S[bp]]) ^ mode) {
				++wins[b][bp];
				scores[b] += (1.0 - scores[b]) / (double)nums[b];
			}
			else {
				scores[b] += (0.0 - scores[b]) / (double)nums[b];
			}
			c[b] = 3.0 * gamma2 * sqrt(log_delta / (double)nums[b]);
			min_ucb = min(min_ucb, scores[b] + c[b]);
			max_lcb = max(max_lcb, scores[b] - c[b]);
		}
		complexity += R.size();

		if (min_ucb <= max_lcb) {
			int min_arm = min_element(scores.begin(), scores.end()) - scores.begin();
			for (int b : R) {
				if (b == min_arm) {
					continue;
				}
				scores[b] = (scores[b] * nums[b] - wins[b][min_arm]) / (double)(nums[b] - compares[b][min_arm]);
				nums[b] -= compares[b][min_arm];
			}
			int index = 0;
			for (index = 0; R[index] != min_arm; ++index);
			R.erase(R.begin() + index);
			scores[min_arm] = 10000;
		}
	}

	return S[R.back()];
}

// From Preference-based rank elicitation using statistical models: The case of mallows
int MallowsMPI(const vector<int> &S, const vector<vector<double>> &P, const double delta, double &complexity, const bool mode) {
	int a = S[0];
	for (int i = 1; i < S.size(); ++i) {
		double bound, estimate = 0.0;
		int num = 0, b = S[i];
		do {
			++num;
			if (compare(P[a][b]) ^ mode) {
				estimate += (1.0 - estimate) / num;
			}
			else {
				estimate += (0.0 - estimate) / num;
			}
			bound = sqrt(log(4 * S.size() * num * num / delta) / (2.0 * num));
		} while (fabs(estimate - 0.5) <= bound);

		complexity += num;

		if (estimate < 0.5) {
			a = b;
		}
	}

	return a;
}
