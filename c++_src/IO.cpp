/*
 * IO.cpp
 *
 *  Created on: Apr 25, 2016
 *      Author: wdk
 */

#include "IO.hpp"
#include <dirent.h>
#include <errno.h>
//using namespace std;

string get_output_filename(string prefix, int &N, int &K, double &EPSILON,
		double &EPSILON2, string suffix) {
	return prefix + "_num_" + to_string(N) + "_out_" + to_string(K) + "_eps_"
			+ to_string(EPSILON) + "_eps2_" + to_string(EPSILON2) + suffix;
}

void get_files_in_dir(string dir, vector<string> &files) {
	DIR *dp;
	struct dirent *dirp;
	if ((dp = opendir(dir.c_str())) == NULL) {
		cerr << "Error(" << errno << ") opening " << dir << endl;
	}

	while ((dirp = readdir(dp)) != NULL) {
		char *filename = dirp->d_name;
		int len = strlen(filename);
		if (len > 4 && strcmp(filename + len - 4, ".txt") == 0
				&& strncmp(filename, GOAL_FILE_PREFIX, strlen(GOAL_FILE_PREFIX))
						== 0)
			files.push_back(string(filename));
	}
	closedir(dp);
}

void read_outliners(set<int> &outliners, char const *filename) {
	int t;
	ifstream fin;
	fin.open(filename);
	while (fin >> t) {
		outliners.insert(t);
	}
}

void read_file(int &N, int &D, int &K, double &EPSILON, mat &DISTANCES,
		char const *filename) {
	FILE *fp;
	fp = fopen(filename, "r");
	fscanf(fp, "%d %d %d %lf ", &N, &D, &K, &EPSILON);

	DISTANCES = zeros<mat>(N, N);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			fscanf(fp, "%lf", &DISTANCES(i, j));
		}
	}
	fclose(fp);
}

void read_file(int &N, int &D, int &K, double &EPSILON, vector<int> &labels,
		mat &DISTANCES, char const *filename) {
	FILE *fp;
	fp = fopen(filename, "r");
	fscanf(fp, "%d %d %d %lf ", &N, &D, &K, &EPSILON);

	labels.resize(N, 0);
	for (int i = 0; i < N; i++) {
		fscanf(fp, "%d", &labels[i]);
	}

	DISTANCES = zeros<mat>(N, N);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			fscanf(fp, "%lf", &DISTANCES(i, j));
		}
	}
	fclose(fp);
}

void write_miscellaneous(int N, vector<double> norms, Embedding &ret,
		const string &filename) {
	ofstream fout;
	bool good[N];
	memset(good, false, sizeof(good));

	fout.open(filename.c_str());
//	fout << "norm2: " << norm2 << endl;
//	fout << "normalized norm2: " << norm2 / sqrt(size * (size - 1) / 2) << endl;
	fout << "norms(norm_max, norm_1,avg_norm1,norm2,avg_norm2):";
	for(auto &norm : norms) {
		fout << " " << norm;
	}
	for (auto &k : ret.embedded_point_indexes) {
		good[k] = true;
	}
	cerr << N << " " << ret.size();
	fout << endl << N - ret.size() << " outliners found: " << endl;
	for (int i = 0; i < N; i++) {
		if (!good[i]) {
			fout << i << endl;
		}
	}
	fout.close();
}

void write_miscellaneous(int N, double norm2, set<int> &real_outliners,
		Embedding &ret, const string &filename, double &a, double &b) {
	ofstream fout;
	int size = ret.size();
	bool good[N];
	memset(good, false, sizeof(good));

	fout.open(filename.c_str());
	fout << "norm2: " << norm2 << endl;
	fout << "normalized norm2: " << norm2 / sqrt(size * (size - 1) / 2) << endl;
	for (auto &k : ret.embedded_point_indexes) {
		good[k] = true;
	}
	double cnt_rout = real_outliners.size();
	double cnt_fout = N - ret.size();
	double cnt_frout = 0.0;
	fout << endl << cnt_fout << " outliners found: " << endl;
	for (int i = 0; i < N; i++) {
		if (!good[i]) {
			if (real_outliners.find(i) != real_outliners.end()) {
				cnt_frout++;
			}
			fout << i << endl;
		}
	}
	a = cnt_frout / cnt_fout;
	b = cnt_frout / cnt_rout;
	fout << "real outliner in found rate : " << a << endl;
	fout << "real_outliner not found rate: " << b << endl;
	fout.close();
}

void write_positions(Embedding &ret, const string &filename) {
	ofstream fout;
	fout.open(filename.c_str());
	// fout << ret.size() << endl;
	cout << "Finished, " << ret.size() << " points embeddable." << endl;
	int row = ret.embedded_positions.n_rows;
	int col = ret.embedded_positions.n_cols;
	for(int i = 0; i<ret.size(); i++){
		fout << ret.embedded_point_indexes[i] + 1; // matlab start from 1
		if(i != ret.size() - 1)
			fout << ",";
	}
	fout << endl;
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			fout << ret.embedded_positions(i, j);
			if (j != col - 1)
				fout << ",";
		}
		fout << endl;
	}
	fout.close();
}

void write_positions(Embedding &ret, const string &filename,
		vector<int> labels) {
	ofstream fout;
	fout.open(filename.c_str());
	// fout << ret.size() << endl;
	cout << "Finished, " << ret.size() << " points embeddable." << endl;
	int row = ret.embedded_positions.n_rows;
	int col = ret.embedded_positions.n_cols;
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			fout << ret.embedded_positions(i, j);
			if (j != col - 1)
				fout << ",";
		}
		fout << endl;
	}
	for (int i = 0; i < col; i++) {
		fout << labels[ret.embedded_point_indexes[i]];
		if (i != col - 1)
			fout << ",";
	}
	fout << endl;
	fout.close();
}

void write_rates(vector<double> &rate_as, vector<double> &rate_bs,
		const char *filename) {
	ofstream fout;
	fout.open(filename);
	int size = rate_as.size();
	for (int i = 0; i < size; i++) {
		fout << rate_as[i] << "," << rate_bs[i] << endl;
	}
}

