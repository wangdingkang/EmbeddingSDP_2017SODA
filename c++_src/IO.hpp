/*
 * IO.hpp
 *
 *  Created on: Apr 25, 2016
 *      Author: wdk
 */

#ifndef IO_HPP_
#define IO_HPP_

#include "structs.hpp"

#define GOAL_FILE_PREFIX "OD_"

string get_output_filename(string prefix, int &N, int &K, double &EPSILON,
		double &EPSILON2, string suffix);
void get_files_in_dir(string dir, vector<string> &files);
void read_outliners(set<int> &outliners, char const *filename);
void read_file(int &N, int &D, int &K, double &EPSILON, mat &DISTANCES,
		char const *filename);
void read_file(int &N, int &D, int &K, double &EPSILON, vector<int> &labels,
		mat &DISTANCES, char const *filename);
void write_miscellaneous(int N, vector<double> norms, Embedding &ret,
		const string &filename);
void write_miscellaneous(int N, double norm2, set<int> &outlines,
		Embedding &ret, const string &filename, double &a, double &b);
void write_positions(Embedding &ret, const string &filename);
void write_positions(Embedding &ret, const string &filename,
		vector<int> labels);
void write_rates(vector<double> &rate_as, vector<double> &rate_bs,
		const char *filename);

#endif /* IO_HPP_ */
