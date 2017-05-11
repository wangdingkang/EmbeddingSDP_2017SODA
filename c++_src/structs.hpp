/*
 * structs.hpp
 *
 *  Created on: Apr 25, 2016
 *      Author: wdk
 */

#ifndef STRUCTS_HPP_
#define STRUCTS_HPP_

#include <bits/stdc++.h>
#include <armadillo>
using namespace std;
using namespace arma;

struct Embedding {
	vector<int> embedded_point_indexes;
	vector<int> embedded_sizes;
	mat embedded_positions;

	Embedding() {
	}

	Embedding(vector<int> sizes, vector<int> &indexes, mat positions) {
		this->embedded_sizes = sizes;
		this->embedded_point_indexes = indexes;
		this->embedded_positions = positions;
	}

	int size() {
		return (int) embedded_point_indexes.size();
	}

};

struct CCF_Element {
	int col, row;
	double v;

	CCF_Element(int _col, int _row, double _v) {
		col = _col;
		row = _row;
		v = _v;
	}

	bool operator <(const CCF_Element &another) const {
		if (col == another.col) {
			return row < another.row;
		}
		return col < another.col;
	}

	void show_element() {
		cout << "Element: row: " << row << " col: " << col << " value: " << v
				<< endl;
	}

};

#endif /* STRUCTS_HPP_ */
