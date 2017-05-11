/*
 * sdp.hpp
 *
 *  Created on: Apr 25, 2016
 *      Author: wdk
 */

#ifndef SDP_HPP_
#define SDP_HPP_

#include "./structs.hpp"
#include "scs.h"
#include "linsys/amatrix.h"

#define SQRT2 1.41421356237

void parse_res(int D, Sol* &sol, mat &save_to);
void parse_q_distances(int D, Sol* &sol, vec &save_to);
void standard_malloc(Data* &data, Cone* &cone, Sol* &sol, Info* &info,
		AMatrix* &A, int m, int n, int nnz);
void init_xij_problem(Data* &data, Cone* &cone, Sol* &sol, Info* &info, int D);
void init_xij_a_problem(Data* &data, Cone* &cone, Sol* &sol, Info* &info,
		int D);
void init_xij_barray(int D, vector<int> &B, double EPSILON, mat &DISTANCES,
		Data* &data);
void init_xij_barray(int D, vector<int> &B, double EPSILON, mat &DISTANCES,
		Data* &data, vector<int> &index_to_change);
void init_xij_a_barray(int D, vector<int> &B, double EPSILON,
		mat &fixed_squared_distances, Data* &data,
		vector<int> &index_to_change);
void init_xij_a_barray(int D, int p_kick, int q, vector<int> &B, double EPSILON,
		mat &fixed_squared_distances, mat &DISTANCES, Data* &data);

inline void update_barray(Data* &data, vector<int> &index_to_change,
		vector<scs_float> &constraint_values) {
	for (int i = 0; i < (int) index_to_change.size(); i++) {
		data->b[index_to_change[i]] = constraint_values[i];
	}
}

void standard_release(Data* &data, Cone* &cone, Sol* &sol, Work* &work);
void standard_release(Data* &data, Cone* &cone, Sol* &sol);
void sdp_squared_distances_xi_xj(int N, double EPSILON, vector<int> &B,
		mat &DISTANCES, mat &fixed_square_distances);

#endif /* SDP_HPP_ */
