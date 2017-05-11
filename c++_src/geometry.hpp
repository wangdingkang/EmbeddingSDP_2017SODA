/*
 * geometry.hpp
 *
 *  Created on: Apr 25, 2016
 *      Author: wdk
 */

#ifndef GEOMETRY_HPP_
#define GEOMETRY_HPP_

#include <bits/stdc++.h>
#include <armadillo>
using namespace std;
using namespace arma;

inline double pairwise_distance(const vec &v1, const vec &v2) {
	return norm(v1 - v2, 2);
}

void solve_quadratic(double a, double b, double c, double &r1, double &r2);
double get_distance(mat &position_matrix, vec &q);
void get_embedding(vec &v1, vec &v2, const vec &vb, vec &v_ret, double distance,
		double &d_ret);
void embed_fixed_points(mat &squared_distances, mat &ret_positions);

vec embed_new_point(mat &fixed_positions,
		mat &fixed_square_distances, vec &q_squared_distances);

pair<vec, vec> embed_new_point(int p_kick, mat &fixed_positions,
		mat &fixed_square_distances, vec &q_squared_distances);

#endif /* GEOMETRY_HPP_ */
