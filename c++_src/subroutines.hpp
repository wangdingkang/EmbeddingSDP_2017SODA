/*
 * subroutines.hpp
 *
 *  Created on: Apr 25, 2016
 *      Author: wdk
 */

#ifndef SUBROUTINES_HPP_
#define SUBROUTINES_HPP_

#include "structs.hpp"
#include "geometry.hpp"
#include "sdp.hpp"

void gen_random_perm(int N, vector<int> &B, vector<int> &random_permutation);
vector<pair<int, int> > get_conflict_edges(mat &positions, vector<int> &indexes,
		mat &distances, double epsilon);
vector<int> two_approximate_vertex_cover(vector<pair<int, int> > &edges);
void first_step(vector<int> &perm, vector<int> &B, vector<vector<int> > &U, double EPSILON, int K, int D,
		mat &DISTANCES);
Embedding second_step(vector<int> &B, vector<vector<int> > &U, int N, double EPSILON, double EPSILON2,
		mat &DISTANCES);
double parse_distance(mat &positions, vec &embedded_point, Sol *sol);
void norm_distance(Embedding &ret, vector<double> &norms, mat &DISTANCES);

#endif /* SUBROUTINES_HPP_ */
