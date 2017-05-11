/*
 * geometry.cpp
 *
 *  Created on: Apr 25, 2016
 *      Author: wdk
 */

#include "geometry.hpp"

//#define EPS 1e-5

// there must be a solution
void solve_quadratic(double a, double b, double c, double &r1, double &r2) {
	double delta = b * b - 4 * a * c;
	delta = max(0.0, delta);
	delta = sqrt(delta);
	r1 = (-b + delta) / (2 * a);
	r2 = (-b - delta) / (2 * a);
}

double get_distance(mat &position_matrix, vec &q) {
	if (position_matrix.n_rows == 1) {
		return q(0);
	}
	mat difference = diff(position_matrix, 1, 1);
	mat null_space = null(difference.t());
	colvec normal = null_space.col(0);
	return abs(dot(normal, q) / norm(normal, 2));
}

void get_embedding(vec &v1, vec &v2, const vec &vb, vec &v_ret, double distance,
		double &d_ret) {
	double d1 = abs(distance - pairwise_distance(v1, vb));
	double d2 = abs(distance - pairwise_distance(v2, vb));

	if (d1 < d2) {
		v_ret = v1;
		d_ret = d1;
	} else {
		v_ret = v2;
		d_ret = d2;
	}
}

void embed_fixed_points(mat &squared_distances, mat &ret_positions) {
	int n = squared_distances.n_rows;
	ret_positions = zeros<mat>(n - 1, n);
	ret_positions(0, 1) = sqrt(squared_distances(0, 1));
	for (int i = 2; i < n; i++) {
		mat A = zeros<mat>(i - 1, i - 1);
		for (int j = 0; j < i - 1; j++) {
			for (int k = 0; k < i - 1; k++) {
				A(j, k) = -2 * ret_positions(k, j + 1);
			}
		}
		vec b = zeros<vec>(i - 1);
		for (int j = 1; j < i; j++) {
			b(j - 1) = squared_distances(j, i) - squared_distances(0, i)
					- squared_distances(j, 0);
		}
		vec x = solve(A, b);
		double sqx = norm(x, 2);
		sqx = sqx * sqx;
		for (int j = 0; j < i - 1; j++) {
			ret_positions(j, i) = x(j);
		}
		ret_positions(i - 1, i) = sqrt(max(0.0, squared_distances(i, 0) - sqx));
	}
}

vec embed_new_point(mat &fixed_positions, mat &fixed_square_distances,
		vec &q_squared_distances) {
	int n = fixed_square_distances.n_rows;

	vec ret = zeros<vec>(n - 1);
	mat A = zeros<mat>(n - 1, n - 1);
	for (int j = 0; j < n - 1; j++) {
		for (int k = 0; k < n - 1; k++) {
			A(j, k) = -2 * fixed_positions(k, j + 1);
		}
	}
	vec b = zeros<vec>(n - 1);
	for (int j = 1; j < n; j++) {
		b(j - 1) = q_squared_distances(j) - q_squared_distances(0)
				- fixed_square_distances(j, 0);
	}
	vec x = solve(A, b);
	return x;
}

pair<vec, vec> embed_new_point(int p_kick, mat &fixed_positions,
		mat &fixed_square_distances, vec &q_squared_distances) {
	int n = fixed_square_distances.n_rows;
	if (n == 2) {
		vec q1(1), q2(1);
		q1(0) = sqrt(q_squared_distances(0));
		q2(0) = -q1(0);
		if (p_kick == 0) {
			q1(0) += fixed_positions(0, 1);
			q2(0) += fixed_positions(0, 1);
		}
		return make_pair(q1, q2);
	} else {
		mat A = zeros<mat>(n - 2, n - 2);
		vec b = zeros<vec>(n - 2);
		int pivot = (p_kick == 0 ? 1 : 0);
		for (int j = 0, p = pivot + 1; j < n - 2; p++) {
			if (p != p_kick) {
				for (int k = 0; k < n - 2; k++) {
					A(j, k) = -2 * fixed_positions(k, p)
							+ 2 * fixed_positions(k, pivot);
					b(j) = q_squared_distances(j + 1) - q_squared_distances(0)
							- fixed_square_distances(p, 0)
							+ fixed_square_distances(pivot, 0);
				}
				j++;
			}
		}
		vec q1(n - 1), q2(n - 1);
		// solve once
		if (p_kick == n - 1) {
			vec x = solve(A, b);
			double sqx = 0.0;
			for (int i = 0; i < n - 2; i++) {
				sqx += x(i) * x(i);
			}

			double qd = sqrt(max(0.0, q_squared_distances(0) - sqx));

			for (int i = 0; i < n - 2; i++) {
				q1(i) = q2(i) = x(i);
			}
			q1(n - 2) = qd;
			q2(n - 2) = -qd;
			return make_pair(q1, q2);
		}
		// sad story, have to solve qudratic system.
		else {
			vec bs = solve(A, b);
			b(n - 3) += 2 * fixed_positions(n - 2, n - 1);
			vec ks = solve(A, b);
			vec tbs = bs;
			double ta = 1.0, tb = 0.0, tc = 0.0;
			for (int i = 0; i < n - 2; i++) {
				ks(i) -= tbs(i);
				tbs(i) -= fixed_positions(i, pivot);
				ta += (ks(i) * ks(i));
				tb += 2 * ks(i) * tbs(i);
				tc += tbs(i) * tbs(i);
			}
			tc -= q_squared_distances(0);
			double ret1, ret2;
			solve_quadratic(ta, tb, tc, ret1, ret2);
			for (int i = 0; i < n - 2; i++) {
				q1(i) = ks(i) * ret1 + bs(i);
				q2(i) = ks(i) * ret2 + bs(i);
			}
			q1(n - 2) = ret1;
			q2(n - 2) = ret2;
			return make_pair(q1, q2);
		}
	}
}
