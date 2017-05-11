/*
 * subroutines.cpp
 *
 *  Created on: Apr 25, 2016
 *      Author: wdk
 */

#include "subroutines.hpp"

#define OLD_METHOD

void gen_random_perm(int N, vector<int> &B, vector<int> &random_permutation) {
	random_permutation = vector<int>(N);
	for (int i = 0; i < N; i++) {
		random_permutation[i] = i;
	}
	srand((unsigned) time(NULL));
	random_shuffle(random_permutation.begin(), random_permutation.end());
//	for(auto n : random_permutation) {
//		cerr << n << " ";
//	}
//	cerr << endl;
	B.push_back(random_permutation[N - 1]);
	random_permutation.erase(random_permutation.begin() + N - 1);
}

vector<pair<int, int> > get_conflict_edges(mat &positions, vector<int> &indexes,
		mat &distances, double epsilon) {
	vector<pair<int, int>> ret;
	int size = indexes.size();
	double tdist;
	for (int i = 0; i < size; i++) {
		for (int j = i + 1; j < size; j++) {
			tdist = pairwise_distance(positions.col(i), positions.col(j));
			int ti = indexes[i], tj = indexes[j];
			if (tdist > distances(ti, tj) + epsilon
					|| tdist < distances(ti, tj) - epsilon) {
				ret.push_back(make_pair(ti, tj));
			}
		}
	}
	return ret;
}

vector<int> two_approximate_vertex_cover(vector<pair<int, int> > &edges) {
	set<int> my_set;
	for (auto &p : edges) {
		if (my_set.find(p.first) == my_set.end()
				&& my_set.find(p.second) == my_set.end()) {
			my_set.insert(p.first);
			my_set.insert(p.second);
		}
	}
	vector<int> ret(my_set.begin(), my_set.end());
	return ret;
}

double parse_distance(mat &positions, vec &embedded_point, Sol *sol) {
	int D = positions.n_cols;
	for (int i = 0, k = 0; i < D; i++) {
		for (int j = i; j < D; j++) {
			positions(i, j) = positions(j, i) = sol->x[k++];
		}
	}
	mat U, V;
	vec s;
	svd(U, s, V, positions);
	vec ssqrt = arma::sqrt(s);
	mat S = diagmat(ssqrt);
	positions = S * V;
	embedded_point = positions.col(D - 1);
	positions.col(D - 1) = zeros<vec>(D);
	return get_distance(positions, embedded_point);
}

void first_step(vector<int> &perm, vector<int> &B, vector<vector<int> > &U,
		double EPSILON, int K, int D, mat &DISTANCES) {
	double sqrt_epsilon = sqrt(EPSILON);
	bool outliner[perm.size() + 10];
	memset(outliner, false, sizeof(outliner));
	for (int i = 1; i <= D; i++) {
		vector<pair<double, int> > H;
		vector<int> ui;

		Data *data;
		Sol *sol;
		Cone *cone;
		Info *info;
		double d_low, d_high;

		vector<int> index_to_change;
		init_xij_problem(data, cone, sol, info, i);
		init_xij_barray(i, B, EPSILON, DISTANCES, data, index_to_change);
		Work *work = scs_init(data, cone, info);

		mat positions = zeros<mat>(i, i);
		vec embedded_point = zeros<vec>(i);
		vector<scs_float> value_to_change(2 * i);
		for (int j = 0; j < (int) perm.size(); j++) {
			// prepare the b array;
			if (!outliner[perm[j]]) {
				for (int k = 0; k <= i - 1; k++) {
					d_low = DISTANCES(B[k], perm[j]) - EPSILON;
					d_high = DISTANCES(B[k], perm[j]) + EPSILON;
					value_to_change[k << 1] = -d_low * d_low;
					value_to_change[(k << 1) + 1] = d_high * d_high;
				}
				update_barray(data, index_to_change, value_to_change);
				scs_solve(work, data, cone, sol, info);
				if (info->statusVal > 0) {
					//cerr << positions << endl;
					//cerr << sol << endl;
					double distance_to_plane = parse_distance(positions,
							embedded_point, sol);
					//cerr << embedded_point << endl;
					H.push_back(make_pair(distance_to_plane, perm[j]));
				} else {
					ui.push_back(perm[j]);
				}
			}
		}
		standard_release(data, cone, sol, work);

		sort(H.begin(), H.end());
		int h = H.size();
		int upper = h < (2 * K) ? h : (2 * K);
		//cerr << sqrt_epsilon << " " << h << " " << H[h - 1].first << " " << H[0].first << endl;
		if (upper != 0 && H[h - 1].first > sqrt_epsilon) {
			srand((unsigned) time(NULL));
			int r = H.size() - 1 - rand() % upper;
			B.push_back(H[r].second);
			for (int k = r; k < h; k++) {
				outliner[H[k].second] = true;
			}
			for (int k = r + 1; k < h; k++) {
				ui.push_back(H[k].second);
			}
			U.push_back(ui);
		} else {
			cout
					<< "First step stopped before finding out d + 1 good points, try to embed in lower dimension."
					<< endl;
		}
	}
}

Embedding second_step(vector<int> &B, vector<vector<int> > &U, int N,
		double EPSILON, double EPSILON2, mat &DISTANCES) {
//	vector<int> IMAP(N);
//	for (int i = 0; i < N; i++) {
//		IMAP[i] = i;
//	}

	vector<int> res_by_dim(B.size() + 1, 0);
	vector<int> optimal_embeddable_set;
	mat optimal_positions;
	int map_back[N];
	bool outliners[N];
	memset(outliners, false, sizeof(outliners));
	outliners[B[0]] = true;

	for (int i = 1; i < (int) B.size(); i++) {

		mat fixed_square_distances = zeros<mat>(i + 1, i + 1);
		sdp_squared_distances_xi_xj(i, EPSILON, B, DISTANCES,
				fixed_square_distances);
		mat fixed_positions;
		embed_fixed_points(fixed_square_distances, fixed_positions);
//		vector<int> Q(N);
//		vector<int> SORTED_B(B.begin(), B.begin() + i + 1);
//		sort(SORTED_B.begin(), SORTED_B.end());
//		auto it = set_difference(IMAP.begin(), IMAP.end(), SORTED_B.begin(),
//				SORTED_B.end(), Q.begin());
//		Q.resize(it - Q.begin());
		outliners[B[i]] = true;
		if ((int) U.size() >= i) {
			for (auto u : U[i - 1]) {
				outliners[u] = true;
			}
		}
		vector<int> embeddable_points;
		mat embeddable_positions = zeros<mat>(i, N);

		int cnt_embeddable_points = 0;

		Data *data;
		Cone *cone;
		Sol *sol;
		Info *info;

//	 	old method, try to use different i points to embed.
#ifdef OLD_METHOD
		vec q_distances = zeros<vec>(i);
		init_xij_a_problem(data, cone, sol, info, i);
#else
//		new method, use all i + 1 points, then project to R^i space.
		vec q_squared_distances = zeros<vec>(i + 1);
		vector<int> index_2change;
		init_xij_a_problem(data, cone, sol, info, i + 1);
		init_xij_a_barray(i + 1, B, EPSILON, fixed_square_distances, data,
				index_2change);
#endif
// -----------------------

		Work *work = scs_init(data, cone, info);
		for (int q = 0; q < N; q++) {
			if (!outliners[q]) {
				vec q_opvec;

// old method
#ifdef OLD_METHOD
				double q_opdist = 10000.0;
				for (int p_kick = 0; p_kick <= i; p_kick++) {
					init_xij_a_barray(i, p_kick, q, B, EPSILON,
							fixed_square_distances, DISTANCES, data);
					scs_solve(work, data, cone, sol, info);
					if (info->statusVal > 0) {
						parse_q_distances(i, sol, q_distances);
						auto qp_pair = embed_new_point(p_kick, fixed_positions,
								fixed_square_distances, q_distances);
						vec qpr;
						double temp_d;
						get_embedding(qp_pair.first, qp_pair.second,
								fixed_positions.col(p_kick), qpr,
								DISTANCES(q, B[p_kick]), temp_d);
						if (temp_d < q_opdist) {
							q_opdist = temp_d;
							q_opvec = qpr;
						}
					}
				}

				if (q_opdist < EPSILON2) {
					embeddable_points.push_back(q);
					map_back[q] = cnt_embeddable_points;
					embeddable_positions.col(cnt_embeddable_points) = q_opvec;
					cnt_embeddable_points++;
				}
#else
// new method
				double d_low, d_high;
				for (int i = 0; i < (int) index_2change.size(); i += 2) {
					d_low = DISTANCES(B[i >> 1], q) - EPSILON2;
					d_high = DISTANCES(B[i >> 1], q) + EPSILON2;
					if (d_low < 0)
					d_low = 0;
					data->b[index_2change[i]] = -d_low * d_low;
					data->b[index_2change[i + 1]] = d_high * d_high;
				}

				scs_solve(work, data, cone, sol, info);
				if (info->statusVal > 0) {
					parse_q_distances(i + 1, sol, q_squared_distances);
//					cerr << q_squared_distances << endl;
//					cerr << fixed_positions << endl;
					q_opvec = embed_new_point(fixed_positions,
							fixed_square_distances, q_squared_distances);
//					cerr << q_opvec << endl;

					bool flag = true;
					for (int r = 0; r < i + 1; r++) {
						double td = pairwise_distance(fixed_positions.col(r),
								q_opvec);

						if (abs(td - DISTANCES(q, B[r])) > EPSILON2) {
							flag = false;
							break;
						}
					}

					if (flag) {
						embeddable_points.push_back(q);
						map_back[q] = cnt_embeddable_points;
						embeddable_positions.col(cnt_embeddable_points) =
						q_opvec;
						cnt_embeddable_points++;
					}
				}
#endif

// ----------------

			}
		}

		standard_release(data, cone, sol, work);

		vector<pair<int, int> > conflict_edges = get_conflict_edges(
				embeddable_positions, embeddable_points, DISTANCES,
				2 * EPSILON2 + 2 * EPSILON);
		vector<int> vertex_cover = two_approximate_vertex_cover(conflict_edges);

		vector<int> final_points(N);
		sort(vertex_cover.begin(), vertex_cover.end());
		auto it = set_difference(embeddable_points.begin(),
				embeddable_points.end(), vertex_cover.begin(),
				vertex_cover.end(), final_points.begin());
		final_points.resize(it - final_points.begin());
		int es = final_points.size();
		for (int j = 0; j <= i; j++) {
			final_points.push_back(B[j]);
		}

		cout << i << " dimensional space " << final_points.size()
				<< " points embeddable." << endl;

		res_by_dim[i] = max(res_by_dim[i], (int) final_points.size());

		if (final_points.size() > optimal_embeddable_set.size()) {
			optimal_embeddable_set = final_points;
			mat temp_positions = zeros<mat>(i, es);
			for (int j = 0; j < es; j++) {
				temp_positions.col(j) = embeddable_positions.col(
						map_back[final_points[j]]);
			}
			optimal_positions = join_rows(temp_positions, fixed_positions);
		}
	}
	return Embedding(res_by_dim, optimal_embeddable_set, optimal_positions);
}

void norm_distance(Embedding &ret, vector<double> &norms, mat &DISTANCES) {
	int size = ret.size();
	double temp;
	double norm2 = 0.0;
	double norm1 = 0.0;
	double norm_max = 0.0;
//	int ii, ij;
	for (int i = 0; i < size; i++) {
		for (int j = i + 1; j < size; j++) {
			temp = pairwise_distance(ret.embedded_positions.col(i),
					ret.embedded_positions.col(j));
			temp = temp
					- DISTANCES(ret.embedded_point_indexes[i],
							ret.embedded_point_indexes[j]);
			norm2 += temp * temp;
			norm1 += abs(temp);
			if (abs(temp) > norm_max) {
				norm_max = abs(temp);
//				ii = ret.embedded_point_indexes[i];
//				ij = ret.embedded_point_indexes[j];
			}
		}
	}
//	cout << "max: (" << ii << ", " << ij << ") " << norm_max << "." << endl;
	double avg_norm1 = norm1 / (size * size);
	norm2 = sqrt(norm2);
	double avg_norm2 = norm2 / size;
	norms.push_back(norm_max);
	norms.push_back(norm1);
	norms.push_back(avg_norm1);
	norms.push_back(norm2);
	norms.push_back(avg_norm2);
}

