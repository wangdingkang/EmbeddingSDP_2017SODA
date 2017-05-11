/*
 * sdp.cpp
 *
 *  Created on: Apr 25, 2016
 *      Author: wdk
 *      The meaning of a whole bunch of parameter explained on scs/github.
 */

#include "sdp.hpp"

void parse_res(int D, Sol* &sol, mat &save_to) {
	for (int i = 0, k = 0; i < D; i++) {
		for (int j = i; j < D; j++) {
			save_to(i, j) = save_to(j, i) = sol->x[k++];
		}
	}
}

void parse_q_distances(int D, Sol* &sol, vec &save_to) {
	mat temp = zeros<mat>(D, D);
	for (int i = 0, k = 0; i < D; i++) {
		for (int j = i; j < D; j++) {
			temp(i, j) = temp(j, i) = sol->x[k++];
		}
	}
	for (int i = 0; i < D - 1; i++) {
		save_to(i + 1) = temp(i, i) + temp(D - 1, D - 1) - 2 * temp(D - 1, i);
	}
	save_to(0) = temp(D - 1, D - 1);
}

void standard_malloc(Data* &data, Cone* &cone, Sol* &sol, Info* &info,
		AMatrix* &A, scs_int m, scs_int n, scs_int nnz, int ep, int qsize,
		int ssize, vector<scs_int> s_values, int l, int f) {
	data = (Data *) scs_calloc(1, sizeof(Data));
	sol = (Sol *) scs_calloc(1, sizeof(Sol));
	cone = (Cone *) scs_calloc(1, sizeof(Cone));
	info = (Info *) scs_calloc(1, sizeof(Info));
	A = data->A = (AMatrix *) scs_calloc(1, sizeof(AMatrix));
	data->stgs = (Settings *) scs_calloc(1, sizeof(Settings));
	setDefaultSettings(data);
	data->stgs->verbose = 0;

	data->m = m;
	data->n = n;
	A->m = m;
	A->n = n;

	scs_float *b = data->b = (scs_float *) scs_calloc(m, sizeof(scs_float));
	scs_float *c = data->c = (scs_float *) scs_calloc(n, sizeof(scs_float));
	for (int i = 0; i < n; i++) {
		c[i] = 0.0;
	}
	for (int i = 0; i < m; i++) {
		b[i] = 0.0;
	}

	A->i = (scs_int *) scs_calloc(nnz, sizeof(scs_int));
	A->p = (scs_int *) scs_calloc((n + 1), sizeof(scs_int));
	A->x = (scs_float *) scs_calloc(nnz, sizeof(scs_float));
	A->p[0] = 0;

	cone->ep = ep;
	cone->qsize = qsize;
	cone->s = (scs_int *) scs_calloc(ssize, sizeof(scs_int));
	for (int i = 0; i < ssize; i++) {
		cone->s[i] = s_values[i];
	}
	cone->ssize = ssize;
	cone->l = l;
	cone->f = f;

	sol->x = (scs_float *) scs_calloc(n, sizeof(scs_float));
	sol->y = (scs_float *) scs_calloc(m, sizeof(scs_float));
	sol->s = (scs_float *) scs_calloc(m, sizeof(scs_float));
}

/*
 D * D semidefinite matrix
 l, length of inequalities = D * D + D.
 f, length of equalities, 0
 ssize = 1, s[0] = D.
 m = D * D + D + D * (D + 1) / 2, n = D * (D + 1) / 2; nnz = 3 * D * D - D + D * (D + 1) / 2;
 set the elements in matrix A is a pretty tedious and frustrating work.
 lucky we only need to initialize it once.
 */
void init_xij_problem(Data* &data, Cone* &cone, Sol* &sol, Info* &info, int D) {

	AMatrix *A;
	int m = D * D + D + D * (D + 1) / 2;
	int n = D * (D + 1) / 2;
	int nnz = 3 * D * D - D + D * (D + 1) / 2;
	int ep = 0, qsize = 0, l = D * D + D, f = 0, ssize = 1;
	vector<scs_int> s_values(ssize);
	s_values[0] = D;
	standard_malloc(data, cone, sol, info, A, m, n, nnz, ep, qsize, ssize,
			s_values, l, f);

	int get_new_cols[D][D];
	for (int i = 0, k = 0; i < D; i++) {
		for (int j = i; j < D; j++) {
			get_new_cols[j][i] = k++;
		}
	}

	vector<CCF_Element> elements;
	int t_row = 0;
	for (int i = 0; i < D; i++) {
		elements.push_back(CCF_Element(get_new_cols[i][i], t_row++, -1.0));
		elements.push_back(CCF_Element(get_new_cols[i][i], t_row++, 1.0));
	}
	for (int i = 0; i < D; i++) {
		for (int j = i + 1; j < D; j++) {
			elements.push_back(CCF_Element(get_new_cols[i][i], t_row, -1.0));
			elements.push_back(CCF_Element(get_new_cols[j][j], t_row, -1.0));
			elements.push_back(CCF_Element(get_new_cols[j][i], t_row, 2.0));
			t_row++;
			elements.push_back(CCF_Element(get_new_cols[i][i], t_row, 1.0));
			elements.push_back(CCF_Element(get_new_cols[j][j], t_row, 1.0));
			elements.push_back(CCF_Element(get_new_cols[j][i], t_row, -2.0));
			t_row++;
		}
	}

	for (int i = 0, k = 0; i < D; i++) {
		for (int j = i; j < D; j++) {
			elements.push_back(
					CCF_Element(k++, t_row++, (i == j ? -1.0 : -SQRT2)));
		}
	}

	sort(elements.begin(), elements.end());

	int t_col = 1;
	for (int j = 0; j < (int) elements.size(); j++) {
		A->x[j] = elements[j].v;
		A->i[j] = elements[j].row;
		if (j != 0 && elements[j].col != elements[j - 1].col)
			A->p[t_col++] = j;
	}
	A->p[t_col] = elements.size();

}

/*
 *   embed new point a in D-dim space.
 l, length of inequalities, D * 2.
 f, length of equalities, (D - 1) * D / 2
 ssize = 1, s[0] = D;
 m = D * D + 2 * D; D * (D + 1) / 2; nnz = 6 * D - 4 + D * D;
 annoying too...
 */
void init_xij_a_problem(Data* &data, Cone* &cone, Sol* &sol, Info * &info,
		int D) {
	AMatrix *A;
	int m = D * D + 2 * D;
	int n = D * (D + 1) / 2;
	int nnz = 6 * D - 4 + D * D;
	int f = (D - 1) * D / 2, l = D * 2, ssize = 1, ep = 0, qsize = 0;
	vector<scs_int> s_values(ssize);
	s_values[0] = D;
	standard_malloc(data, cone, sol, info, A, m, n, nnz, ep, qsize, ssize,
			s_values, l, f);

	int get_new_cols[D][D];
	for (int i = 0, k = 0; i < D; i++) {
		for (int j = i; j < D; j++) {
			get_new_cols[j][i] = k++;
		}
	}

	vector<CCF_Element> elements;
	// add equalities first X[i, j] = constant.
	int t_row = 0;
	for (int i = 0; i < D - 1; i++) {
		for (int j = i; j < D - 1; j++) {
			elements.push_back(CCF_Element(get_new_cols[j][i], t_row++, 1.0));
		}
	}

	// add constraints about X[D - 1, D - 1].
	elements.push_back(CCF_Element(get_new_cols[D - 1][D - 1], t_row++, -1.0));
	elements.push_back(CCF_Element(get_new_cols[D - 1][D - 1], t_row++, 1.0));

	// add inequalities d- <= X[D - 1, i] <= d+ (0 <= i < D - 1)
	for (int i = 0; i < D - 1; i++) {
		elements.push_back(
				CCF_Element(get_new_cols[D - 1][D - 1], t_row, -1.0));
		elements.push_back(CCF_Element(get_new_cols[i][i], t_row, -1.0));
		elements.push_back(CCF_Element(get_new_cols[D - 1][i], t_row, 2.0));
		t_row++;
		elements.push_back(CCF_Element(get_new_cols[D - 1][D - 1], t_row, 1.0));
		elements.push_back(CCF_Element(get_new_cols[i][i], t_row, 1.0));
		elements.push_back(CCF_Element(get_new_cols[D - 1][i], t_row, -2.0));
		t_row++;
	}

	// add SDP constraints.
	for (int i = 0, k = 0; i < D; i++) {
		for (int j = i; j < D; j++) {
			elements.push_back(
					CCF_Element(k++, t_row++, (i == j ? -1.0 : -SQRT2)));
		}
	}

	sort(elements.begin(), elements.end());

	int t_col = 1;
	for (int j = 0; j < (int) elements.size(); j++) {
		A->x[j] = elements[j].v;
		A->i[j] = elements[j].row;
		if (j != 0 && elements[j].col != elements[j - 1].col)
			A->p[t_col++] = j;
	}
	A->p[t_col] = elements.size();

}

// initialize the b array for step one, only need once.
// hard to implement....
void init_xij_barray(int D, vector<int> &B, double EPSILON, mat &DISTANCES,
		Data* &data, vector<int> &index_to_change) {

	int t_index = 0;
	for (int i = 1; i < D; i++) {
		double d_low = DISTANCES(B[i], B[0]) - EPSILON;
		double d_high = DISTANCES(B[i], B[0]) + EPSILON;
		data->b[t_index++] = -d_low * d_low;
		data->b[t_index++] = d_high * d_high;
	}

	index_to_change.push_back(t_index++);
	index_to_change.push_back(t_index++);

	for (int i = 1; i <= D; i++) {
		for (int j = i + 1; j <= D; j++) {
			if (j != D) {
				double d_low = DISTANCES(B[i], B[j]) - EPSILON;
				double d_high = DISTANCES(B[i], B[j]) + EPSILON;
				// if the pertubated distance is less than 0, then we set
				// it to 0, since there can't be negative distance
				if (d_low < 0)
					d_low = 0;
				data->b[t_index++] = -d_low * d_low;
				data->b[t_index++] = d_high * d_high;
			} else {
				index_to_change.push_back(t_index++);
				index_to_change.push_back(t_index++);
			}
		}
	}
	// cout << "init_xij_barray finished." << endl;
}

// initialize the b array for step two, only need once.
void init_xij_barray(int D, vector<int> &B, double EPSILON, mat &DISTANCES,
		Data* &data) {
	int t_index = 0;
	for (int i = 1; i <= D; i++) {
		double d_low = DISTANCES(B[i], B[0]) - EPSILON;
		double d_high = DISTANCES(B[i], B[0]) + EPSILON;
		// if the pertubated distance is less than 0, then we set
		// it to 0, since there can't be negative distance
		if (d_low < 0)
			d_low = 0;
		data->b[t_index++] = -d_low * d_low;
		data->b[t_index++] = d_high * d_high;
	}
	for (int i = 1; i <= D; i++) {
		for (int j = i + 1; j <= D; j++) {
			double d_low = DISTANCES(B[i], B[j]) - EPSILON;
			double d_high = DISTANCES(B[i], B[j]) + EPSILON;
			// if the pertubated distance is less than 0, then we set
			// it to 0, since there can't be negative distance
			if (d_low < 0)
				d_low = 0;
			data->b[t_index++] = -d_low * d_low;
			data->b[t_index++] = d_high * d_high;
		}
	}
}

// initialzie the b array for embed new point problem.
void init_xij_a_barray(int D, vector<int> &B, double EPSILON,
		mat &fixed_squared_distances, Data* &data,
		vector<int> &index_to_change) {
	int t_index = 0;
	for (int i = 0; i < D - 1; i++) {
		for (int j = i; j < D - 1; j++) {
			if (i == j) {
				data->b[t_index++] = fixed_squared_distances(i + 1, 0);
			} else {
				data->b[t_index++] = (fixed_squared_distances(i + 1, 0)
						+ fixed_squared_distances(j + 1, 0)
						- fixed_squared_distances(i + 1, j + 1)) / 2;
			}
		}
	}

	// add constraints about X[D - 1, D - 1].
	// add inequalities d- <= X[D - 1, i] <= d+ (0 <= i < D - 1)
	for (int i = 0; i < 2 * D; i++) {
		index_to_change.push_back(t_index++);
	}

}

// initialzie the b array for embed new point problem.
void init_xij_a_barray(int D, int p_kick, int q, vector<int> &B, double EPSILON,
		mat &fixed_squared_distances, mat &DISTANCES, Data* &data) {
	int t_index = 0;
	int pivot = (p_kick == 0 ? 1 : 0);
	for (int i = pivot + 1; i <= D; i++) {
		for (int j = i; j <= D; j++) {
			if (i != p_kick && j != p_kick) {
				if (i == j) {
					data->b[t_index++] = fixed_squared_distances(i, pivot);
				} else {
					data->b[t_index++] = (fixed_squared_distances(i, pivot)
							+ fixed_squared_distances(j, pivot)
							- fixed_squared_distances(i, j)) / 2;
				}
			}
		}
	}
	double d_low, d_high;
	for (int i = 0; i <= D; i++) {
		if (i != p_kick) {
			d_low = DISTANCES(B[i], q) - EPSILON;
			d_high = DISTANCES(B[i], q) + EPSILON;
			// if the pertubated distance is less than 0, then we set
			// it to 0, since there can't be negative distance
			if (d_low < 0)
				d_low = 0;
			data->b[t_index++] = -d_low * d_low;
			data->b[t_index++] = d_high * d_high;
		}
	}
}

// N * N semidefinite matrix, B[1 ... N];
void sdp_squared_distances_xi_xj(int N, double EPSILON, vector<int> &B,
		mat &DISTANCES, mat &fixed_square_distances) {
	Data *data;
	Cone *cone;
	Sol *sol;
	Info *info;
	init_xij_problem(data, cone, sol, info, N);
	init_xij_barray(N, B, EPSILON, DISTANCES, data);
	scs(data, cone, sol, info);
	mat ret = zeros<mat>(N, N);
	parse_res(N, sol, ret);
	for (int i = 1; i <= N; i++) {
		fixed_square_distances(i, 0) = fixed_square_distances(0, i) = ret(i - 1,
				i - 1);
	}
	for (int i = 1; i <= N; i++) {
		for (int j = i + 1; j <= N; j++) {
			fixed_square_distances(i, j) = fixed_square_distances(j, i) = ret(
					i - 1, i - 1) + ret(j - 1, j - 1) - 2 * ret(i - 1, j - 1);
		}
	}
}

// Free memory
void standard_release(Data* &data, Cone* &cone, Sol* &sol, Work* &work) {
	scs_finish(work);
	freeData(data, cone);
	freeSol(sol);
}

void standard_release(Data* &data, Cone* &cone, Sol* &sol) {
	freeData(data, cone);
	freeSol(sol);
}

