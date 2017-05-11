//============================================================================
// Name        : Embeddings.cpp
// Author      : Dingkang Wang
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "subroutines.hpp"
#include "IO.hpp"

//#define INPUT_MNIST_FILE
#define INPUT_DISTANCE_FILE "./dataset/POLY_KER_swiss_roll_normal_5000.txt" //FACE_PCA_norm2"
// "./dataset/PERTURBATED_DIGIT5_HAM_60k_RATIO=4.000000.txt"

#define OUTPUT_POSITION_PREFIX "./output/POSITIONS_POLY_KER_swiss_roll_normal_5000"
#define OUTPUT_MISCELL_PREFIX "./output/MISCELLANEOUS_POLY_KER_swiss_roll_normal_5000"

int N, D, K;
double EPSILON;
set<int> real_outliners;
double rate_a, rate_b;
vector<double> rate_as, rate_bs;
mat DISTANCES;
vector<vector<int> > ret_sizes;

int main() {

	ios_base::sync_with_stdio(false);

	time_t temp;
	time_t start = time(0);

	temp = time(0);

	cout << "Catch you." << endl;
	read_file(N, D, K, EPSILON, DISTANCES, INPUT_DISTANCE_FILE);

	cout << "Input file read, take " << difftime(time(0), temp) << " seconds."
			<< endl;
//	cerr << log2(N) << endl;
	int iterations = (1 << (D));
//	for (EPSILON = 0.05; EPSILON <= 0.10; EPSILON += 0.01) {
//	for (double EPSILON2 = 2 * EPSILON; EPSILON2 <= sqrt(EPSILON);
//			EPSILON2 += 0.02) { // (sqrt(EPSILON) - EPSILON) / 10)
	// EPSILON = 0.02;
	K = 100;
	EPSILON = 0.01;
	double EPSILON2 = 0.05;

	Embedding ret;

	cout << N << " points, with " << K << " outliners, and epsilon = "
			<< EPSILON << ", epsilon2 = " << EPSILON2 << " in " << D
			<< " dimensional space." << endl;
	cout << "Need " << iterations << " iterations." << endl;

	for (int i = 1; i <= iterations; i++) {
		cout << "Iteration " << i << "." << endl;
		vector<int> random_permutation;
		vector<int> B;
		vector<vector<int> > U;
		gen_random_perm(N, B, random_permutation);

		temp = time(0);
		first_step(random_permutation, B, U, EPSILON, K, D, DISTANCES);
		cout << "Step 1 finished, take " << difftime(time(0), temp)
				<< " seconds." << endl;

		temp = time(0);
		Embedding temp_ret = second_step(B, U, N, EPSILON, EPSILON2, DISTANCES);
		cout << "Step 2 finished, take " << difftime(time(0), temp)
				<< " seconds." << endl;
		vector<double> norms;
		norm_distance(temp_ret, norms, DISTANCES);
//			cerr << temp_ret.size() << endl;
//			for(auto noo : norms) {
//				cerr << noo << " ";
//			}
//			cerr << endl;
		if (temp_ret.size() > ret.size()) {
			ret = temp_ret;
		}
	}
	cerr << ret.size() << endl;
	if (ret.size() > D) {
		vector<double> norms;
		norm_distance(ret, norms, DISTANCES);

		ret_sizes.push_back(ret.embedded_sizes);
		write_miscellaneous(N, norms, ret,
				get_output_filename(OUTPUT_MISCELL_PREFIX, N, K, EPSILON,
						EPSILON2, ".txt"));
		write_positions(ret,
				get_output_filename(OUTPUT_POSITION_PREFIX, N, K, EPSILON,
						EPSILON2, ".csv"));

	} else {
		cout << "Fail to find out a good embedding." << endl;
	}

	cout << "Spend " << difftime(time(0), start) << " seconds." << endl;

//		}
//	}
	cout << INPUT_DISTANCE_FILE << endl;
	for (auto s : ret_sizes) {
		for (int i = 1; i <= D; i++) {
			cout << s[i];
			if (i != D)
				cout << ", ";
		}
		cout << endl;
	}
	return 0;
}

