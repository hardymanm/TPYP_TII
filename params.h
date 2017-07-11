/*
 * params.h
 *
 */

#ifndef PARAMS_H_
#define PARAMS_H_

typedef struct dpyp_params {
	int I;
	int *K_i;
	int *W_i;
	double *a;
	double *b;
	double *alpha;
	double gamma;
	double tr_percent;
	int gibbs_max_iter;
	int burn_in;
	int sampling_lag;
	int *top_words;
	char** train_data;
} dpyp_params;

void new_params(int I, int K, dpyp_params* params);
void free_params(dpyp_params* params);
void read_test_params(char*);
void read_pmi_params(char*);
void read_params(char*);
void write_params(char*);
void print_params();
void default_params();
void set_paras(int I, int *K_i, int* W_i, double *a, double *b, double *alpha, double gamma,
		double tr_per, int gibbs_max_iter, int burn_in, int sampling_lag, int *top_words);
#endif /* PARAMS_H_ */
