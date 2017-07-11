/*
 * model.h
 *
 */

#ifndef MODEL_H_
#define MODEL_H_

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "corpus.h"

typedef struct vlist{
	int len;
	int v;
	struct vlist* next;
}vlist;

typedef struct Model {
	int I; //number of classes of the data
	int v; // total vocabulary length
	int *W_I;// sub-vocabulary length for each class
	int *K_I; //number of topics for each class data
	double *a; // PDP parameters a for each class, K
	double *b; // PDP parameters b for each class, K
	double gamma_sum;
	gsl_vector** alpha;
	double* alpha_sum;
	gsl_vector* gamma;
	gsl_matrix* phi;
	gsl_matrix** phi_i;
	double*** P; // projection matrices, sparse representation
	int*** P_ind; // indexes for non-zero elements
	int** P_size; // index sizes for each row of the projection matrices
	gsl_vector** P_k; // row sums of the projection matrices
	gsl_matrix* rest_P;

	double **P_all;
	int **P_idx_all;
	int * P_size_all;
	int maxV;	// max P_size
	int ****q;
	vlist**** vL;
} Model;

/*
 * count tables
 */
typedef struct Top_cts {
	int** m; // K_{i} * W_i
	int* M; //K_i
	int** t; // (K_{i}+1) * W_i
	int* T_W; //W_i + 1
	int* T_K; // K_{i}+1
} Top_cts;

typedef struct Cts {
	int I;
	double* sum_T;	// \sum_i\sum_w t_{ikw}, K
	double* B; //K
	double** T_KQ; //K * V
	int*** n; // d_i * K_i
	int** N; // d_i
	Top_cts* top_cts; // I
} Cts;

/*
 * topic assignments
 */
typedef struct assignment {
	int*** topic_ass; //I * D_i * L_{D_i}
} Assignment;

/*
 * estimator
 */
typedef struct Estimator{
	gsl_matrix** mu;
	gsl_matrix* phi;
	gsl_matrix** phi_i;
} Estimator;
/*
 * functions
 */
Model* new_model(int I, int v, int *K_i, int * W_i, int inf);
Cts* new_cts(int I, int * K_i, int *W_i, int v, Corpus* c);
Assignment* new_assignment(Corpus* c);

Model* random_init(int I, int *K_i, int *W_i, int v, double *alpha, double gamma_sum,
		double *a, double *b, int** l2g_vocabulary, char* dir_p);//double***p, int*** p_ind, int** p_size);
Model* read_model(char* file);
void save_model(Model* model, Corpus* c, Cts* cts, gsl_matrix* mphi, gsl_matrix** mphi_i, char* file);
/*
 * estimator for computing topic proportions
 */
Estimator* new_estimator(Corpus* c, int I, int *K_i, int *W_i, int v, int do_mu, int do_phi);
void printf_estimator(Estimator* est, Corpus* c, char* file, int do_mu, int do_phi);
void save_estimator(Estimator* est, Corpus* c, char* file, int do_mu, int do_phi);
Estimator* read_estimator(Corpus* c, int v, int *K_i, int *W_i, char* file, int do_mu, int do_phi);
void free_estimator(Estimator* est, Corpus* c, int do_mu, int do_phi);

void calPMI(char* dir, char** model_dir, int numM, int topk, char* pmifile);

//void save_table(Cts* cts, Assignment* ass, Corpus* c, char* file);
//void read_table(Cts* cts, Assignment* ass, Corpus* c, char* file);

Assignment* read_topic_assignmnet(Corpus* c, char* file);
void save_topic_assignmnet(Corpus* c, Assignment* ass, char* file);
void print_top_words(int num_words, int begin_k, int end_k, gsl_matrix* M, Cts* cts, vocabulary* v,
		char* file);
void print_local_top_words(int* num_words, int begin_k, int end_k, gsl_matrix**M_i,
		Cts* cts, Corpus* c, vocabulary* v, char* file);
void print_word_asso(Model* model, vocabulary* voc, char* filename);

void free_assignment(Assignment*, Corpus* c);
void free_cts(Cts*, Corpus* c, int K);
void free_model(Model*);
void free_all(Corpus*, Model*, Cts*, Assignment*, vocabulary*);

void save_P(char* file, Model* model, int i);
void load_P(char* file, Model* model, int i);
void load_P_all(char* file, Model* model);
void save_P_ind(char* file, Model* model, int i);
void load_P_ind(char* file, Model* model, int i);
void load_P_ind_all(char* file, Model* model);
void save_P_size(char* file, Model* model);
void load_P_size(char* file, Model* model);
void load_P_size_all(char* file, Model* model);

void save_customer_table(Cts* cts, int I, int K, int v, char* file);

void deleteList(vlist* list);
vlist* addListNode(vlist* list, int v);
vlist* removeListNode(vlist* list, int i);
int getNodeV(vlist* list, int i);

#endif /* MODEL_H_ */
