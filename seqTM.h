/*
 * seqTM.h
 *
 */

#ifndef SEQTM_H_
#define SEQTM_H_

#include "corpus.h"
#include "model.h"

void estimate(Corpus* c, Corpus* c_te, vocabulary* v, char* root, char* data_dir, int init, char* z_root);

double spyp_model_lglihood(Corpus* c, Cts* cts, Model* model);
double spyp_marg_lglihood(Corpus* c, Cts* cts, Assignment* ass, Model* model, int inf);
double spyp_marg_lglihood1(Corpus* c, Cts* cts, Model* model, gsl_matrix** mu, int train);

void gamma_opt(Cts* cts, Model* model, Corpus* c);
void gamma_opt_newton(Cts* cts, Model* model, Corpus* c);
void alpha_opt(Cts* cts, Model* model, int* nDoc, int K);
void alpha_opt_newton(Cts* cts, Model* model, int* nDoc, int K);
void sample_alpha(Cts* cts, Model* model, Corpus* c);

void sample_b(int rep, Cts* cts, Model* model, Corpus* c);
void sample_b_all(int rep, Cts* cts, Model* model, Corpus* c);

void sample_ab(int rep, Cts* cts, Model* model, Corpus* c);

#endif /* SEQTM_H_ */
