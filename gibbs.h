/*
 * gibbs.h
 *
 */

#ifndef GIBBS_H_
#define GIBBS_H_

#include <time.h>
#include <math.h>
#include <assert.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>

#include "model.h"
#include "corpus.h"
#include "seqTM.h"
#include "util.h"

void init_state_gibbs(Model* model, Cts* cts, Assignment* ass, Corpus* c, int est_inf, int init_z);

void spyp_e_gibbs(Model* model, Cts* cts, Assignment* ass, Corpus* c);
int remove_old_topic(int i, int d, int w, int k, Cts* cts, Model* model);

void add_new_topic (int i, int d, int w, int v, int k, Cts* cts, Model* model);

void held_out_gibbs(Model* model, Cts* cts, Assignment* ass, Corpus* c);

double left2right_predict(Model* model, Cts* cts, Assignment* ass, Corpus* c, int nSamples, int burnin);

void compute_phi(Cts* cts, Model* model, Corpus* c);
void compute_mean_phi(gsl_matrix* mphi, gsl_matrix** mphi_i, Model* model, int n);
void compute_topic_proportion(Estimator* estimator, Cts* cts, Model* model, Corpus* c);

void compute_mu(Cts* cts, Model* model, Corpus* c, gsl_matrix** mu);

#endif /* GIBBS_H_ */
