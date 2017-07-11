/*
 * sample_alpha.c
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "util.h"
#include "seqTM.h"
#include "arms.h"

#define ALPHA_MIN 0.001
#define ALPHA_MAX 1

extern gsl_rng* glob_r;

static int ti, tK;
static Corpus* tc;
static Cts* tcts;

static double alphaterms(double alpha, void* data)
{
        double val;
	int d, k;
	double lga = gsl_sf_lngamma(alpha);
	double lgaK = gsl_sf_lngamma(alpha*tK);

	/*
	 *   use a Gamma prior with mean 1/K (or ALPHA_MIN if greater)
	 */
	val = - tK*alpha;
	if ( tK*ALPHA_MIN<1 ) 
	  val = - alpha / ALPHA_MIN;
	for (d=0; d<tc->ndocs[ti]; d++) {
	  for (k=0; k<tK; k++) {
	    val += gsl_sf_lngamma(alpha + tcts->n[ti][d][k]) - lga;
	  } 
	  val -=  gsl_sf_lngamma(alpha*tK + tcts->N[ti][d]) - lgaK;
	}
	printf("ARMS ALPHA:  %lf -> %lf\n", alpha, val);
	return val;
}

void sample_alpha(Cts* cts, Model* model, Corpus* c)
{
  int i;

	double lv = ALPHA_MIN;
	double rv = ALPHA_MAX;

	tc = c;
	tcts = cts;

	printf("<<<<<<<<Begin reject sampling for alpha......\n");

	for(i = 0; i < cts->I; i++){
	  double newalpha[5];
	  double alpha = gsl_vector_get(model->alpha[i], 0);
	  ti = i;
	  tK = model->K_I[i];
	  arms_simple(4, &lv, &rv, alphaterms, NULL, 0, &alpha, newalpha);  
	  if ( newalpha[0]<lv ) {
	    printf("ARMS ALPHA OUT: %lf, reset to %lf\n", newalpha[0], lv);
	    newalpha[0] = lv;
	  } else if ( newalpha[0]>rv ) {
	    printf("ARMS ALPHA OUT: %lf, reset to %lf\n", newalpha[0], rv);
	    newalpha[0] = rv;
	  } else {
	    printf("ARMS ALPHA OUT: %lf\n", newalpha[0]);
	  }
	  gsl_vector_set_all(model->alpha[i], newalpha[0]);
	  model->alpha_sum[i] = newalpha[0] * tK;
	}  
	tc = NULL;
	tcts = NULL;
}

