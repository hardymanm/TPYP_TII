/*
 * sample_b.c
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
#include "model.h"
#include "arms.h"

#define B_MIN 0.01
#define B_MAX 5000 //5000000

extern gsl_rng* glob_r;

static double a;
static double lgq;
static int k_id;
static Corpus* tc;
static Cts* tcts;

#define arms_b


#ifdef arms_b
static double bterms_1(double b, void* data)
{
	double val;
	int j;

	val = lgq*b;
	double lba =  gsl_sf_lngamma(b/a);

	for(j = 0; j < tc->I; j++){
		double T = tcts->top_cts[j].T_K[k_id];
		assert(b/a + T > 0);
		val+= gsl_sf_lngamma(b/a+T) - lba;
	}

	return val;
}
#endif

void sample_b(int rep, Cts* cts, Model* model, Corpus* c)
{
	int i, j;

	double lv = B_MIN;
	double rv = B_MAX;

	tc = c;
	tcts = cts;

	//printf("<<<<<<<<Begin reject sampling for b......\n");

	while((rep--) > 0){
		for(i = 0; i < model->K_I[0]; i++){
			double newb[10];
			double b = cts->B[i], tot_t = 0;
			a = model->a[i];
			lgq = 0;
			k_id = i;
			for(j = 0; j < cts->I; j++){
				double tot_cust = cts->top_cts[j].M[i];
				tot_t += cts->top_cts[j].T_K[i];
				if(tot_cust > 0){
					double q = 0;
					int try = 20;
					while(q <= 0 && try >= 0){
						q = gsl_ran_beta(glob_r, b, tot_cust);
						try--;
					}
					if(q <= 0){
						fprintf(stderr, "q in sample_b(b = %f, tot_cust = %f) went zero!!!\n", b, tot_cust);
						exit(1);
					}
					if(a > 0)
						lgq += gsl_sf_log(q);
					else
						lgq -= gsl_sf_log(q);
				}
			}

			if(a > 0){
#ifdef arms_b
				arms_simple (8, &lv, &rv, bterms_1, NULL, 0, &b, newb);
#endif
				cts->B[i] = newb[0];
				model->b[i] = newb[0];
			}else{
				if ( tot_t > 400 ) {
				  //   gamma is near enough a Gaussian with tiny std dev.
				  newb[0] = tot_t / lgq + gsl_ran_gaussian(glob_r, sqrt(tot_t) / lgq);
				  if ( newb[0] <= 0 )
					  newb[0] = gsl_ran_gamma(glob_r, tot_t, 1.0 / lgq);
				} else
				  // GSL has beta parameter as inverse to some others
					newb[0] = gsl_ran_gamma(glob_r, tot_t, 1.0 / lgq);
				//if(newb[0] < B_MIN)
				//	newb[0] = B_MIN;
				//if(newb[0] > B_MAX)
				//	newb[0] = B_MAX;
				cts->B[i] = newb[0];
				model->b[i] = newb[0];
			}
#ifndef SAB
			if(rep == 1)// && (i%10 == 0))
				printf("b[%d] = %lf, ", i, newb[0]);
#endif
		}
	}
#ifndef SAB
	printf("\n");
#endif
	tc = NULL;
	tcts = NULL;
}

void sample_ab(int rep, Cts* cts, Model* model, Corpus* c)
{
	int i, j, t, w;//, ave_c, k;

	//printf("<<<<<<<<Begin reject sampling for b......\n");

	while((rep--) > 0){
		for(i = 0; i < model->K_I[0]; i++){
			double b = cts->B[i], log_x = 0, y = 0, y_1 = 0, z = 0, z_1 = 0;

			double a = model->a[i];
			for(j = 0; j < cts->I; j++){
				if(cts->top_cts[j].T_K[i] > 1){
					log_x += gsl_sf_log(gsl_ran_beta(glob_r, b + 1, cts->top_cts[j].M[i] - 1));
					for(t = 1; t < cts->top_cts[j].T_K[i]; ++t){
						double yy = ((next_uniform() >= b / (b + a * t))?1:0);
						y += yy;
						y_1 += (1 - yy);
					}
				}
				for(w = 0; w < model->W_I[j]; ++w){
					///////////////////////////////////////////////////////////////////
					/*ave_c = (int)(cts->top_cts[j].m[i][w] / cts->top_cts[j].t[i][w] + 0.5);
					if(ave_c > 1){
						for(k = 0; k < cts->top_cts[j].t[i][w]; ++k){
							for(t = 1; t < ave_c; ++t){
								double zz = ((next_uniform() >= (t - 1) / (t - a))?1:0);
								z += zz;
								z_1 += (1 - zz);
							}
						}
					}else{
						for(t = 1; t < cts->top_cts[j].m[i][w] - cts->top_cts[j].t[i][w]; ++t){
							double zz = ((next_uniform() >= (t - 1) / (t - a))?1:0);
							z += zz;
							z_1 += (1 - zz);
						}
					}*/
					//////////////////////////////////////////////////////////////////////////
					if(cts->top_cts[j].m[i][w] > 1){
						for(t = 1; t < cts->top_cts[j].m[i][w] - cts->top_cts[j].t[i][w]; ++t){
							double zz = ((next_uniform() >= (t - 1) / (t - a))?1:0);
							z += zz;
							z_1 += (1 - zz);
						}
					}
					//////////////////////////////////////////////////////////////////////////
				}
			}
			a = gsl_ran_beta(glob_r, 1 + y_1, 1 + z_1);
			model->a[i] = a;
			sample_b(1, cts, model, c);
			//b = gsl_ran_gamma(glob_r, 1 + y, 1 / (1 - log_x));
			//cts->B[i] = b; model->b[i] = b;

			if(rep == 1)
				printf("a[%d] = %lf, b[%d] = %lf, ", i, a, i, model->b[i]);
		}
	}
	printf("\n");
}

/***
* all b the same
*/

int K;
static double bterms_all(double b, void* data)
{
	double val;
	int j, k;

	val = lgq*b;
	double lba =  gsl_sf_lngamma(b/a);

	for(k = 0; k < K; k++){
		for(j = 0; j < tc->I; j++){
			double T = tcts->top_cts[j].T_K[k];
			assert(b/a + T > 0);
			val+= gsl_sf_lngamma(b/a+T) - lba;
		}
	}
	return val;
}

void sample_b_all(int rep, Cts* cts, Model* model, Corpus* c)
{
	int i, j;

	double lv = B_MIN;
	double rv = B_MAX;

	tc = c;
	tcts = cts;
	K = model->K_I[0];

	//printf("<<<<<<<<Begin reject sampling for b......\n");

	while((rep--) > 0){
		double newb[10];
		double b = cts->B[0], tot_t = 0;
		a = model->a[0];
		lgq = 0;

		for(i = 0; i < model->K_I[0]; i++){
			k_id = i;
			for(j = 0; j < cts->I; j++){
				double tot_cust = cts->top_cts[j].M[i];
				tot_t += cts->top_cts[j].T_K[i];
				if(tot_cust > 0){
					double q = 0;
					int try = 20;
					while(q <= 0 && try >= 0){
						q = gsl_ran_beta(glob_r, b, tot_cust);
						try--;
					}
					if(q <= 0){
						fprintf(stderr, "q in sample_b(b = %f, tot_cust = %f) went zero!!!\n", b, tot_cust);
						exit(1);
					}
					if(a > 0)
						lgq += gsl_sf_log(q);
					else
						lgq -= gsl_sf_log(q);
				}
			}
		}

		arms_simple (8, &lv, &rv, bterms_all, NULL, 0, &b, newb);

		for(i = 0; i < model->K_I[0]; i++){
			cts->B[i] = newb[0];
			model->b[i] = newb[0];
		}
		if(rep == 1)// && (i%10 == 0))
			printf("b = %lf, ", newb[0]);
	}
	printf("\n");
	tc = NULL;
	tcts = NULL;
}
/**********************************/
