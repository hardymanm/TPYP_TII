#include <gsl/gsl_sf.h>
#include <math.h>

#include "sym_polya_fit.h"
#include "util.h"

#define maxIt 100
#define threshold 1e-4
#define GAMMA_MAX 10
#define GAMMA_MIN 0.001

double sym_polya_fit(double**counts, double* count_sum, int size1, int size2, double initV)
{
	int k, v, iter;
	double sum1, sum2, gamma_old, gamma;
	iter = 0;
	gamma = gamma_old = initV;
	while(iter++ < maxIt){
		sum1 = 0;
		sum2 = 0;
		gamma_old = gamma;
		for(k = 0; k < size1; k++){
			for(v = 0; v < size2; v++){
				sum1 += gsl_sf_psi(counts[k][v] + gamma_old);
			}
			sum2 += gsl_sf_psi(count_sum[k] + size2 * gamma_old);
		}
		sum1 -= size1 * size2 * gsl_sf_psi(gamma_old);
		sum2 *= size2;
		sum2 -= size1 * size2 * gsl_sf_psi(size2 * gamma_old);
		gamma = gamma_old * sum1 / sum2;
		if(gamma > GAMMA_MAX){
			gamma = 0.5;
			break;
		}
		if(gamma < GAMMA_MIN){
			gamma = GAMMA_MIN;
			break;
		}
		if(iter > 20 && fabs(gamma - gamma_old) / gamma_old < threshold){
			break;
		}
	}
	return gamma;
}

double sym_polya_fit_newton(double**counts, double* count_sum, int size1, int size2, double initV)
{
	int k, v, iter;
	double sum1, sum2, gamma_old, gamma;
	iter = 0;
	gamma = gamma_old = initV;
	while(iter++ < maxIt){
		sum1 = 0;
		sum2 = 0;
		gamma_old = gamma;
		for(k = 0; k < size1; k++){
			for(v = 0; v < size2; v++){
				sum1 += gsl_sf_psi(counts[k][v] + gamma_old);
				sum2 += gsl_sf_psi_1(counts[k][v] + gamma_old);
			}
			sum1 -= size2 * gsl_sf_psi(count_sum[k] + size2 * gamma_old);
			sum2 -= (double)size2 * (double)size2 * gsl_sf_psi_1(count_sum[k] + size2 * gamma_old);
		}
		sum1 += size1 * size2 * (gsl_sf_psi(size2 * gamma_old) - gsl_sf_psi(gamma_old));
		sum2 += size1 * size2 * (size2 * gsl_sf_psi_1(size2 * gamma_old) - gsl_sf_psi_1(gamma_old));
		gamma -= sum1 / sum2;
		if(gamma > GAMMA_MAX){
			gamma = 0.5;
			break;
		}
		if(gamma < GAMMA_MIN){
			gamma = GAMMA_MIN;
			break;
		}
		if(iter > 20 && fabs(gamma - gamma_old) / gamma_old < threshold){
			break;
		}
	}
	return gamma;
}

double sym_polya_fit1(int**counts, int* count_sum, int size1, int size2, double initV)
{
	int k, v, iter;
	double sum1, sum2, gamma_old, gamma;
	iter = 0;
	gamma = gamma_old = initV;
	while(iter++ < maxIt){
		sum1 = 0;
		sum2 = 0;
		gamma_old = gamma;
		for(k = 0; k < size1; k++){
			for(v = 0; v < size2; v++){
				sum1 += gsl_sf_psi(counts[k][v] + gamma_old);
			}
			sum2 += gsl_sf_psi(count_sum[k] + size2 * gamma_old);
		}
		sum1 -= size1 * size2 * gsl_sf_psi(gamma_old);
		sum2 *= size2;
		sum2 -= size1 * size2 * gsl_sf_psi(size2 * gamma_old);
		gamma = gamma_old * sum1 / sum2;
		if(gamma > GAMMA_MAX){
			gamma = 0.5;
			break;
		}
		if(gamma < GAMMA_MIN){
			gamma = GAMMA_MIN;
			break;
		}
		if(iter > 20 && fabs(gamma - gamma_old) / gamma_old < threshold){
			break;
		}
	}
	return gamma;
}

double sym_polya_fit1_newton(int**counts, int* count_sum, int size1, int size2, double initV)
{
	int k, v, iter;
	double sum1, sum2, gamma_old, gamma;
	iter = 0;
	gamma = gamma_old = initV;
	while(iter++ < maxIt){
		sum1 = 0;
		sum2 = 0;
		gamma_old = gamma;
		for(k = 0; k < size1; k++){
			for(v = 0; v < size2; v++){
				sum1 += gsl_sf_psi(counts[k][v] + gamma_old);
				sum2 += gsl_sf_psi_1(counts[k][v] + gamma_old);
			}
			sum1 -= size2 * gsl_sf_psi(count_sum[k] + size2 * gamma_old);
			sum2 -= (double)size2 * (double)size2 * gsl_sf_psi_1(count_sum[k] + size2 * gamma_old);
		}
		sum1 += size1 * size2 * (gsl_sf_psi(size2 * gamma_old) - gsl_sf_psi(gamma_old));
		sum2 += size1 * size2 * (size2 * gsl_sf_psi_1(size2 * gamma_old) - gsl_sf_psi_1(gamma_old));
		gamma -= sum1 / sum2;
		if(gamma > GAMMA_MAX){
			gamma = 0.5;
			break;
		}
		if(gamma < GAMMA_MIN){
			gamma = GAMMA_MIN;
			break;
		}
		if(iter > 20 && fabs(gamma - gamma_old) / gamma_old < threshold){
			break;
		}
	}
	return gamma;
}

void gamma_opt(Cts* cts, Model* model, Corpus* c)
{
	double gamma = sym_polya_fit(cts->T_KQ, cts->sum_T, model->K_I[0], model->v, vget(model->gamma, 0));
	gsl_vector_set_all(model->gamma, gamma);
	model->gamma_sum = gamma * model->v;
}

void gamma_opt_newton(Cts* cts, Model* model, Corpus* c)
{
	double gamma = sym_polya_fit_newton(cts->T_KQ, cts->sum_T, model->K_I[0], model->v, vget(model->gamma, 0));
	gsl_vector_set_all(model->gamma, gamma);
	model->gamma_sum = gamma * model->v;
}

void alpha_opt(Cts* cts, Model* model, int* nDoc, int K)
{
	int i;
	double alpha;
	for(i = 0; i < model->I; i++){
		alpha = sym_polya_fit1(cts->n[i], cts->N[i], nDoc[i] , K, vget(model->alpha[i], 0));
		gsl_vector_set_all(model->alpha[i], alpha);
		model->alpha_sum[i] = alpha * K;
	}
}

void alpha_opt_newton(Cts* cts, Model* model, int* nDoc, int K)
{
	int i;
	double alpha;
	for(i = 0; i < model->I; i++){
		alpha = sym_polya_fit1_newton(cts->n[i], cts->N[i], nDoc[i] , K, vget(model->alpha[i], 0));
		gsl_vector_set_all(model->alpha[i], alpha);
		model->alpha_sum[i] = alpha * K;
	}
}
