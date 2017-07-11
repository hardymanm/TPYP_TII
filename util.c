/*
 * util.c
 *
 */

#include <time.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_rng.h>

#include "util.h"

// global variable
gsl_rng* glob_r;

/*
 * wrapped vector operations
 *
 */
double vget(const gsl_vector* v, int i)
{
	return (gsl_vector_get(v, i));
}

void vset(gsl_vector* v, int i, double x)
{
	gsl_vector_set(v, i, x);
}

void vinc(gsl_vector* v, int i, double x)
{
	vset(v, i, vget(v, i) + x);
}

void vprint(gsl_vector * v)
{
	int i;
	for (i = 0; i < v->size; i++)
		printf("%0.9lf ", vget(v, i));
	printf("\n");
}

void scanf_vector(char* filename, gsl_vector* v)
{
	FILE* fileptr;
	fileptr = fopen(filename, "r");
	gsl_vector_fscanf(fileptr, v);
	fclose(fileptr);
}

void save_vector(char* filename, gsl_vector* v)
{
	FILE* fileptr;
	fileptr = fopen(filename, "w");
	gsl_vector_fprintf(fileptr, v, "%g");
	fclose(fileptr);
}

void normalize(gsl_vector* v)
{
	int size = v->size;
	double sum_v = vsum(v);
	int i;
	for (i = 0; i < size; i++)
		vset(v, i, vget(v, i) / sum_v);
}

double vsum(gsl_vector* v)
{
	double *data = v->data, val = 0;
	int size = v->size;
	int i;
	for (i = 0; i < size; i++)
		val += data[i];
	return (val);
}

int veq(const gsl_vector* v1, const gsl_vector* v2)
{
	int i;
	if(v1->size != v2->size){
		fprintf(stderr, "V_EQUAL: v1->size != v2->size !!!\n");
		exit(1);
	}
	double *data1 = v1->data;
	double *data2 = v2->data;
	for(i = 0; i < v1->size; i++){
		if(data1[i] != data2[i]){
			return 0;
		}
	}
	return 1;
}

int v_contain(const gsl_vector* v, double k)
{
	int i;
	for (i = 0; i < v->size; i++) {
		if (vget(v, i) == k) {
			return 1;
		}
	}
	return 0;
}

int v_non_zero(const gsl_vector* v)
{
	int i;
	int count = 0;
	for(i = 0; i < v->size; i++){
		if(vget(v, i) != 0){
			count++;
		}
	}
	return count;
}

/*
 * wrapped matrix operations
 *
 */
double mget(const gsl_matrix* m, int i, int j)
{
	return (gsl_matrix_get(m, i, j));
}

void mset(gsl_matrix* m, int i, int j, double x)
{
	gsl_matrix_set(m, i, j, x);
}

void minc(gsl_matrix* m, int i, int j, double x)
{
	mset(m, i, j, mget(m, i, j) + x);
}

void col_sum(gsl_matrix* m, gsl_vector* v)
{
	int i, j;

	if(m->size2 != v->size){
		fprintf(stderr, "COL_SUM: v->size != m->size2 !!!\n");
		exit(1);
	}

	gsl_vector_set_all(v, 0);
	for (i = 0; i < m->size1; i++)
		for (j = 0; j < m->size2; j++)
			vinc(v, j, mget(m, i, j));
}

void row_sum(gsl_matrix* m, gsl_vector* v)
{
	int i, j;
	if(m->size1 != v->size){
		fprintf(stderr, "ROW_SUM: v->size != m->size1 !!!\n");
		exit(1);
	}

	gsl_vector_set_zero(v);
	for (i = 0; i < m->size1; i++)
	{
		for( j = 0; j < m->size2; j++)
			vinc(v, i, mget(m, i, j));
	}
}

double msum(gsl_matrix* m)
{
	int i, j;
	double sum = 0;

	for(i = 0; i < m->size1; i++)
	{
		for(j = 0; j < m->size2; j++)
		{
			sum += mget(m, i, j);
		}
	}
	return sum;
}

void mprint(gsl_matrix * m)
{
	int i, j;
	for (i = 0; i < m->size1; i++)
	{
		for (j = 0; j < m->size2; j++)
			printf("%0.9lf ", mget(m, i, j));
		printf("\n");
	}
}

void scanf_matrix(char* filename, gsl_matrix * m)
{
	FILE* fileptr;
	fileptr = fopen(filename, "r");
	gsl_matrix_fscanf(fileptr, m);
	fclose(fileptr);
}

void save_matrix(char* filename, gsl_matrix* m)
{
	FILE* fileptr;
	fileptr = fopen(filename, "w");
	gsl_matrix_fprintf(fileptr, m, "%g");
	fclose(fileptr);
}

void printf_matrix(char* file, gsl_matrix * m)
{
	FILE* fileptr;
	int i, j;
	fileptr = fopen(file, "w");
	for(i = 0; i < m->size1; i++)
	{
		for(j = 0; j < m->size2; j++)
		{
			fprintf(fileptr, "%.9lf ", mget(m, i, j));
		}
		fprintf(fileptr, "\n\n");
	}
	fclose(fileptr);
}

/*
 * initial the random number generator
 *
 */
void initial_rng()
{
//	long seed = (long)1115574245;
	glob_r = gsl_rng_alloc(gsl_rng_taus);
	time_t seed;
	time(&seed);
	gsl_rng_set(glob_r, (long)seed);
	printf("The seed for glob_r: %ld\n", (long)seed);
}

/*
 * return a random double in the range 0 to 1, exclusive
 *
 */
double next_uniform()
{
	return (gsl_rng_uniform(glob_r));
}
/*
 * daw a single sample from (unnormalized) multinormial v, with normalizing factor m
 *
 */
int next_discrete_unnormalised(gsl_vector* v, double sum)
{
	int i;
	double b = 0;
	double r = next_uniform() * sum;
//	printf(">>>>>>>>>>>>>>>r = %lf\n", r);
	for (i = 0; i < v->size; i++)
	{
		b += vget(v, i);
		if (b > r) {
			return (i);
		}
	}
	return (v->size - 1);
}
/*
 * daw a single sample from multinormial v
 *
 */
int next_discrete_normalised(double* v, int k)
{
	int i;
	double b = 0;
	double r = next_uniform();
	for (i = 0; i < k; i++){
		b += v[i];
		if (b > r) {
			return (i);
		}
	}
	return (k - 1);
}

double max_val(double* vec, int k)
{
	int i;
	double max = vec[0];
	for(i = 1; i < k; ++i){
		if(vec[i] > max){
			max = vec[i];
		}
	}
	return max;
}

void free_rng()
{
	gsl_rng_free(glob_r);
}

double safe_log(double x)
{
	if (x == 0)
		return (-1000);
	else
		return (log(x));
}

double log_sum(double V, double lp)
{
	if (lp > V) {
		// swap so V is bigger
		double t = lp;
		lp = V;
		V = t;
	}
	return V + log(1.0 + exp(lp - V));
}

/*
 * (x|y)_n
 */
double poch_sym(double x, double y, int n)
{
	int i;
	double tmp = 1.0;
	for(i = 0; i < n; i++){
		tmp *= (x + i*y);
	}
	return tmp;
}
/*
 * (x|y)_n
 */
double log_poch_sym(double x, double y, int n)
{
	double log_poch = 0;

	if(y == 0){
		log_poch = n * gsl_sf_log(x);
	}else{
		double tmp = x/y;
		log_poch = gsl_sf_lngamma(tmp + n) - gsl_sf_lngamma(tmp) + n * gsl_sf_log(y);
	}
	return log_poch;
}

/*
 * (x|y)_n
 */
double log_poch_sym_1(double x, double y, int n)
{
	int i;
	double log_poch = 0;
	for(i = 0; i < n; i++)
	{
		log_poch += gsl_sf_log(x + i*y);
	}
	return log_poch;
}
/* int main(){

 int i, k;
 double vsum;
 gsl_vector* test_vector = gsl_vector_calloc(10);

 initial_rng();
 for(i = 0; i < test_vector->size; i ++){
 vset(test_vector, i, next_uniform());
 }
 vprint(test_vector);
 vsum = sum(test_vector);

 printf("sum = %lf\n", vsum);

 k = next_discrete_unnormalised(test_vector, vsum);
 printf("k: %d\n", k);
 k = next_discrete_unnormalised(test_vector, vsum);
 printf("k: %d\n", k);
 k = next_discrete_unnormalised(test_vector, vsum);
 printf("k: %d\n", k);
 k = next_discrete_unnormalised(test_vector, vsum);
 printf("k: %d\n", k);
 k = next_discrete_unnormalised(test_vector, vsum);
 printf("k: %d\n", k);
 k = next_discrete_unnormalised(test_vector, vsum);
 printf("k: %d\n", k);
 k = next_discrete_unnormalised(test_vector, vsum);
 printf("k: %d\n", k);
 k = next_discrete_unnormalised(test_vector, vsum);
 printf("k: %d\n", k);
 return 0;
 }*/
