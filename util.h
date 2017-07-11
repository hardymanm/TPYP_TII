/*
 * util.h
 *
 */

#ifndef UTIL_H_
#define UTIL_H_

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <math.h>

#define MAX_ITE_RAN 5


typedef struct pair{
    int key;
    double value;
} pair;

/*
 * vector operations
 */
double vget(const gsl_vector*, int);
void vset(gsl_vector*, int, double);
void vinc(gsl_vector*, int, double);
void vprint(gsl_vector*);
void scanf_vector(char*, gsl_vector*);
void save_vector(char*, gsl_vector*);
double vsum(gsl_vector*);
void normalize(gsl_vector* v);

int veq(const gsl_vector*, const gsl_vector*);
int v_contain(const gsl_vector* v, double);
int v_non_zero(const gsl_vector* v);
/*
 * matrix operations
 */
double mget(const gsl_matrix* m, int, int);
void mset(gsl_matrix*, int, int, double);
void minc(gsl_matrix*, int, int, double);
void col_sum(gsl_matrix*, gsl_vector*);
void row_sum(gsl_matrix*, gsl_vector*);
double msum( gsl_matrix*);

void mprint(gsl_matrix*);
void scanf_matrix(char*, gsl_matrix*);
void save_matrix(char*, gsl_matrix* m);
void printf_matrix(char*, gsl_matrix*);

/*
 * sampling operations
 */
void initial_rng();
void free_rng();
double next_uniform();
int next_discrete_unnormalised(gsl_vector*, double);
int next_discrete_normalised(double*, int);
double max_val(double*, int);

void sort(gsl_vector*);
void sort1(gsl_vector*, gsl_vector*);
double safe_log(double);
double log_sum(double, double);
/*
 * stirling number generation
 */
void make_stirling_table(int maxN, int maxM, double* a, int I);
void update_stirling_table(double* a, int I);
void make_stirling_table_I(int maxN, int maxM, double a, int i);
double stirling(int n, int m, double a, int i);
void free_stirling_table(int I);
/*
 * Pochhammer Symbal
 */
double poch_sym(double, double, int);
double log_poch_sym(double, double, int);

#endif /* UTIL_H_ */
