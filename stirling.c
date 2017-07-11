/*
 * stirling.c
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>

#include "util.h"
#include "stable.h"

static stable_t** S_m;
static int *usedM, *usedN;

/*
 *
 * Note: 0 <= a < 1
 * 		 maxN >= maxM
 */
void make_stirling_table(int maxN, int maxM, double* a, int I) {
	int i;
	usedM = (int*)malloc(sizeof(int) * I);
	usedN = (int*)malloc(sizeof(int) * I);
	for(i = 0; i < I; ++i){
		usedM[i] = maxM;
		usedN[i] = maxN;
	}
	S_m = (stable_t**)malloc(sizeof(stable_t*) * I);
	for(i = 0; i < I; ++i){
		S_m[i] = S_make(maxN, maxM, a[i]);
	}
}

void make_stirling_table_I(int maxN, int maxM, double a, int i) {
	usedM[i] = maxM;
	usedN[i] = maxN;

	S_m[i] = S_make(maxN, maxM, a);
}

void update_stirling_table(double* a, int I)
{
	int i;
	for(i = 0; i < I; ++i){
		S_remake(S_m[i], usedM[i], a[i]);
	}
}

double stirling(int n, int m, double a, int i) {
	return S_safe(S_m[i], n, m);
}

void free_stirling_table(int I)
{
	int i;
	for(i = 0; i < I; ++i)
		S_free(S_m[i]);
	free(S_m);
	free(usedM);
	free(usedN);
	printf("Stirling table S_m is freed!!!\n");
}
