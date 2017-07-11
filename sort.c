/*
 * sort.c
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>

#include "util.h"
/*
 *
 *
 */
int choose_pivot(int left, int right) {
	return ((left + right) / 2);
}

/*
 *
 *
 */
void swap(pair* x, pair* y) {
	pair tmp;

	tmp.key = x->key;
	tmp.value = x->value;
	x->key = y->key;
	x->value = y->value;
	y->key = tmp.key;
	y->value = tmp.value;
}

/*
 *
 *
 */
void quick_sort(pair* pairs, int m, int n) {

	int pivot_index;
	double val;
	int i, j;

	val = 0;
	if (m < n) {
		pivot_index = choose_pivot(m, n);
		swap(&pairs[m], &pairs[pivot_index]);
		val = pairs[m].value;
		i = m + 1;
		j = n;
		while (i <= j) {
			//			while( (i <= n) && (pairs[i].value <= val)){
			while ((i <= n) && (pairs[i].value >= val)) {
				i++;
			}
			//			while((j >= m) && (pairs[j].value > val)){
			while ((j >= m) && (pairs[j].value < val)) {
				j--;
			}
			if (i < j)
				swap(&pairs[i], &pairs[j]);
		}
		swap(&pairs[m], &pairs[j]);
		quick_sort(pairs, m, j - 1);
		quick_sort(pairs, j + 1, n);
	}
}

/*
 *
 *
 */
void sort(gsl_vector* vector) {
	int i;
	pair* pairs = (pair*) malloc(sizeof(pair) * vector->size);
	for (i = 0; i < vector->size; i++) {
		pairs[i].key = i;
		pairs[i].value = gsl_vector_get(vector, i);
	}
////		printf("before quick sort ...\n");
////		for (i = 0; i < vector->size; i++) {
////			printf("%d : %d\n", pairs[i].key, pairs[i].value);
////		}
//
	quick_sort(pairs, 0, vector->size - 1);

////		 	printf("after quick sort ...\n");
	for (i = 0; i < vector->size; i++) {
////				printf("%d : %d\n", pairs[i].key, pairs[i].value);
		vset(vector, i, pairs[i].key);
	}
	free(pairs);
}

void sort1(gsl_vector* vector, gsl_vector* value)
{
	int i;
	pair* pairs = (pair*) malloc(sizeof(pair) * vector->size);
	for (i = 0; i < vector->size; i++) {
		pairs[i].key = i;
		pairs[i].value = gsl_vector_get(vector, i);
	}

	quick_sort(pairs, 0, vector->size - 1);

	for (i = 0; i < vector->size; i++) {
		vset(vector, i, pairs[i].key);
		vset(value, i, pairs[i].value);
	}
	free(pairs);
}

/*int main() {
	int i;
	gsl_vector* vector = gsl_vector_alloc(10);
	for (i = 0; i < vector->size; i++) {
		gsl_vector_set(vector, i, rand()%100);
	}

	printf("before sort ...\n");
	for (i = 0; i < vector->size; i++) {
		printf("%d ", (int) gsl_vector_get(vector, i));
	}
	printf("\n");
	sort(vector);
	printf("after sort ...\n");
	for (i = 0; i < vector->size; i++) {
			printf("%d ", (int) gsl_vector_get(vector, i));
		}
	printf("\n");
	gsl_vector_free(vector);
	return 0;
}*/
