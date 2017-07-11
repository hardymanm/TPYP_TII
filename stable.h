/*
 * Stirling Number table handling for Pitman-Yor processing
 * Copyright (C) 2009-2012 Wray Buntine 
 * All rights reserved.
 *
 * The contents of this file are subject to the Mozilla Public License
 * Version 1.1 (the "License"); you may not use this file except in
 * compliance with the License. You may obtain a copy of the License at
 * http://www.mozilla.org/MPL/
 *
 * Software distributed under the License is distributed on an "AS IS"
 * basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
 * License for the specific language governing rights and limitations
 * under the License.
 *
 * Author: Wray Buntine (wray.buntine@nicta.com.au)
 *
 *   This is a simplified version of earlier code
 *   released to accompany this simple software tester.
 *   For definitions and details, see
 *   	http://arxiv.org/abs/1007.0296
 *
 */

#ifndef __STABLE_H
#define __STABLE_H

/*
 *  data stored here
 */
typedef struct stable_s {
  int usedM, usedN;        /*  current dimensions */
  double **S;
  double a;
} stable_t;

/*
 *  access a Stirling Number safely, but gives the log() of it,
 *   i.e.,  does bounds checking and extends table if needed
 */
double S_safe(stable_t *sp, int N, int M);
/*
 *  fills S table with maxN and maxM values;
 *  return NULL on error
 */
stable_t *S_make(int maxN, int maxM, double a);
/*
 *  fill with new values for different "a"
 *  extends S table using larger maxM value,
 *  or leave same size (if maxM==0);
 *  return non-zero on error
 */
int S_remake(stable_t *sp, int maxM, double a);

void S_free(stable_t *sp);

#endif
