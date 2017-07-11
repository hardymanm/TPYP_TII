/*
 * Stirling Number table handling
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
 *
 *   The full package allows one to compute derivatives as well as
 *   implementing various tricks to support sparse tables.
 *     
 */
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "stable.h"

static double logadd(double V, double lp) {
  if ( lp>V ) {
    // swap so V is bigger
    double t = lp;
    lp = V;
    V = t;
  }
  return V + log(1.0+exp(lp-V));
}

/*
 *   stores the S table as a float
 */
#define tblSNM(N,M) sp->S[M][N]

/*
 *  fills S table with arrays of size maxN, which is assumed
 *  to be a constant for the full program;
 *  fills M to maxM and stores M and a for later
 *
 *  return NULL on error (memory allocation trouble)
 */
stable_t *S_make(int maxN, int maxM, double a) {
  int N, M;
  stable_t *sp = malloc(sizeof(stable_t));
  if ( !sp ) 
    return NULL;

  if ( maxN<maxM ) {
    /*  Should have maxN>maxM in S_make */
    maxM = maxN;
  }

  sp->usedN = maxN;
  sp->usedM = maxM;
  
  sp->S = malloc(sizeof(sp->S[0])*(maxM+1));
  if ( !sp->S ) {
    free(sp);
    return NULL;
  }
  sp->S[0] = NULL;
  for (M=1; M<=maxM; M++) {
    sp->S[M] = malloc(sizeof(sp->S[0][0])*(1+maxN));
    if ( !sp->S[M] ) {
      int M2;
      for (M2=1; M2<M; M2++)
	free(sp->S[M2]);
      free(sp->S);
      free(sp);
      return NULL;
    }
  }

  /*
   *  all values outside bounds to log(0);
   *  when N=M set to log(1)
   */
  for (N=1; N<maxN; N++) {
    if ( N<=maxM )
        tblSNM(N,N) = 0;
    for (M=N+1; M<=maxM; M++)
      tblSNM(N,M) =  -HUGE_VAL;
  }

  /*
   *  this is where we actually build the Stirling numbers
   */
  sp->S[1][1] = 0;
  S_remake(sp,0,a);
  return sp;
}

/*
 *    assume tables filled to usedM;
 *    fiddle maxM to make it non-trivial
 *    create extra space first;
 *    then extend
 *    return non-zero on error
 */
int S_extend(stable_t *sp, int maxM) {
  int N, M;
  /*
   *  shouldn't be too big 
   */
  if ( sp->usedN<maxM ) {
    /*  Should have maxN>maxM in S_extend */
    maxM = sp->usedN;
  }  
  /*
   *   increase should not be trivial
   */
  if ( maxM<sp->usedM*1.1 )
    maxM = sp->usedM*1.1;
  if ( maxM<sp->usedM+10 )
    maxM = sp->usedM+10;
  /*
   *  reset if made too big
   */
  if ( sp->usedN<maxM ) {
    maxM = sp->usedN;
  }

  sp->S = realloc(sp->S,sizeof(sp->S[0])*(maxM+1));
  if ( !sp->S ) {
    free(sp);
    return 1;
  }
  for (M=sp->usedM+1; M<=maxM; M++) {
    sp->S[M] = malloc(sizeof(sp->S[0][0])*(sp->usedN+1));
    if ( !sp->S[M] ) {
      int M2;
      for (M2=1; M2<M; M2++)
	free(sp->S[M2]);
      free(sp->S);
      free(sp);
      return 1;
    }
  }
  
  /*
   *  all values outside bounds to log(0);
   *  when N=M set to log(1)
   */
  for (N=sp->usedM+1; N<=maxM; N++) {
    tblSNM(N,N) = 0;
  }
  for (M=sp->usedM+1; M<=maxM; M++) {
    for (N=1; N<M; N++) 
      tblSNM(N,M) =  -HUGE_VAL;
  }
  
  for (N=sp->usedM+1; N<sp->usedN; N++) {
    for (M=sp->usedM+1; M<=maxM && M<N; M++) {
      tblSNM(N,M) = logadd(tblSNM(N-1,M-1),log(N-(M*sp->a)-1.0)+tblSNM(N-1,M));
      assert(isfinite(sp->S[M][N]));
    }
  }    

  sp->usedM = maxM;
  return 0;
}

/*
 *   just rebuild to maxM
 *   if maxM==0, rebuild to current usedM
 *   if maxM>usedM, then extend as well;
 *   then completely refill the tables with the new a;
 *   return non-zero on error
 */
int S_remake(stable_t *sp, int maxM, double a) {
  int N, M;
  int save_maxM = 0;
  if ( maxM<=0 )
    maxM = sp->usedM;
  if ( sp->usedN<maxM ) {
    maxM = sp->usedN;
  }
  if ( maxM>sp->usedM ) {
    /*
     *  bigger than it was, so build to the
     *  old usedM first, then extend
     */
    save_maxM = maxM;
    maxM = sp->usedM;
  }
  assert(maxM<=sp->usedN);

  sp->a = a;
  for (N=2; N<sp->usedN; N++) {
    sp->S[1][N] = log(N-a-1.0) + sp->S[1][N-1];
    for (M=2; M<=maxM && M<N; M++) {
      tblSNM(N,M) = logadd(tblSNM(N-1,M-1),log(N-(M*a)-1.0)+tblSNM(N-1,M));
      assert(isfinite(sp->S[M][N]) );
    }
  }   
  /*
   *   if usedM was > maxM, shrink
   *   unused stuff
   */
  for (M=maxM+1; M<=sp->usedM; M++) {
    free(sp->S[M]);
    sp->S[M] = NULL;
  }
  sp->usedM = maxM;

  if ( save_maxM>0 )
    if ( S_extend(sp, save_maxM) )
      return 1;
  return 0;
}
  
double S_safe(stable_t *sp, int N, int T) {
  if ( N>=sp->usedN  ) {
    return -HUGE_VAL;
  }
  if ( T>sp->usedM ) {
    S_extend(sp,T);
  }
  if ( N==T )
    return 0;
  if ( N<T || T==0 )
    return -HUGE_VAL;
  return tblSNM(N,T);
}

void S_free(stable_t *sp) {
  if ( !sp )
    return;
  if ( sp->S ) {
    int m;
    for (m=1; m<=sp->usedM; m++)
      free(sp->S[m]);
    free(sp->S);
  }
  free(sp);
}

