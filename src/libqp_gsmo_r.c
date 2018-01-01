#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "libqp.h"



#define INDEX(ROW,COL,NUM_ROWS) ((COL*NUM_ROWS)+ROW)
#define MIN(A,B) ((A < B) ? A : B)
#define MAX(A,B) ((A > B) ? A : B)


/* ------------------------------------------------------------*/
/* Declaration of global variables                             */
/* ------------------------------------------------------------*/
double *mat_H;
uint32_t nVar;
libqp_state_T state;

/* ------------------------------------------------------------
  Returns pointer at the a-th column of the matrix H.
------------------------------------------------------------ */
const double *get_col( uint32_t i )
{
  return( &mat_H[ nVar*i ] );
}


libqp_state_T libqp_gsmo_solver(const double* (*get_col)(uint32_t),
            double *diag_H,
            double *f,
            double *a,
            double b,
            double *LB,
            double *UB,
            double *x,
            uint32_t n,
            uint32_t MaxIter,
            double TolKKT,
            void (*print_state)(libqp_state_T state));

SEXP gsmo(SEXP mat_H_r, SEXP f_r, SEXP a_r, SEXP b_r,  SEXP LB_r, SEXP UB_r,  SEXP vec_x0_r, SEXP MaxIter_r, SEXP TolKKT_r, SEXP verb_r)
{           
  uint32_t i;         
  uint32_t MaxIter;   
  double TolKKT; 
  double *vec_x;         /* output arg -- solution*/ 
  double *vec_x0;         
  double *diag_H;    /* diagonal of matrix H */
  double *f;         /* vector f */
  double *a;
  double b;
  double *LB;
  double *UB;
  /*uint32_t verb;*/
  SEXP Rdim, temp;
  
  mat_H = REAL(mat_H_r);
  Rdim = getAttrib(mat_H_r, R_DimSymbol);
  nVar = INTEGER(Rdim)[1];
  f = REAL(f_r);   
  a = REAL(a_r);
  b = Rf_asReal(b_r);
  LB = REAL(LB_r);
  UB = REAL(UB_r);
  vec_x0 = REAL(vec_x0_r);
  MaxIter = (long)Rf_asInteger(MaxIter_r);
  TolKKT = Rf_asReal(TolKKT_r); 
  /*verb = Rf_asInteger(verb_r);*/
  
  const char *names[] = {"alpha", "obj_value", "exitflag", "nIter", ""};
  
  SEXP res = PROTECT(mkNamed(VECSXP, names));
  PROTECT(temp = allocVector(REALSXP, nVar));
  vec_x = REAL(temp);
  memcpy( vec_x, vec_x0, sizeof(double)*nVar );
  
  diag_H = calloc(nVar, sizeof(double));
  for(i = 0; i < nVar; i++ ) { 
    diag_H[i] = mat_H[nVar*i+i];
  }
  
  state = libqp_gsmo_solver( &get_col, diag_H, f, a, b, LB, UB, vec_x, 
                        nVar, MaxIter, TolKKT, NULL );
						
  SET_VECTOR_ELT(res, 0, temp);
  SET_VECTOR_ELT(res, 1, ScalarReal(state.QP));
  SET_VECTOR_ELT(res, 2, ScalarReal(state.exitflag));
  SET_VECTOR_ELT(res, 3, ScalarInteger(state.nIter)); 
  UNPROTECT(2);
  return(res);
}



