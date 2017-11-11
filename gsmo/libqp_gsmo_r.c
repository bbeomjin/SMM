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

SEXP gsmo(SEXP mat_H1, SEXP f1, SEXP a1, SEXP b1,  SEXP LB1, SEXP UB1,  SEXP vec_x01, SEXP MaxIter1, SEXP TolKKT1, SEXP verb1)
{ 
  int verb;          
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
  double fval;    
  SEXP Rdim, temp;
  
  mat_H = REAL(mat_H1);
  Rdim = getAttrib(mat_H1, R_DimSymbol);
  nVar = INTEGER(Rdim)[1];
  f = REAL(f1);   
  a = REAL(a1);
  b = REAL(b1)[0];
  LB = REAL(LB1);
  UB = REAL(UB1);
  vec_x0 = REAL(vec_x01);
  MaxIter = (long)REAL(MaxIter1)[0];
  TolKKT = REAL(TolKKT1)[0]; 
  verb = (int)REAL(verb1)[0]; 
  
  PROTECT(temp = allocVector(REALSXP, nVar));
  vec_x = REAL(temp);
  memcpy( vec_x, vec_x0, sizeof(double)*nVar );
  
  diag_H = calloc(nVar, sizeof(double));
  for(i = 0; i < nVar; i++ ) 
    diag_H[i] = mat_H[nVar*i+i];

  state = libqp_gsmo_solver( &get_col, diag_H, f, a, b, LB, UB, vec_x, 
                        nVar, MaxIter, TolKKT, NULL );
						
  UNPROTECT(1);
  return(temp);
}



