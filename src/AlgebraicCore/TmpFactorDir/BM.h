#ifndef BMdeg_h
#define BMdeg_h

#include <gmp.h>
#include "jaa.h"

/***************************************************************************/

struct BuchbergerMoeller_struct
{
  int npoints;   /* number of points                                       */
  int nvars;     /* number of variables = dimension of the space           */
  int npp;       /* number of power products currently in use              */
  int npp_max;   /* maximum number of power products we can manage         */
                 /* this value can be increased by calling BM_grow         */
  int maxdeg;    /* highest degree PP used in projective case; 0 for affine*/
  int *det;      /* determinant of the separators accumulates here         */
  int GBsize;    /* number of elements in the Groebner basis               */
  int *GB;       /* indices of the rows of M which corr to GB elems        */
  int *sep;      /* indices of the rows of M which corr to separators      */
  int *projective_sep; /* indices of the rows of M which corr to separators*/
  int *issep;    /* flag, true iff the corr PP is a LPP of a separator     */
  int **pp;      /* the power product which gave rise to this row of M     */
  FFelem **M;    /* matrix of the values of the PPs at the given points,   */
                 /* subsequently row-reduced.                              */
  int (*cmp)(const void*,const void*); /* power product order              */
};
typedef struct BuchbergerMoeller_struct *BM;


/***************************************************************************/

struct BMGB_struct
{
  int npoints;
  int GBsize;
  int compute_GB, compute_separators; /* booleans */
  int last_fail; /* used only in BM_check and BM_check_with_det */
  int target_height; /* height to reach before checking, log base 2 */
  mpz_t modulus;
  mpz_t *det;
  mpz_t **GB_mod;
  mpq_t **GB;
  mpz_t **sep_mod;
  mpq_t **sep;
};
typedef struct BMGB_struct *BMGB;


/***************************************************************************/
/* This is the public interface.                                           */
/* THESE FUNCTIONS REQUIRE AT LEAST 1 INPUT POINT: the case of 0 input     */
/* points is not handled gracefully.                                       */

void BM_affine(BMGB *char0_ptr, BM *modp_ptr, int nvars, int npoints, mpz_t **points, int (*cmp)(const void*, const void*));
void BM_projective(BMGB *char0_ptr, BM *modp_ptr, int nvars, int npoints, mpz_t **points, int (*cmp)(const void*, const void*));
void BM_dtor(BM self);
void BMGB_dtor(BMGB self);
BM BM_affine_mod_p(int nvars, int npoints, FFelem **points, int (*cmp)(const void*, const void*));
BM BM_projective_mod_p(int nvars, int npoints, FFelem **points, int (*cmp)(const void*, const void*));

#if 0
void BM_print_GBasis(BMGB char0, BM modp);
void BM_print_separators(BMGB char0, BM modp);
#endif

#endif
