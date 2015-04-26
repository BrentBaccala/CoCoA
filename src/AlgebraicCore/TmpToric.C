//   Copyright (c)  2011  Anna Bigatti

//   This file is part of the source of CoCoALib, the CoCoA Library.

//   CoCoALib is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.

//   CoCoALib is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.

//   You should have received a copy of the GNU General Public License
//   along with CoCoALib.  If not, see <http://www.gnu.org/licenses/>.


// Source code for toric ideals

#include "CoCoA/MatrixView.H"
#include "CoCoA/PPMonoid.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/TmpHilbert.H"
#include "CoCoA/TmpToric.H"
#include "CoCoA/BigInt.H"
#include "CoCoA/ideal.H"
#include "CoCoA/matrix.H"
#include "CoCoA/ring.H"
#include "CoCoA/utils.H"
#include "TmpHilbertDir/AnnaUtils.h"
#include "TmpHilbertDir/IVectors.h"
#include "TmpHilbertDir/eterms.h"
#include "TmpHilbertDir/unipoly.h" // for PoincareMaxPower
#include "TmpHilbertDir/poincare.h"
#include "TmpHilbertDir/toric.h"

#include <iostream> // for debugging only
#include <vector>
using std::vector;
#include <memory>
using std::auto_ptr;

namespace CoCoA
{

  //----------------------------------------------------------------------
  namespace // anonymous
  {
    ints IndicesForToric(const std::vector<long> v);
    ideal NewIdeal(const SparsePolyRing& P, const BitermList& BL);
    ideal NewIdeal(const SparsePolyRing& P, const BitermList& BL, const ints Weights);
    biterm PPs2Binom(ConstRefPPMonoidElem t1, ConstRefPPMonoidElem t2);
    BitermList NewBitermList(const ideal& I);
    void MatKerToBListAndIndices(ConstMatrixView TrM, BitermList *BL, ints *Indices, ints *Weights);
  }
  //----------------------------------------------------------------------


  ideal SequentialToric_C(const ideal& I, const std::vector<long>& indices)
  {
    // check input
    const SparsePolyRing P = RingOf(I);
    //    ::StartPoincare(0);
    ::StartToric(NumIndets(P));
    BitermList BLIn = NewBitermList(I);
    ints IndicesIn = IndicesForToric(indices);
    BitermList BLOut = SequentialToric(BLIn, IndicesIn);
    // BLIn and IndicesIn are freed by SequentialToric
    ideal T = NewIdeal(P, BLOut);
    EraseAndFreeBList(BLOut);
    return T;
  }


  ideal SequentialToric_C(const SparsePolyRing& P, ConstMatrixView M)
  {
    // check input
    //    ::StartPoincare(0);
    ::StartToric(NumIndets(P));
    BitermList BLIn;
    ints IndicesIn;
    ints WeightsIn;
    MatKerToBListAndIndices(transpose(M), &BLIn, &IndicesIn, &WeightsIn);    
    BitermList BLOut = SequentialToric(BLIn, IndicesIn);
    // BLIn and IndicesIn are freed by SequentialToric
    ideal T = NewIdeal(P, BLOut, WeightsIn);
    EraseAndFreeBList(BLOut);
    ints_free(WeightsIn);
    return T;
  }


  void EndToric_C()
  {
    EndPoincare_C(); // clear *C* global variable
  }
  


  namespace // anonymous
  {

    ints IndicesForToric(const std::vector<long> v)
    {
      long l = len(v);
      ints res=ints_malloc(l);
      IntsSetSize(res,l);
      IntsSetLen(res,0);
      for (long i=0; i<l; ++i)
        IntsPutLast(res, v[i]+1);
      return res;
    }


    RingElem BitermToRingElem(const SparsePolyRing& P, biterm b, ints weights)
    {
      ivec v = b->Vect;
      long n = ivec_len(v);
      if (n>NumIndets(P))
        CoCoA_ERROR("wrong number of indeterminates", "BitermToRingElem");
      PPMonoid M = PPM(P);
//   if ( ToricAlg == EATI )

      PPMonoidElem pp1(M);
      PPMonoidElem pp2(M);
      for ( long i=1 ; i<=n ; ++i )
        if (ivec_nth(v,i) > 0)
          M->myMulIndetPower(raw(pp1), i-1, ivec_nth(v,i)/weights[i]);
        else if (ivec_nth(v,i) < 0)
          M->myMulIndetPower(raw(pp2), i-1, -ivec_nth(v,i)/weights[i]);
      return monomial(P, 1, pp1) - monomial(P, 1, pp2);
    }


    ideal NewIdeal(const SparsePolyRing& P, const BitermList& BL)
    {
      biterms Bs = Biterms(BL);
      long n = BListLen(BL);
      std::vector<RingElem> res;

      ints weights1;
      weights1 = ints_malloc(NumIndets(P));
      IntsSetSize(weights1, NumIndets(P));
      IntsSetLen(weights1, NumIndets(P));
      for (long i=1 ; i<=NumIndets(P) ; i++ )  weights1[i] = 1;
      for (long i=1 ; i<=n ; i++ )
        //        if ( ToricAlg != EATI || ivec_nth(Bs[i]->Vect, IndNo)==0 )
        res.push_back(BitermToRingElem(P, Bs[i], weights1));
      ints_free(weights1);
      return ideal(P, res);
    }
    
    
    ideal NewIdeal(const SparsePolyRing& P, const BitermList& BL, const ints weights)
    {
      biterms Bs = Biterms(BL);
      long n = BListLen(BL);
      std::vector<RingElem> res;
      for (long i=1 ; i<=n ; i++ )
        //        if ( ToricAlg != EATI || ivec_nth(Bs[i]->Vect, IndNo)==0 )
        res.push_back(BitermToRingElem(P, Bs[i], weights));
      return ideal(P, res);
    }
    
    
    biterm PPs2Binom(ConstRefPPMonoidElem t1, ConstRefPPMonoidElem t2)
    {
      ivec Vect;
      ivec_elem n=NumIndets(owner(t1));
      int ord;
 
//   if ( ToricAlg == EATI )
//   {  
//     Vect = ivec_init(i+1);
//     ivec_set_nth(Vect, i+1, 0);
//     ord = i+1;
//   }
//   else
      {
        Vect = ivec_init(n);
        ord = 1;
      }
      vector<long> exps1, exps2;
      exponents(exps1, t1);
      exponents(exps2, t2);
      for ( long i=0 ; i<n ; ++i )
        ivec_set_nth(Vect, (unsigned long)i+1, exps1[i]-exps2[i]);
//       if ( term_degree(t1)!=term_degree(t2) )
//       {
//         ivec_free(Vect);
//         return NULL;
//       }
      return BitermMake(Vect, ord);
    }


    BitermList NewBitermList(const ideal& I)
    {
      const SparsePolyRing P = RingOf(I);
      const std::vector<RingElem>& GensI = gens(I);
      BitermList BL = BLNew(len(GensI), NumIndets(P));;
      biterms Bs = Biterms(BL);
      long MaxDeg = 0, BLLen = 0;

      RingElem g(P);
      for ( long i=0 ; i<len(GensI) ; ++i )
      {
        g = GensI[i];
        Bs[++BLLen] = PPs2Binom(LPP(g), LPP(g-monomial(P, LC(g), LPP(g))));
//       if ( Bs[BLLen]==NULL )
//       {
// 	EraseAndFreeBList(BL);
// 	math_error = cocoa_error("BAD_PARAMS");
// 	return NULL;
//       }
        BListSetLen(BL, BLLen);
        if ( BitermLtDeg(Bs[BLLen]) > MaxDeg)  MaxDeg = BitermLtDeg(Bs[BLLen]);
      }
      BListSetMaxDeg (BL, MaxDeg);
      
      return BL;
    }
    
    
    void MatKerToBListAndIndices(ConstMatrixView M, BitermList *BL, ints *Indices, ints *Weights)
    {
      int dim;
      int **small_basis, matrices =8, CurrLen =0, IndNo;
      ints Indices1, weights;  
      biterm B, MaxDegB;
      biterms Bs;
      ivec aux_V; 
      int Index, MaxDeg=0;

      long nrows = NumRows(M);
      long ncols = NumCols(M);
      
//   if ( ToricAlg == EATI )
//     IndNo = nrows+1;
//   else
    IndNo = nrows;

  /*  TODO CONTROLLA!!!!!!!!!!!!!!!!!!!!!!!!!!!
  rum_init(sizeof(expType)*(IndNo+1), RUM_STD_SIZE);
  rum_init(sizeof(expType)*(IndNo+2), RUM_STD_SIZE);
  */

  /* INDICES */
  Indices1 = ints_malloc(ncols);
  IntsSetSize(Indices1, ncols);
  IntsSetLen(Indices1, 0);
  *Indices = Indices1;
  INDICES = Indices1;
  
  /* M --> small_matrix */
  int **Msmall;
  Msmall = (int**)malloc(NumRows(M)*sizeof(int*));
  int tmp_int;
  BigInt tmp_BigInt;
  for (long i=0; i<NumRows(M); ++i)
  {
    Msmall[i] = (int*)malloc(ncols*sizeof(int));
    for (long j=0; j<NumCols(M); ++j)
    {
      if (!IsInteger(tmp_BigInt, M(i,j)))
        CoCoA_ERROR(ERR::BadArg, "Toric/MatKerToBListAndIndices");
      if (!IsConvertible(tmp_int, tmp_BigInt))
        CoCoA_ERROR(ERR::BadArg, "Toric/MatKerToBListAndIndices");
      Msmall[i][j] = tmp_int;
    }
  }
  /* WEIGHTS */
  weights = ints_init(nrows);
  IntsSetLen(weights, nrows);
  int i, j;
  for (i=0; i<nrows; ++i)
  {
    weights[i+1]=0;
    for (j=0; j<ncols; ++j)  weights[i+1] += Msmall[i][j];
  }
  i = 1;
  while ( weights[i+1]==weights[i] )  if ((++i)==nrows) break;
  if ((i)==nrows)
    for (j=1; j <= nrows; j++)  weights[j] =1;
  else
  {
    i = 0;
    while ( Msmall[i][0]!=0 )  if ((++i)==nrows) break;    
    if ((i)==nrows)
      for (j=0; j<nrows; j++)  weights[j+1] = Msmall[j][0];
  }
  *Weights = weights;  
  /* NULL_SPACE TO BLIST */
  dim = null_space( &small_basis, Msmall, nrows, ncols);  
  *BL = BLNew (matrices*dim, IndNo);
  BListSetLen (*BL, matrices*dim);
  Bs = Biterms(*BL);
  for ( i=0 ; i<dim ; i++)
  {
    aux_V =ivec_init(IndNo);
    for ( Index=0 ; Index<IndNo ; Index++ )
      if ( Index>=nrows )
        ivec_set_nth(aux_V, Index+1, 0);
      else
        ivec_set_nth(aux_V, Index+1, weights[Index+1]*small_basis[i][Index]);
    B = BitermMake(aux_V,1); // frees aux_V
    Bs[++CurrLen] = B;
    if ( BitermLtDeg(B) > MaxDeg) { MaxDeg = BitermLtDeg(B); MaxDegB = B; }
  }
  for (i=0; i<dim; i++) free(small_basis[i]);  free(small_basis);
  for ( matrices-- ; matrices>0 ; matrices-- )
  {    
    dim = null_space( &small_basis, Msmall, nrows, ncols);  
    for ( i=0 ; i<dim ; i++ )
    {
      aux_V =ivec_init(IndNo);
      for ( Index=0 ; Index<IndNo ; Index++ )
        if ( Index>=nrows )
          ivec_set_nth(aux_V, Index+1, 0);
        else
          ivec_set_nth(aux_V, Index+1, weights[Index+1]*small_basis[i][Index]);
      B = BitermMake(aux_V,1); // frees aux_V
      Bs[++CurrLen] = B; 
      if ( BitermLtDeg(B) > MaxDeg) { MaxDeg = BitermLtDeg(B); MaxDegB = B; }
    }
    for (i=0; i<dim; i++) free(small_basis[i]);  free(small_basis);
  }
  for (i=0; i<nrows ; i++) free(Msmall[i]);  free(Msmall);
  BListSetMaxDeg (*BL, MaxDeg);

  /*
  if (trouble) return -dim;
  return dim;
  */
}

  } // anonymous namespace



} // end of namespace CoCoA

// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/TmpToric.C,v 1.12 2014/07/30 14:11:54 abbott Exp $
// $Log: TmpToric.C,v $
// Revision 1.12  2014/07/30 14:11:54  abbott
// Summary: Changed name AmbientRing --> RingOf
// Author: JAA
//
// Revision 1.11  2014/07/07 13:16:30  abbott
// Summary: Removed AsSparsePolyRing
// Author: JAA
//
// Revision 1.10  2014/04/17 13:39:39  bigatti
// -- MatrixViews --> MatrixView
//
// Revision 1.9  2013/02/04 14:25:35  bigatti
// -- free weights1 (memory leak found with valgrind)
//
// Revision 1.8  2013/02/01 17:57:29  bigatti
// -- added EraseAndFreeBList in SequentialToric_C
//
// Revision 1.7  2013/01/31 17:11:37  bigatti
// -- renamed a variable
//
// Revision 1.6  2013/01/31 13:14:02  bigatti
// -- removed comments
//
// Revision 1.5  2013/01/31 13:13:32  bigatti
// -- added free
//
// Revision 1.4  2013/01/31 12:57:01  bigatti
// -- investigating memory leak - clarifying code
//
// Revision 1.3  2013/01/25 13:49:06  bigatti
// -- nrows->ncols in  Msmall[i] = (int*)malloc(ncols*sizeof(int));  (segv)
//
// Revision 1.2  2011/08/14 15:52:16  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.1  2011/05/09 14:48:11  bigatti
// -- first import
//
