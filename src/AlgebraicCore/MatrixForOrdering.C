//   Copyright (c)  2008-2012  John Abbott, Anna M. Bigatti

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


#include "CoCoA/MatrixForOrdering.H"
#include "CoCoA/apply.H"
#include "CoCoA/assert.H"
#include "CoCoA/error.H"
#include "CoCoA/matrix.H"
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/IntOperations.H"
#include "CoCoA/MatrixOperations.H"
#include "CoCoA/MatrixView.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/utils.H"

#include <iostream>
using std::ostream;
#include <vector>
using std::vector;

namespace CoCoA
{

  bool IsTermOrdering(ConstMatrixView M)
  {
    if  (!IsZZ(RingOf(M)) && !IsQQ(RingOf(M)))
      CoCoA_ERROR("matrix must be over RingZZ() or RingQQ()", "IsTermOrdering");
    for (long col=0; col<NumCols(M); ++col)
    {
      long row=0;
      for (; row<NumRows(M); ++row)
      {
        if ( IsZero(M(row,col)) ) continue;
        if ( sign(M(row,col))<0 ) return false;
        else break;
      }
      if (row==NumRows(M)) return false;
    }
    return true;
  }


/**
    expects a matrix with entries in an ordered ring: very likely only Z (maybe Q?)
    returns a matrix with positive entries which defines an equivalent ordering
*/

  matrix NewPositiveMat(ConstMatrixView M)
  {
    if  (!IsZZ(RingOf(M)) && !IsQQ(RingOf(M)))
      CoCoA_ERROR("matrix must be over RingZZ() or RingQQ()", "NewPositiveMatrix");
    if ( !IsTermOrdering(M) )
      CoCoA_ERROR(ERR::NotTermOrdering, "PositiveMatrix");

    matrix PM(NewDenseMat(M));
    const RingElem& unity(one(RingOf(PM)));
    for (long row=0; row < NumRows(PM); ++row)
      for (long col=0; col < NumCols(PM); ++col)
        if (sign(PM(row,col)) < 0)
        {
          long PosRow=0;
          while (IsZero(PM(PosRow,col))) ++PosRow; // find positive entry in col
          while (sign(PM(row,col)) < 0) PM->myAddRowMul(row, PosRow, unity);
        }
    return PM;
  }


  RingElem LCMofDenInRow(ConstMatrixView M, long row)
  {
    if (!IsQQ(RingOf(M)))
      CoCoA_ERROR("matrix must be over RingQQ()", "LCMofDenInRow");
    RingElem D = one(RingZZ());
    for (long col=0; col < NumCols(M); ++col)
      D = lcm(D, den(M(row,col)));
    return EmbeddingHom(RingQQ())(D);
  }
  

  matrix NewIntegerOrdMat(ConstMatrixView M)
  {
    if  (IsZZ(RingOf(M)))  return NewDenseMat(M);
    if (!IsQQ(RingOf(M)))
      CoCoA_ERROR("matrix must be over RingZZ() or RingQQ()", "NewIntegerOrdMat");
    if (!IsTermOrdering(M)) CoCoA_ERROR(ERR::NotTermOrdering, "NewIntegerOrdMat");

    matrix IM(NewDenseMat(M));
    for (long row=0; row < NumRows(IM); ++row)
      for (long col=0; col < NumCols(IM); ++col)
        if (!IsOne(den(IM(row,col))))
        {
          IM->myRowMul(row, LCMofDenInRow(IM,row));
          break;
        }
    return IM;
  }


/**
    expects a matrix with entries in an ordered ring
*/

  bool IsPositiveGrading(ConstMatrixView M)
  { return IsPositiveGrading(M, NumRows(M)); }
  

  bool IsPositiveGrading(ConstMatrixView M, long GradingDim)
  {
    if (!IsTermOrdering(M))  return false; // seems excessive but it isn't
    // check there are no null columns
    long i;
    for (long j=0; j < NumCols(M); ++j)
    {
      for (i=0; i < GradingDim; ++i) if (M(i,j) > 0) break;
      if (i==GradingDim) return false;
    }
    return true;
  }


  matrix NewMatMinimize(ConstMatrixView M)
  {
    if  (!IsZZ(RingOf(M)) && !IsQQ(RingOf(M)))
      CoCoA_ERROR("matrix must be over RingZZ() or RingQQ()", "NewMatMinimize");
    if ( NumRows(M)<NumCols(M) )
      CoCoA_ERROR(ERR::BadMatrixSize, "MatrixMinimize");
    if ( rank(M)!=NumCols(M) )
      CoCoA_ERROR(ERR::NotFullRank, "MatrixMinimize");

    matrix MM(NewDenseMat(RingOf(M), NumCols(M), NumCols(M)));
    long row=0;
    for (long r=0; r<NumRows(M); ++r)
    {
      for (long col=0; col<NumCols(M); ++col)
        SetEntry(MM, row, col, M(r, col));
      if (rank(MM) > row) ++row;
      if (row==NumRows(MM)) return MM;
    }
    CoCoA_ERROR(ERR::SERIOUS, "MatrixMinimize: should never get here");
    return MM;
  }


  matrix NewDenseMatXel(long n)
  {
    matrix RL(NewDenseMat(RingQQ(),n,n));
    for (long r=0; r<n; ++r)
      SetEntry(RL, r, n-r-1, 1);
    return RL;
  }


  matrix NewDenseMatRevLex(long n)
  {
    matrix RL(NewDenseMat(RingQQ(),n,n));
    for (long r=0; r<n; ++r)
      SetEntry(RL, r, n-r-1, -1);
    return RL;
  }


  matrix NewDenseMatStdDegLex(long n)
  {
    matrix RL(NewDenseMat(RingQQ(),n,n));
    for (long r=0; r<n; ++r)
      SetEntry(RL, 0, r, 1);
    for (long r=1; r<n; ++r)
      SetEntry(RL, r, r-1, 1);
    return RL;
  }


  matrix NewDenseMatStdDegRevLex(long n)
  {
    matrix RL(NewDenseMat(RingQQ(),n,n));
    for (long r=0; r<n; ++r)
      SetEntry(RL, 0, r, 1);
    for (long r=1; r<n; ++r)
      SetEntry(RL, r, n-r, -1);
    return RL;
  }


  class MatXelImpl: public MatrixViewBase
  {
  private:
    friend ConstMatrixView MatXel(long dim); // pseudo-ctor, uses RingQQ
    MatXelImpl(const ring& R, long dim);
    // default dtor is fine
  private: // disable default copy ctor and assignment
    MatXelImpl(const MatXelImpl&);            // NEVER DEFINED -- copy ctor disabled
    MatXelImpl& operator=(const MatXelImpl&); // NEVER DEFINED -- copy ctor disabled

  public:
    typedef std::vector<RingElem> vec;
    const ring& myRing() const;
    virtual long myNumRows() const;
    virtual long myNumCols() const;
    virtual RingElemAlias myEntry(long i, long j) const;
    virtual void myMulByRow(vec& lhs, const vec& v) const;
    virtual void myMulByCol(vec& lhs, const vec& v) const;
    virtual bool IamEqual(ConstMatrixView M) const;
    virtual bool IamSymmetric() const {return false;}
    virtual bool IamAntiSymmetric() const {return myNumRows()==0;}
    virtual bool IamDiagonal() const {return false;}
    virtual bool myIsZeroRow(long i) const;
    virtual bool myIsZeroCol(long j) const;
    virtual void myDet(RingElem& d) const;
    virtual long myRank() const;

    virtual void myAssignZero();
    virtual bool myIsWritable(long i, long j);
    virtual RingElemRawPtr myRawEntry(long i, long j);
    virtual void mySetEntry(long i, long j, ConstRefRingElem r);
    virtual void mySetEntry(long i, long j, const MachineInt& n);
    virtual void mySetEntry(long i, long j, const BigInt& N);
    virtual void mySetEntry(long i, long j, const BigRat& Q);

  private: // data members
    const ring myR;
    long myDim;
  };


//   /***************************************************************************/
//   // Functions for MatXelImpl

//   // bool IsMatXelImpl(ConstMatrixView M)
//   // {return dynamic_cast<const MatXelImpl*>(M.mySmartPtr.myRawPtr()) != 0;}


//   MatXelImpl::MatXelImpl(const ring& R, long dim):
//       MatrixViewBase(),
//       myR(R),
//       myDim(dim)
//   {}


//   const ring& MatXelImpl::myRing() const
//   {
//     return myR;
//   }


//   long MatXelImpl::myNumRows() const
//   {
//     return myDim;
//   }


//   long MatXelImpl::myNumCols() const
//   {
//     return myDim;
//   }


//   RingElemAlias MatXelImpl::myEntry(long i, long j) const
//   {
//     CoCoA_ASSERT(i >= 0 && j >= 0);
//     CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
//     if (i+j == myDim-1) return one(myR);
//     return zero(myR);
//   }


//   //BUG: should use FreeModule elems!
//   void MatXelImpl::myMulByRow(vec& lhs, const vec& v) const
//   {
//     // ASSUMES NO ALIASING
//     CoCoA_ASSERT(len(lhs) == myNumCols() && len(v) == myNumRows());
//     CoCoA_ASSERT(myNumRows() == 0 || myRing() == owner(v[0]));
//     CoCoA_ASSERT(myNumRows() == 0 || myRing() == owner(lhs[0]));
// //BUG???    CoCoA_ASSERT(myRing() == RingOf(v));
//     const long n = len(v);
//     for (long i=0; i < n; ++i)
//       lhs[i] = v[n-i-1];
//   }


//   void MatXelImpl::myMulByCol(vec& lhs, const vec& v) const
//   {
//     CoCoA_ASSERT(len(lhs) == myNumRows() && len(v) == myNumCols());
// //BUG???    CoCoA_ASSERT(myRing() == RingOf(v));
//     const long n = len(v);
//     for (long i=0; i < n; ++i)
//       lhs[i] = v[n-i-1];
//   }

//   // bool MatXelImpl::IamEqual(ConstMatrixView M) const
//   // {
//   //   if (myRing() != RingOf(M)) return false;
//   //   if (myNumRows() != NumRows(M)) return false;
//   //   if (myNumCols() != NumCols(M)) return false;
//   //   if (IsZeroMatImpl(M)) return NumCols(M) == 0;
//   //   if (IsMatXelImpl(M)) return true;
//   //   if (IsDiagMatImpl(M) || IsConstDiagMatImpl(M))
//   //   {
//   //     //      std::cout << "MatXelImpl::IamEqual - diag" << std::endl;
//   //     for (long i=0; i < myNumRows(); ++i) if (!IsOne(M(i,i))) return false;
//   //     return true;
//   //   }
//   //   return ConstMatrixViewBase::IamEqual(M);
//   // }


//   bool MatXelImpl::myIsZeroRow(long i) const
//   {
//     (void)(i); // to avoid compiler warning about unused parameter
//     CoCoA_ASSERT(0 <= i && i < myNumRows());
//     return false;
//   }


//   bool MatXelImpl::myIsZeroCol(long j) const
//   {
//     (void)(j); // to avoid compiler warning about unused parameter
//     CoCoA_ASSERT(0 <= j && j < myNumCols());
//     return false;
//   }


//   void MatXelImpl::myDet(RingElem& d) const
//   {
//     CoCoA_ASSERT(owner(d) == myRing());
//     if (IsEven(myDim/2))
//       d = 1;
//     else
//       d = -1;
//   }


//   long MatXelImpl::myRank() const
//   {
//     return myNumRows();
//   }


//   // Always gives error!
//   void MatXelImpl::myAssignZero()
//   {
// //???    if (myDim == 0) return;
//     CoCoA_ERROR(ERR::ConstMatEntry, "MatXel->myAssignZero()");
//   }


//   bool MatXelImpl::myIsWritable(long i, long j)
//   {
//     (void)(i); (void)(j); // to avoid compiler warning about unused parameter
//     CoCoA_ASSERT(i >= 0 && j >= 0);
//     CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
//     return false;
//   }


//   RingElemRawPtr MatXelImpl::myRawEntry(long /*i*/, long /*j*/)
//   {
//     CoCoA_ERROR(ERR::ConstMatEntry, "MatXel->myRawEntry(i,j)");
//     return RingElemRawPtr(0); // never reached, just to keep compiler quiet
//   }


//   void MatXelImpl::mySetEntry(long /*i*/, long /*j*/, ConstRefRingElem /*r*/)
//   {
//     CoCoA_ERROR(ERR::ConstMatEntry, "SetEntry(MatXel,i,j,RingElem)");
//   }

//   void MatXelImpl::mySetEntry(long /*i*/, long /*j*/, const MachineInt& /*n*/)
//   {
//     CoCoA_ERROR(ERR::ConstMatEntry, "SetEntry(MatXel,i,j,n)");
//   }

//   void MatXelImpl::mySetEntry(long /*i*/, long /*j*/, const BigInt& /*N*/)
//   {
//     CoCoA_ERROR(ERR::ConstMatEntry, "SetEntry(MatXel,i,j,BigInt)");
//   }


//   void MatXelImpl::mySetEntry(long /*i*/, long /*j*/, const BigRat& /*Q*/)
//   {
//     CoCoA_ERROR(ERR::ConstMatEntry, "SetEntry(MatXel,i,j,BigRat)");
//   }


//   // MatrixView MatXel(const ring& R, long dim)
//   // {
//   //   if (dim < 0)  CoCoA_ERROR(ERR::BadRowIndex, "MatXel(R,dim)");
//   //   return MatrixView(new MatXelImpl(R, dim));
//   // }

//   ConstMatrixView MatXel(long n)
//   {
//     return ConstMatrixView(new MatXelImpl(RingQQ(), n));
//   }



//   // matrix NewDenseMatRevLex(long n)
//   // {
//   //   matrix RL(NewDenseMat(RingQQ(),n,n));
//   //   for (long r=0; r<n; ++r)
//   //     SetEntry(RL, r, n-r-1, -1);
//   //   return RL;
//   // }


//   // matrix NewDenseMatStdDegLex(long n)
//   // {
//   //   matrix RL(NewDenseMat(RingQQ(),n,n));
//   //   for (long r=0; r<n; ++r)
//   //     SetEntry(RL, 0, r, 1);
//   //   for (long r=1; r<n; ++r)
//   //     SetEntry(RL, r, r-1, 1);
//   //   return RL;
//   // }


//   // matrix NewDenseMatStdDegRevLex(long n)
//   // {
//   //   matrix RL(NewDenseMat(RingQQ(),n,n));
//   //   for (long r=0; r<n; ++r)
//   //     SetEntry(RL, 0, r, 1);
//   //   for (long r=1; r<n; ++r)
//   //     SetEntry(RL, r, n-r, -1);
//   //   return RL;
//   // }


  matrix NewMatCompleteOrd(ConstMatrixView M)
  {
    matrix TORow=NewDenseMat(RingOf(M), 1, NumCols(M));
    for (long c=0; c<NumCols(M); ++c)
      if (IsZeroCol(M,c)) SetEntry(TORow, 0, c, 1);
    if (!IsTermOrdering(ConcatVer(M,TORow)))
      CoCoA_ERROR(ERR::NotTermOrdering, "NewMatCompleteOrd");
    return NewMatCompleteOrd(ConcatVer(M,TORow),NewDenseMatRevLex(NumCols(M)));
  }


  matrix NewMatCompleteOrd(ConstMatrixView M, ConstMatrixView NewRows)
  {
    if  (!IsZZ(RingOf(M)) && !IsQQ(RingOf(M)))
      CoCoA_ERROR("matrix must be over RingZZ() or RingQQ()", "NewMatCompleteOrd");
    if  (!IsZZ(RingOf(NewRows)) && !IsQQ(RingOf(NewRows)))
      CoCoA_ERROR("matrix must be over RingZZ() or RingQQ()", "NewMatCompleteOrd");
    if (RingOf(M) == RingOf(NewRows))
      return NewMatMinimize(NewDenseMat(ConcatVer(M, NewRows)));
    RingHom phi = ZZEmbeddingHom(RingQQ());
    if (IsZZ(RingOf(M)))
      return NewMatMinimize(NewDenseMat(ConcatVer(apply(phi,M), NewRows)));
    return NewMatMinimize(NewDenseMat(ConcatVer(M, apply(phi,NewRows))));
  }

  //--- matrices for elimination -----------------------------

  namespace // anonymous namespace for file local auxiliary funcs and defs
  {
    
    matrix ElimRow(const std::vector<long>& IndetsToElim, long NumIndets)
    {
      matrix M = NewDenseMat(RingQQ(), 1, NumIndets);
      long s = len(IndetsToElim);
      for (long i=0; i < s; ++i)
      {
        if (IndetsToElim[i]<0 || IndetsToElim[i]>=NumIndets)
          CoCoA_ERROR(ERR::BadIndetIndex, "ElimRow");
        SetEntry(M,  0, IndetsToElim[i],  1);
      }
      return M;
    }
    
    
    matrix AuxElimAndM(const vector<long>& IndetsToElim, ConstMatrixView M)
    {
      if  (!IsZZ(RingOf(M)) && !IsQQ(RingOf(M)))
        CoCoA_ERROR("matrix must be over RingZZ() or RingQQ()", "NewMatAElim");
      if (IsQQ(RingOf(M)))
        return NewMatCompleteOrd(ConcatVer(ElimRow(IndetsToElim,NumCols(M)),M));
      return NewMatCompleteOrd(ConcatVer(ElimRow(IndetsToElim, NumCols(M)),
                                         apply(ZZEmbeddingHom(RingQQ()),M)));
    }
    

    matrix AuxMAndElim(ConstMatrixView M, const vector<long>& IndetsToElim)
    {
      if  (!IsZZ(RingOf(M)) && !IsQQ(RingOf(M)))
        CoCoA_ERROR("matrix must be over RingZZ() or RingQQ()", "NewMatAElim");
      if (IsQQ(RingOf(M)))
        return NewMatCompleteOrd(ConcatVer(M, ElimRow(IndetsToElim,NumCols(M))));
      return NewMatCompleteOrd(ConcatVer(apply(ZZEmbeddingHom(RingQQ()),M),
                                         ElimRow(IndetsToElim, NumCols(M)) ));
    }

  }
  
  
  matrix ElimMat(long NumIndets, const std::vector<long>& IndetsToElim)
  {
    // return AuxElimAndM(IndetsToElim, FilledMat(RingQQ(),1,NumIndets, 1));
    return AuxElimAndM(IndetsToElim,
                       RowMat(vector<RingElem>(NumIndets, one(RingQQ()))));
  }


  matrix ElimMat(ConstMatrixView GradM,
                 const std::vector<long>& IndetsToElim)
  {
    if (NumRows(GradM)!=0) return AuxElimAndM(IndetsToElim, GradM);
    return ElimMat(NumCols(GradM), IndetsToElim); // with row of 1's
  }


  matrix HomogElimMat(ConstMatrixView GradM,
                      const std::vector<long>& IndetsToElim)
  {
    if (NumRows(GradM)==0)
      CoCoA_ERROR(ERR::ZeroGradingDim, "HomogElimMat");
    return AuxMAndElim(GradM, IndetsToElim);
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/MatrixForOrdering.C,v 1.21 2014/07/30 14:06:06 abbott Exp $
// $Log: MatrixForOrdering.C,v $
// Revision 1.21  2014/07/30 14:06:06  abbott
// Summary: Changed BaseRing into RingOf
// Author: JAA
//
// Revision 1.20  2014/04/17 13:38:41  bigatti
// -- MatrixViews --> MatrixView
//
// Revision 1.19  2014/04/11 15:44:27  abbott
// Summary: Renamed MatrixArith to MatrixOperations (in includes)
// Author: JAA
//
// Revision 1.18  2014/04/08 12:56:12  bigatti
// -- cambiata ElimMat per non usare FilledMat
//
// Revision 1.17  2013/02/14 17:34:35  bigatti
// -- cleaned up code for elimination matrices
//
// Revision 1.16  2012/11/21 09:46:33  bigatti
// -- NewMatCompleteOrd(M) now returns a term-ordering (if M suitable)
//
// Revision 1.15  2012/03/30 17:28:09  bigatti
// -- added NewIntegerOrdMat
// -- accepting and returning matrices over QQ
//
// Revision 1.14  2012/02/10 10:28:08  bigatti
// -- changed RingZ.H, RingQ.H --> RingZZ.H, RingQQ.H
//
// Revision 1.13  2012/02/08 17:20:43  bigatti
// -- changed: Z,Q -> ZZ,QQ
// -- code reorganization
//
// Revision 1.12  2011/05/26 11:58:05  bigatti
// -- added IsPositiveGrading with one arg
//
// Revision 1.11  2011/04/26 10:11:06  bigatti
// -- added NewMatCompleteOrd
//
// Revision 1.10  2011/03/23 17:29:54  bigatti
// -- added NewDenseMatStdDegLex
//
// Revision 1.9  2011/03/21 07:50:51  bigatti
// -- added NewDenseMatXel, NewDenseMatStdDegRevLex
//
// Revision 1.8  2011/03/10 11:26:26  bigatti
// -- using len(v) instead of v.size()
//
// Revision 1.7  2011/03/09 09:08:25  bigatti
// -- removing signed/unsigned warnings
// -- added samity checks for indets to elim
//
// Revision 1.6  2011/03/08 16:10:16  abbott
// Changed size_t into long.
//
// Revision 1.5  2011/02/16 14:58:40  bigatti
// -- fixed typo
//
// Revision 1.4  2009/12/03 17:40:36  abbott
// Added include directives for ZZ.H (as a consequence of removing
// the directive from ring.H).
//
// Revision 1.3  2009/09/22 13:35:55  bigatti
// -- following coding conventions in function names Matrix --> Mat
// -- forced all matrices to be over RingZZ
//
// Revision 1.2  2008/05/30 12:44:14  abbott
// Moved "ordering matrices" into their ownn special file.
//
// Revision 1.1  2008/04/21 11:23:11  abbott
// Separated functions dealing with matrices and PPOrderings into a new file.
// Added matrix norms, and completed adjoint.
//
//
