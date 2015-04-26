//   Copyright (c)  2004-2011  John Abbott

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


#include "CoCoA/PPOrdering.H"
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/MatrixOperations.H"
#include "CoCoA/MatrixForOrdering.H"
#include "CoCoA/OpenMath.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/BigRat.H"
#include "CoCoA/assert.H"
#include "CoCoA/error.H"
#include "CoCoA/matrix.H"
#include "CoCoA/MatrixView.H"  // for submat
#include "CoCoA/utils.H"  // for LongRange

#include <iostream>
using std::ostream;
//#include <vector>
using std::vector;

namespace CoCoA
{

  std::ostream& operator<<(std::ostream& out, const PPOrdering& PPO)
  {
    PPO->myOutputSelf(out);
    return out;
  }


  OpenMathOutput& operator<<(OpenMathOutput& OMOut, const PPOrdering& PPO)
  {
    PPO->myOutputSelf(OMOut);
    return OMOut;
  }


  //---------------------------------------------------------------------------

  PPOrderingBase::PPOrderingBase(long NumIndets, long GradingDim):
    IntrusiveReferenceCount(),
    myNumIndets(NumIndets),
    myGradingDim(GradingDim)
  {
    CoCoA_ASSERT(NumIndets > 0);
    CoCoA_ASSERT(0 <= GradingDim && GradingDim <= NumIndets);
  }


  PPOrderingBase::~PPOrderingBase()
  {}



  //---------------------------------------------------------------------------

  // ??? Compiler bug: the friend declarations below need explicit namespace
  //     scoping with g++ 3.2.2, otherwise g++ thinks NewLexOrdering etc are
  //     in the CoCoA::PPO namespace.  Very odd!

  namespace PPOrd
  {
PPOrdering lex(long NumIndets){return NewLexOrdering(NumIndets);}
PPOrdering StdDegLex(long NumIndets){return NewStdDegLexOrdering(NumIndets);}
PPOrdering StdDegRevLex(long NumIndets){return NewStdDegRevLexOrdering(NumIndets);}

    class LexImpl: public PPOrderingBase
    {
    private:
      friend PPOrdering CoCoA::NewLexOrdering(long NumIndets);
      LexImpl(long NumIndets);
      LexImpl(const LexImpl&);            ///< NEVER DEFINED -- copy ctor disabled
      LexImpl& operator=(const LexImpl&); ///< NEVER DEFINED -- assignment disabled
    public:
      virtual void myOutputSelf(std::ostream& out) const;
      virtual void myOutputSelf(OpenMathOutput& OMOut) const;
      virtual void myOrdMatCopy(MatrixView& M) const;
      virtual bool IamStdGraded() const {return false;}
    };


    class StdDegLexImpl: public PPOrderingBase
    {
    private:
      friend PPOrdering CoCoA::NewStdDegLexOrdering(long NumIndets);
      StdDegLexImpl(long NumIndets);
      StdDegLexImpl(const StdDegLexImpl&);            ///< NEVER DEFINED -- copy ctor disabled
      StdDegLexImpl& operator=(const StdDegLexImpl&); ///< NEVER DEFINED -- assignment disabled
    public:
      virtual void myOutputSelf(std::ostream& out) const;
      virtual void myOutputSelf(OpenMathOutput& OMOut) const;
      virtual void myOrdMatCopy(MatrixView& M) const;
      virtual bool IamStdGraded() const {return true;}
    };


    class StdDegRevLexImpl: public PPOrderingBase
    {
    private:
      friend PPOrdering CoCoA::NewStdDegRevLexOrdering(long NumIndets);
      StdDegRevLexImpl(long NumIndets);
      StdDegRevLexImpl(const StdDegRevLexImpl&);            ///< NEVER DEFINED -- copy ctor disabled
      StdDegRevLexImpl& operator=(const StdDegRevLexImpl&); ///< NEVER DEFINED -- assignment disabled
    public:
      virtual void myOutputSelf(std::ostream& out) const;
      virtual void myOutputSelf(OpenMathOutput& OMOut) const;
      virtual void myOrdMatCopy(MatrixView& M) const;
      virtual bool IamStdGraded() const {return true;}
    };


    class MatrixOrderingImpl: public PPOrderingBase
    {
    private:
      friend PPOrdering CoCoA::NewMatrixOrdering(long NumIndets, long GradingDim, ConstMatrixView OrderMatrix);
      MatrixOrderingImpl(long NumIndets, long GradingDim, ConstMatrixView OrderMatrix);
      MatrixOrderingImpl(const MatrixOrderingImpl&);            ///< NEVER DEFINED -- copy ctor disabled
      MatrixOrderingImpl& operator=(const MatrixOrderingImpl&); ///< NEVER DEFINED -- assignment disabled
    private: ///< data members (in addition to those inherited)
      matrix myDefiningMatrix;
    public:
      virtual void myOutputSelf(std::ostream& out) const;
      virtual void myOutputSelf(OpenMathOutput& OMOut) const;
      virtual void myOrdMatCopy(MatrixView& M) const;
      virtual bool IamStdGraded() const;
    };


  } // end of namespace PPO


  //----------------------------------------------------------------------
  // Here are the pseudo-constructors:

  PPOrdering NewLexOrdering(long NumIndets)
  {
    if (NumIndets <= 0) CoCoA_ERROR(ERR::BadArg, "NewLexOrdering(n)");
    return PPOrdering(new PPOrd::LexImpl(NumIndets));
  }


  PPOrdering NewStdDegLexOrdering(long NumIndets)
  {
    if (NumIndets <= 0) CoCoA_ERROR(ERR::BadArg, "NewStdDegLexOrdering(n)");
    return PPOrdering(new PPOrd::StdDegLexImpl(NumIndets));
  }


  PPOrdering NewStdDegRevLexOrdering(long NumIndets)
  {
    if (NumIndets <= 0) CoCoA_ERROR(ERR::BadArg, "NewStdDegRevLexOrdering(n)");
    return PPOrdering(new PPOrd::StdDegRevLexImpl(NumIndets));
  }


  PPOrdering NewMatrixOrdering(long NumIndets, long GradingDim, ConstMatrixView OrderMatrix)
  {
    if (NumIndets <= 0 || GradingDim < 0 || GradingDim > NumRows(OrderMatrix))
      CoCoA_ERROR(ERR::BadArg, "NewMatrixOrdering(n,gd,M)");
    return PPOrdering(new PPOrd::MatrixOrderingImpl(NumIndets, GradingDim, OrderMatrix));
  }


//   PPOrdering NewMatrixOrderingMod32749(long NumIndets, long GradingDim, ConstMatrix OrderMatrix, const matrix& InverseOrderMatrixMod32749)
//   {
//    if (NumIndets <= 0 || GradingDim < 0 || GradingDim > NumRows(OrderMatrix))
//       CoCoA_ERROR(ERR::BadArg, "NewMatrixOrdering(n,gd,M)"); ///??? IsZero(NumIndets)???
//     return PPOrdering(new PPOrd::MatrixOrderingMod32749Impl(NumIndets, GradingDim, OrderMatrix, InverseOrderMatrixMod32749));
//   }




  namespace PPOrd
  {

    LexImpl::LexImpl(long NumIndets):
      PPOrderingBase(NumIndets, 0)
    {}


    void LexImpl::myOutputSelf(std::ostream& out) const
    {
      out << "PPOrderingLex(" << myNumIndets << ")";
    }


    void LexImpl::myOutputSelf(OpenMathOutput& OMOut) const
    {
      OMOut->mySendApplyStart();
      OMOut << OpenMathSymbol("cocoa", "PPOrderingLex");
      OMOut << myNumIndets;
      OMOut->mySendApplyEnd();
    }


    void LexImpl::myOrdMatCopy(MatrixView& M) const
    {
      CoCoA_ASSERT(NumRows(M) == myNumIndets);
      CoCoA_ASSERT(NumCols(M) == myNumIndets);
      //??? There should be a better way than this!
      AssignZero(M);
      for (long i=0; i < myNumIndets; ++i)
        SetEntry(M, i, i, 1);
    }




    StdDegLexImpl::StdDegLexImpl(long NumIndets):
      PPOrderingBase(NumIndets, 1)
    {}


    void StdDegLexImpl::myOutputSelf(std::ostream& out) const
    {
      out << "PPOrderingStdDegLex(" << myNumIndets << ")";
    }


    void StdDegLexImpl::myOutputSelf(OpenMathOutput& OMOut) const
    {
      OMOut->mySendApplyStart();
      OMOut << OpenMathSymbol("cocoa", "PPOrderingStdDegLex");
      OMOut << myNumIndets;
      OMOut->mySendApplyEnd();
    }


    void StdDegLexImpl::myOrdMatCopy(MatrixView& M) const
    {
      CoCoA_ASSERT(NumRows(M) == myNumIndets);
      CoCoA_ASSERT(NumCols(M) == myNumIndets);

      AssignZero(M);
      for (long i=0; i < myNumIndets; ++i)
      {
        SetEntry(M, 0, i, 1);
        if (i > 0)
          SetEntry(M, i, i-1, 1);
      }
    }


    //------------------------------------------------------------//


    StdDegRevLexImpl::StdDegRevLexImpl(long NumIndets):
      PPOrderingBase(NumIndets, 1)
    {}


    void StdDegRevLexImpl::myOutputSelf(std::ostream& out) const
    {
      out << "PPOrderingStdDegRevLex(" << myNumIndets << ")";
    }


    void StdDegRevLexImpl::myOutputSelf(OpenMathOutput& OMOut) const
    {
      OMOut->mySendApplyStart();
      OMOut << OpenMathSymbol("cocoa", "PPOrderingStdDegRevLex");
      OMOut << myNumIndets;
      OMOut->mySendApplyEnd();
    }


    void StdDegRevLexImpl::myOrdMatCopy(MatrixView& M) const
    {
      CoCoA_ASSERT(NumRows(M) == myNumIndets);
      CoCoA_ASSERT(NumCols(M) == myNumIndets);

      AssignZero(M);
      for (long i=0; i < myNumIndets; ++i)
      {
        SetEntry(M, 0, i, 1);
        if (i > 0)
          SetEntry(M, i, myNumIndets-i, -1);
      }
    }


    //------------------------------------------------------------//


    MatrixOrderingImpl::MatrixOrderingImpl(long NumIndets, long GradingDim, ConstMatrixView OrderMatrix):
      PPOrderingBase(NumIndets, GradingDim),
      myDefiningMatrix(NewDenseMat(OrderMatrix))
    {
      const char* const fn = "MatrixOrderingImpl ctor";
      if (NumCols(OrderMatrix) != NumIndets)
        CoCoA_ERROR(ERR::BadMatrixSize, fn);
      // ??? ANNA: for the time being only square matrices are accepted.
      //           Later may allow more rows, provided rank = NumIndets.
      if (NumRows(OrderMatrix) != NumIndets)
        CoCoA_ERROR(ERR::BadMatrixSize, fn);
      if (GradingDim > NumRows(OrderMatrix))
        CoCoA_ERROR(ERR::BadRowIndex, fn);
      if (NumIndets == 0) return;

      if (rank(OrderMatrix) != NumIndets)
        CoCoA_ERROR(ERR::NotFullRank, fn);
    }


    void MatrixOrderingImpl::myOutputSelf(std::ostream& out) const
    {
      out << "PPOrderingMatrix(GrDim=" << myGradingDim << ", " << myDefiningMatrix << ")";
    }


    void MatrixOrderingImpl::myOutputSelf(OpenMathOutput& OMOut) const
    {
      OMOut->mySendApplyStart();
      OMOut << OpenMathSymbol("cocoa", "PPOrderingMatrix");
      OMOut << myNumIndets;
      OMOut << myGradingDim;
      OMOut->mySendApplyEnd();
    }


    void MatrixOrderingImpl::myOrdMatCopy(MatrixView& M) const
    {
      CoCoA_ASSERT(NumRows(M) == myNumIndets);
      CoCoA_ASSERT(NumCols(M) == myNumIndets);

      BigRat tmp;
      for (long i=0; i < myNumIndets; ++i)
        for (long j=0; j < myNumIndets; ++j)
        {
          IsRational(tmp, myDefiningMatrix(i,j));
          SetEntry(M, i, j, tmp);
        }
    }


    bool MatrixOrderingImpl::IamStdGraded() const
    {
      if (myGradingDim!=1) return false;
      for (long j=0; j < myNumIndets; ++j)
        if (myDefiningMatrix(0,j)!=1) return false;
      return true;
    }
    

    //------------------------------------------------------------//


  } // end of namespace PPO


  bool IsLex(const PPOrdering& PPO)
  {
    if (dynamic_cast<const PPOrd::LexImpl*>(PPO.myRawPtr())) return true;
    if (!dynamic_cast<const PPOrd::MatrixOrderingImpl*>(PPO.myRawPtr())) return false;
    // must decide whether the matrix is Lex, possibly in disguise
    return false; // just to keep compiler quiet for the moment
  }


  bool IsStdDegLex(const PPOrdering& PPO)
  {
    if (dynamic_cast<const PPOrd::StdDegLexImpl*>(PPO.myRawPtr())) return true;
    if (!dynamic_cast<const PPOrd::MatrixOrderingImpl*>(PPO.myRawPtr())) return false;
    // must decide whether the matrix is StdDegLex, possibly in disguise
    return false; // just to keep compiler quiet for the moment
  }


  bool IsStdDegRevLex(const PPOrdering& PPO)
  {
    if (dynamic_cast<const PPOrd::StdDegRevLexImpl*>(PPO.myRawPtr())) return true;
    if (!dynamic_cast<const PPOrd::MatrixOrderingImpl*>(PPO.myRawPtr())) return false;
    // must decide whether the matrix is StdDegRevLex, possibly in disguise
    return false; // just to keep compiler quiet for the moment
  }


  bool IsMatrixOrdering(const PPOrdering& PPO)
  {
    return (dynamic_cast<const PPOrd::MatrixOrderingImpl*>(PPO.myRawPtr())); // nullptr <--> false
  }


  bool IsTermOrdering(const PPOrdering& PPO)
  {
    return IsLex(PPO) ||
           IsStdDegLex(PPO) ||
           IsStdDegRevLex(PPO) ||
           IsTermOrdering(OrdMat(PPO));
  }


  matrix OrdMat(const PPOrdering& PPO)
  {
    // To which ring should the elements belong???
    matrix ans(NewDenseMat(RingQQ(), NumIndets(PPO), NumIndets(PPO)));
    PPO->myOrdMatCopy(ans);
    return ans;
  }


  matrix GradingMat(const PPOrdering& PPO)
  {
    matrix M(NewDenseMat(OrdMat(PPO)));
    return NewDenseMat(submat(M,
                              LongRange(0,GradingDim(PPO)-1),
                              LongRange(0,NumCols(M)-1)));
  }


  LexCtor lex;
  StdDegLexCtor StdDegLex;
  StdDegRevLexCtor StdDegRevLex;


} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/PPOrdering.C,v 1.16 2014/07/31 13:10:45 bigatti Exp $
// $Log: PPOrdering.C,v $
// Revision 1.16  2014/07/31 13:10:45  bigatti
// -- GetMatrix(PPO) --> OrdMat(PPO)
// -- added OrdMat and GradingMat to PPOrdering, PPMonoid, SparsePolyRing
//
// Revision 1.15  2014/04/11 15:44:27  abbott
// Summary: Renamed MatrixArith to MatrixOperations (in includes)
// Author: JAA
//
// Revision 1.14  2014/01/28 16:44:40  abbott
// Added new fns:  IsMatrixOrdering  and  IsTermOrdering
//
// Revision 1.13  2013/07/30 14:57:45  bigatti
// -- added IamStdGraded
//
// Revision 1.12  2012/03/30 17:29:56  bigatti
// -- accepting and returning matrices over QQ
//
// Revision 1.11  2012/02/10 10:28:08  bigatti
// -- changed RingZ.H, RingQ.H --> RingZZ.H, RingQQ.H
//
// Revision 1.10  2012/02/08 16:13:41  bigatti
// -- changed Z,Q --> ZZ,QQ
//
// Revision 1.9  2011/08/14 15:52:17  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.8  2011/03/10 16:39:34  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.7  2011/03/04 16:37:55  bigatti
// -- changed: matrix members functions args of type long instead of size_t
//
// Revision 1.6  2010/02/03 16:13:52  abbott
// Added new single word tags for specifying the ordering in PPMonoid
// pseudo-ctors.
//
// Revision 1.5  2009/07/30 15:44:53  bigatti
// -- added comment for pseudo-constructors
//
// Revision 1.4  2008/04/18 15:35:57  abbott
// (long overdue) Major revision to matrices
//
// Revision 1.3  2008/04/08 15:26:42  abbott
// Major revision to matrix implementation: added matrix views.
// Lots of changes.
//
// Revision 1.2  2007/10/30 17:14:07  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.11  2007/03/08 18:22:29  cocoa
// Just whitespace cleaning.
//
// Revision 1.10  2007/03/07 14:06:53  bigatti
// -- CoCoA_ASSERT converted into CoCoA_ERROR
//
// Revision 1.9  2007/03/02 10:47:53  cocoa
// First stage of RingZ modifications -- tests do not compile currently, Anna will fix this.
//
// Revision 1.8  2007/01/18 17:05:22  bigatti
// -- changed namespace PPO into PPOrd (as for ModuleTermOrd)
//
// Revision 1.7  2006/11/27 14:25:15  cocoa
// -- reorganised #include files
//
// Revision 1.6  2006/11/22 14:50:33  cocoa
// -- changed: PPOrdering defined as class (instead of typedef)
//
// Revision 1.5  2006/11/16 11:27:20  cocoa
// -- reinserted myRefCountZero(): sometimes really necessary, in general safe
//
// Revision 1.4  2006/11/14 17:38:47  cocoa
// -- deleted commented out code about reference counting
// -- commented out myRefCountZero() (not necessary?)
//
// Revision 1.3  2006/11/02 13:25:44  cocoa
// Simplification of header files: the OpenMath classes have been renamed.
// Many minor consequential changes.
//
// Revision 1.2  2006/10/06 14:04:14  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.6  2006/05/12 16:10:58  cocoa
// Added OpenMathFwd.H, and tidied OpenMath.H.
// Many consequential but trivial changes.
//
// Revision 1.5  2006/04/28 16:33:51  cocoa
// Used SmartPtrIRC for PPOrderings.
//
// Revision 1.4  2006/03/15 18:09:31  cocoa
// Changed names of member functions which print out their object
// into myOutputSelf -- hope this will appease the Intel C++ compiler.
//
// Revision 1.3  2006/03/12 21:28:33  cocoa
// Major check in after many changes
//
// Revision 1.2  2005/12/16 11:24:21  cocoa
// -- fixed  StdDegRevLexImpl::myOrdMatCopy  [1 --> -1]
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.3  2005/06/22 14:47:56  cocoa
// PPMonoids and PPMonoidElems updated to mirror the structure
// used for rings and RingElems.  Many consequential changes.
//
// Revision 1.2  2005/05/04 16:51:47  cocoa
// -- changed: check on matrix ordering computes "rank" (was "det")
//
// Revision 1.1.1.1  2005/05/03 15:47:31  cocoa
// Imported files
//
// Revision 1.4  2005/04/19 15:39:55  cocoa
// Matrices now use reference counts.
//
// Revision 1.3  2005/04/19 14:06:04  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.2  2005/03/02 18:46:41  cocoa
// Added new types ConstRefMatrix, and RefMatrix following along
// the lines of ConstRefRingElem and RefRingElem.  The semantics
// should be a bit clearer now.
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.14  2004/11/29 16:22:35  cocoa
// -- added function for computing adjoint and inverse for DenseMatrix
//    (so adjoint/inverse matrix is computed by OrdvArith and is no
//    longer needed by PPOrdering)
//
// Revision 1.13  2004/11/12 15:49:29  cocoa
// Tidying prior to 0.90 release.
// (a) documentation improved (or marked as poor)
// (b) sundry minor improvements to the code
//
// Revision 1.12  2004/11/11 14:42:01  cocoa
// -- change: cout --> GlobalLogput()
//
// Revision 1.11  2004/11/05 16:44:20  cocoa
// -- deleted MatrixOrderingMod32749Impl (implemented in OrdvArith)
// -- changed C++ matrices into "matrix" over RingZ
//
// Revision 1.10  2004/11/03 17:54:44  cocoa
// -- added implementation of GetMatrix (OrdMat)
// -- added some functions for order matrices modulo 32749:
//    they will be deleted soon
//
// Revision 1.9  2004/11/02 14:49:03  cocoa
// -- new code for matrix orderings
//
// Revision 1.8  2004/10/29 16:09:21  cocoa
// -- added MatrixOrderingMod32749Impl (not tested)
// -- fixed assignment of myAdjointMatrix for MatrixOrderingImpl
//
// Revision 1.7  2004/10/21 17:16:37  cocoa
// Fairly major change: new OrdvArith namspace with various members,
//   new global typedef  SmallExponent_t (defined in config.H).
//
