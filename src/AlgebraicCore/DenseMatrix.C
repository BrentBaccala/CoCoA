//   Copyright (c)  2005-2009,2011  John Abbott

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


// Source code for classes DenseMatBase and DenseMatImpl

#include "CoCoA/DenseMatrix.H"
//#include "CoCoA/MachineInt.H"
#include "CoCoA/BigInt.H"
#include "CoCoA/BigRat.H"
#include "CoCoA/MatrixOperations.H"
#include "CoCoA/assert.H"
#include "CoCoA/error.H"
#include "CoCoA/matrix.H"
#include "CoCoA/ring.H"
#include "CoCoA/utils.H" // for len


//#include <vector>
using std::vector;


namespace CoCoA
{

  class DenseMatBase: public MatrixBase
  {
  protected: // data members
    const ring myR;
    long myNumRowsValue;
    long myNumColsValue;
  protected:
    DenseMatBase(const ring& R, long nrows, long ncols);
  };



  class DenseMatImpl: public DenseMatBase
  {
  public:
    DenseMatImpl(const ring& R, long nrows, long ncols); // creates zero matrix
    ~DenseMatImpl();
  private: // implementation detail
    void myDtorBody();

  public:
    typedef std::vector<RingElem> vec;
    const ring& myRing() const;
    virtual long myNumRows() const;
    virtual long myNumCols() const;
    virtual RingElemAlias myEntry(long i, long j) const;
    virtual DenseMatImpl* myClone() const;
    virtual DenseMatImpl* myZeroClone(const ring& R, long NumRows, long NumCols) const;

    virtual void myAssignZero(); // sets all entries to zero
    bool myIsWritable(long /*i*/, long /*j*/) { return true; }
    virtual RingElemRawPtr myRawEntry(long i, long j);
    virtual void mySetEntry(long i, long j, ConstRefRingElem r);
    virtual void mySetEntry(long i, long j, const MachineInt& n);
    virtual void mySetEntry(long i, long j, const BigInt& N);
    virtual void mySetEntry(long i, long j, const BigRat& Q);
// Use generic versions of the following four functions.
//     virtual void myMulByRow(vec& lhs, const vec& v) const;
//     virtual void myMulByCol(vec& lhs, const vec& v) const;
//     virtual bool myIsZeroRow(long i) const; ///< tests whether row i is zero
//     virtual bool myIsZeroCol(long j) const; ///< tests whether col j is zero
    virtual void myRowMul(long i, ConstRefRingElem c); ///< row(i) = c*row(i)
    virtual void myColMul(long j, ConstRefRingElem c); ///< col(j) = c*col(j)
    virtual void myAddRowMul(long i1, long i2, ConstRefRingElem c); ///< row(i1) += c*row(i2)
    virtual void myAddColMul(long j1, long j2, ConstRefRingElem c); ///< col(j1) += c*col(j2)
    virtual void mySwapRows(long i1, long i2);
    virtual void mySwapCols(long j1, long j2);
    virtual void myResize(long NumRows, long NumCols);
    virtual void myDet(RingElem& d) const;
    virtual long myRank() const;

   private: // data members
    vector< vector< RingElemRawPtr > > myEntries;
  };



  /////////////////////////////////////////////////////////////////////////////
  // Inline functions

  inline DenseMatBase::DenseMatBase(const ring& R, long NumRows, long NumCols):
      myR(R),
      myNumRowsValue(NumRows),
      myNumColsValue(NumCols)
  {}


//----------------------------------------------------------------------

  DenseMatImpl::DenseMatImpl(const ring& R, long NumRows, long NumCols):
      DenseMatBase(R, NumRows, NumCols)
  {
    // Making this exception safe turns out to be a bit complicated.
    try
    {
      myEntries.reserve(myNumRows()); // not necessary but might help
      for (long i=0; i < myNumRows(); ++i)
      {
        myEntries.push_back(vector<RingElemRawPtr>());
        myEntries[i].reserve(myNumCols()); // not necessary but might help
        for (long j=0; j < myNumCols(); ++j)
        {
          myEntries[i].push_back(myR->myNew());
        }
      }
    }
    catch (...)
    {
      myDtorBody();
    }
  }


  DenseMatImpl::~DenseMatImpl()
  {
    myDtorBody();
  }


  // Ignore the values of myNumRows and myNumCols because this procedure
  // may be called on a partially constructed DenseMat if an exception
  // is thrown during the ctor.
  void DenseMatImpl::myDtorBody()
  {
    // Written this way, we destroy in reverse order of construction.
    // NB:  myEntries.size()-1 is **unsigned**
    for (long i=len(myEntries)-1; i >= 0; --i)
      for (long j=len(myEntries[i])-1; j >= 0; --j)
        myR->myDelete(myEntries[i][j]);
  }


  const ring& DenseMatImpl::myRing() const
  {
    return myR;
  }


  long DenseMatImpl::myNumRows() const
  {
    return myNumRowsValue;
  }


  long DenseMatImpl::myNumCols() const
  {
    return myNumColsValue;
  }


  inline RingElemAlias DenseMatImpl::myEntry(long i, long j) const
  {
    CoCoA_ASSERT(i < myNumRows());
    CoCoA_ASSERT(j < myNumCols());
    return RingElemAlias(myR, myEntries[i][j]);
  }


  DenseMatImpl* DenseMatImpl::myClone() const
  {
    // Making this exception safe turns out to be a bit complicated.
    DenseMatImpl* ans = 0/*nullptr*/;
    try
    {
      ans = new DenseMatImpl(myR, myNumRows(), myNumCols());
      for (long i=0; i < myNumRows(); ++i)
      {
        for (long j=0; j < myNumCols(); ++j)
        {
          myR->myAssign(ans->myEntries[i][j],myEntries[i][j]); // could throw
        }
      }
    }
    catch (...)
    {
      ans->myDtorBody();
      throw;
    }
    return ans;
  }


  DenseMatImpl* DenseMatImpl::myZeroClone(const ring& R, long NumRows, long NumCols) const
  {
    return new DenseMatImpl(R, NumRows, NumCols);
  }


  void DenseMatImpl::myAssignZero()
  {
    // This fn is probably not exception safe -- can RingBase::myAssignZero throw?
    for (long i=0; i < myNumRows(); ++i)
      for (long j=0; j < myNumCols(); ++j)
        myR->myAssignZero(myEntries[i][j]);
  }


  RingElemRawPtr DenseMatImpl::myRawEntry(long i, long j)
  {
    CoCoA_ASSERT(i < myNumRows());
    CoCoA_ASSERT(j < myNumCols());
    return myEntries[i][j];
  }


  void DenseMatImpl::mySetEntry(long i, long j, ConstRefRingElem r)
  {
    CoCoA_ASSERT(owner(r) == myR);
    CoCoA_ASSERT(i < myNumRows());
    CoCoA_ASSERT(j < myNumCols());
    myR->myAssign(myEntries[i][j], raw(r));
  }


  void DenseMatImpl::mySetEntry(long i, long j, const MachineInt& n)
  {
    CoCoA_ASSERT(i < myNumRows());
    CoCoA_ASSERT(j < myNumCols());
    myR->myAssign(myEntries[i][j], n);
  }


  void DenseMatImpl::mySetEntry(long i, long j, const BigInt& N)
  {
    CoCoA_ASSERT(i < myNumRows());
    CoCoA_ASSERT(j < myNumCols());
    myR->myAssign(myEntries[i][j], N);
  }


  void DenseMatImpl::mySetEntry(long i, long j, const BigRat& Q)
  {
    CoCoA_ASSERT(i < myNumRows());
    CoCoA_ASSERT(j < myNumCols());
    myR->myAssign(myEntries[i][j], Q);
  }


//   MOVED UP TO GENERIC IMPLEMENTATIONS IN ConstMatrixBase
//   void DenseMatImpl::myMulByRow(vec& lhs, const vec& v) const
//   {
//     CoCoA_ASSERT(lhs.size() == myNumCols() && v.size() == myNumRows());
//     vector<RingElem> ans;
//     ans.resize(myNumCols(), zero(myR));
//     for (long j=0; j < myNumCols(); ++j)
//       for (long i=0; i < myNumRows(); ++i)
//         ans[j] += v[i]*myEntry(i, j);
//     // We have successfully computed the answer, now swap it in.
//     swap(lhs, ans);
// //     for (long j=0; j < myNumCols(); ++j)
// //       swap(lhs[j], ans[j]);
//   }


//   void DenseMatImpl::myMulByCol(vec& lhs, const vec& v) const
//   {
//     CoCoA_ASSERT(lhs.size() == myNumRows() && v.size() == myNumCols());
//     vector<RingElem> ans;
//     ans.resize(myNumRows(), zero(myR));
//     for (long i=0; i < myNumRows(); ++i)
//       for (long j=0; j < myNumCols(); ++j)
//         ans[i] += myEntry(i, j)*v[j];

//     // We have successfully computed the answer, now swap it in.
//     swap(lhs, ans);
// //     for (long i=0; i < myNumRows(); ++i)
// //       swap(lhs[i], ans[i]);
//   }


//   bool DenseMatImpl::myIsZeroRow(long i) const
//   {
//     if (i >= myNumRows())
//       CoCoA_ERROR(ERR::BadRowIndex, "DenseMatImpl::myIsZeroRow");
//     for (long j=0; j < myNumCols(); ++j)
//       if (!myR->myIsZero(myEntries[i][j])) return false;
//     return true;
//   }


//   bool DenseMatImpl::myIsZeroCol(long j) const
//   {
//     if (j >= myNumCols())
//       CoCoA_ERROR(ERR::BadColIndex, "DenseMatImpl::myIsZeroCol");
//     for (long i=0; i < myNumRows(); ++i)
//       if (!myR->myIsZero(myEntries[i][j])) return false;
//     return true;
//   }


  void DenseMatImpl::myRowMul(long i, ConstRefRingElem c)
  {
    const char* const FnName = "DenseMatImpl::myRowMul";
    myCheckRowIndex(i, FnName);
    if (owner(c) != myR)
      CoCoA_ERROR(ERR::MixedRings, FnName);
    if (IsOne(c)) return;
    vector<RingElem> ans;
    ans.resize(myNumCols(), zero(myR));
    for (long j = 0; j < myNumCols(); ++j)
      myR->myMul(raw(ans[j]), raw(c), myEntries[i][j]);
    // Answer successfully computed in ans, swap it into the i-th row.
    for (long j = 0; j < myNumCols(); ++j)
      myR->mySwap(myEntries[i][j], raw(ans[j]));
  }


  void DenseMatImpl::myColMul(long j, ConstRefRingElem c)
  {
    const char* const FnName = "DenseMatImpl::myColMul";
    myCheckColIndex(j, FnName);
    if (owner(c) != myR)
      CoCoA_ERROR(ERR::MixedRings, FnName);
    if (IsOne(c)) return;
    vector<RingElem> ans;
    ans.resize(myNumRows(), zero(myR));
    for (long i = 0; i < myNumRows(); ++i)
      myR->myMul(raw(ans[i]), raw(c), myEntries[i][j]);
    // Answer successfully computed in ans, swap it into the i-th row.
    for (long i = 0; i < myNumRows(); ++i)
      myR->mySwap(myEntries[i][j], raw(ans[i]));
  }


  void DenseMatImpl::myAddRowMul(long i1, long i2, ConstRefRingElem c)
  {
    const char* const FnName = "DenseMatImpl::myAddRowMul";
    myCheckRowIndex(i1, FnName);
    myCheckRowIndex(i2, FnName);
    if (owner(c) != myR)
      CoCoA_ERROR(ERR::MixedRings, FnName);
    if (IsZero(c)) return;
    const long ncols = myNumCols();
    vector<RingElem> ans(ncols, zero(myR));
    for (long j = 0; j < ncols; ++j)
    {
      myR->myAssign(raw(ans[j]),myEntries[i1][j]);
      if (!myR->myIsZero(myEntries[i2][j]))
        myR->myIsZeroAddMul(raw(ans[j]), raw(c), myEntries[i2][j]);
    }
    // Answer successfully computed in ans, swap it into the i-th row
    for (long j = 0; j < ncols; ++j)
      myR->mySwap(raw(ans[j]), myEntries[i1][j]);
  }


  void DenseMatImpl::myAddColMul(long j1, long j2, ConstRefRingElem c)
  {
    const char* const FnName = "DenseMatImpl::myAddColMul";
    myCheckColIndex(j1, FnName);
    myCheckColIndex(j2, FnName);
    if (owner(c) != myR)
      CoCoA_ERROR(ERR::MixedRings, FnName);
    if (IsZero(c)) return;
    vector<RingElem> ans;
    ans.resize(myNumRows(), zero(myR));
    for (long i = 0; i < myNumRows(); ++i)
    {
      myR->myAssign(raw(ans[i]),myEntries[i][j1]);
      myR->myIsZeroAddMul(raw(ans[i]), raw(c), myEntries[i][j2]);
    }
    // Answer successfully computed in ans, swap it into the j-th col
    for (long i = 0; i < myNumRows(); ++i)
      myR->mySwap(raw(ans[i]), myEntries[i][j1]);
  }


  void DenseMatImpl::mySwapRows(long i1, long i2)
  {
    const char* const FnName = "DenseMatImpl::mySwapRows";
    myCheckRowIndex(i1, FnName);
    myCheckRowIndex(i2, FnName);
    if (i1 == i2) return;
    std::swap(myEntries[i1], myEntries[i2]); // trouble with iterators???
//     for (long j=0; j < myNumCols(); ++j)
//       myR->mySwap(myEntries[i1][j], myEntries[i2][j]);
  }


  void DenseMatImpl::mySwapCols(long j1, long j2)
  {
    const char* const FnName = "DenseMatImpl::mySwapCols";
    myCheckColIndex(j1, FnName);
    myCheckColIndex(j2, FnName);
    if (j1 == j2) return;
    // Must do this the "hard way".
    for (long i=0; i < myNumRows(); ++i)
      myR->mySwap(myEntries[i][j1], myEntries[i][j2]);
  }


  // NOT FULLY EXCEPTION CLEAN (do not see how without taking a performance hit, or making it even harder to read)
  void DenseMatImpl::myResize(long NumRows, long NumCols)
  {
    // First remove any excess rows
    if (NumRows < myNumRowsValue)
    {
      for (long i=NumRows; i < myNumRowsValue; ++i)
        for (long j=0; j < myNumColsValue; ++j)
          myR->myDelete(myEntries[i][j]);
      myEntries.resize(NumRows);
    }
    const long LenEntries = len(myEntries);
    // Now resize each remaining row.
    // Two cases: lengthen or shorten  (or no change)
    if (NumCols > myNumColsValue)
    {
      for (long i=0; i < LenEntries; ++i)
        for (long j=myNumColsValue; j < NumCols; ++j)
          myEntries[i].push_back(myR->myNew()); // could throw std::bad_alloc
    }
    if (NumCols < myNumColsValue)
    {
      RingElem useless=zero(myR);
      for (long i=0; i < LenEntries; ++i)
      {
        for (long j=NumCols; j < myNumColsValue; ++j)
          myR->myDelete(myEntries[i][j]);
        myEntries[i].resize(NumCols, raw(useless)); // second arg not used
      }
    }
    // Add new zero rows if number of rows is to increase.
    if (NumRows > myNumRowsValue)
    {
      for (long i=myNumRowsValue; i < NumRows; ++i)
      {
        myEntries.push_back(vector<RingElemRawPtr>());
        for (long j=0; j < NumCols; ++j)
          myEntries[i].push_back(myR->myNew()); // could throw std::bad_alloc
      }
    }
    myNumRowsValue = NumRows;
    myNumColsValue = NumCols;
  }

  void DenseMatImpl::myDet(RingElem& d) const
  {
    if (myNumRows()==2)
    { det2x2(d, ConstMatrixView(this)); return; }
    if (myNumRows()==3)      
    { det3x3(d, ConstMatrixView(this)); return; }
    if (IsField(myR))
    { DetByGauss(d, ConstMatrixView(this)); return; }
    if (IsIntegralDomain(myR))
    { d = DetByBareiss(ConstMatrixView(this)); return; }
    CoCoA_ERROR(ERR::NYI, "det for non integral domain");
  }


  long DenseMatImpl::myRank() const
  {
    vector<long> discard;
    return RankByGauss(discard, ConstMatrixView(this));
  }


  matrix NewDenseMat(const ring& R, long NumRows, long NumCols)
  {
    if (NumRows < 0)  CoCoA_ERROR(ERR::BadRowIndex, "NewDenseMat(R,n,m,M)");
    if (NumCols < 0)  CoCoA_ERROR(ERR::BadColIndex, "NewDenseMat(R,n,m,M)");
    return matrix(new DenseMatImpl(R, NumRows, NumCols));
  }

  matrix NewDenseMat(const ring& R, const std::vector< std::vector<long> >& VV)
  {
    if (!IsRectangular(VV))
      CoCoA_ERROR(ERR::BadMatrixSize, "NewDenseMat()");
    const long NumRows = len(VV);
    if (NumRows == 0) return NewDenseMat(R, 0, 0);
    const long NumCols = len(VV[0]);
    matrix ans(new DenseMatImpl(R, NumRows, NumCols));
    for (long i=0; i < NumRows; ++i)
      for (long j=0; j < NumCols; ++j)
        SetEntry(ans, i, j, VV[i][j]);
    return ans;
  }


  matrix NewDenseMat(const ring& R, const std::vector< std::vector<BigInt> >& VV)
  {
    if (!IsRectangular(VV))
      CoCoA_ERROR(ERR::BadMatrixSize, "NewDenseMat()");
    const long NumRows = len(VV);
    if (NumRows == 0) return NewDenseMat(R, 0, 0);
    const long NumCols = len(VV[0]);
    matrix ans(new DenseMatImpl(R, NumRows, NumCols));
    for (long i=0; i < NumRows; ++i)
      for (long j=0; j < NumCols; ++j)
        SetEntry(ans, i, j, VV[i][j]);
    return ans;
  }


  matrix NewDenseMat(const ring& R, const std::vector< std::vector<BigRat> >& VV)
  {
    if (!IsRectangular(VV))
      CoCoA_ERROR(ERR::BadMatrixSize, "NewDenseMat()");
    const long NumRows = len(VV);
    if (NumRows == 0) return NewDenseMat(R, 0, 0);
    const long NumCols = len(VV[0]);
    matrix ans(new DenseMatImpl(R, NumRows, NumCols));
    for (long i=0; i < NumRows; ++i)
      for (long j=0; j < NumCols; ++j)
        SetEntry(ans, i, j, VV[i][j]);
    return ans;
  }


  matrix NewDenseMat(const ring& R, const std::vector< std::vector<RingElem> >& VV)
  {
    if (!IsRectangular(VV))
      CoCoA_ERROR(ERR::BadMatrixSize, "NewDenseMat()");
    const long NumRows = len(VV);
    if (NumRows == 0) return NewDenseMat(R, 0, 0);
    const long NumCols = len(VV[0]);
    matrix ans(new DenseMatImpl(R, NumRows, NumCols));
    for (long i=0; i < NumRows; ++i)
      for (long j=0; j < NumCols; ++j)
        SetEntry(ans, i, j, VV[i][j]);
    return ans;
  }


  // I bet this makes lots of wasteful copies...
  matrix NewDenseMat(ConstMatrixView M)
  {
    const long r = NumRows(M);
    const long c = NumCols(M);
    matrix ans(new DenseMatImpl(RingOf(M), r, c));
    for (long i=0; i < r; ++i)
      for (long j=0; j < c; ++j)
        SetEntry(ans, i, j, M(i, j));
    return ans;
  }

} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/DenseMatrix.C,v 1.31 2014/08/26 12:53:52 abbott Exp $
// $Log: DenseMatrix.C,v $
// Revision 1.31  2014/08/26 12:53:52  abbott
// Summary: Improved dispatching in myDet (uses Bareiss for non fields)
// Author: JAA
//
// Revision 1.30  2014/08/16 13:42:26  abbott
// Summary: Added copy ctor for matrix, and myClone for MatrixBase
// Author: JAA
//
// Revision 1.29  2014/07/30 14:03:19  abbott
// Summary: Changed BaseRing into RingOf
// Author: JAA
//
// Revision 1.28  2014/04/11 15:44:27  abbott
// Summary: Renamed MatrixArith to MatrixOperations (in includes)
// Author: JAA
//
// Revision 1.27  2013/06/28 17:01:33  abbott
// Minor improvement to myAddMulRow (slightly faster)
//
// Revision 1.26  2013/06/27 16:54:21  abbott
// Changed DenseMatrix::mySwapRows -- it is much faster.
// There is a comment worrying about invalidating iterators -- don't see how though!
//
// Revision 1.25  2013/06/17 16:26:17  bigatti
// -- fixed really subtle bug in myResize()  (myNew actually created an
//    element which was not to be used)
//
// Revision 1.24  2012/10/24 12:11:15  abbott
// Changed return type of myEntry.
//
// Revision 1.23  2012/10/16 09:56:01  abbott
// Replaced  RefRingElem  by  RingElem&
//
// Revision 1.22  2012/10/11 14:57:20  abbott
// Replaced mem fn  myRefEntry  by equivalent new mem fn  myRawEntry;
// this way we remove the need for RefRingElem.
// Also the new name/semantics should discourage casual use.
//
// Revision 1.21  2012/01/30 12:50:49  abbott
// Added comment about bug in myDet.
//
// Revision 1.20  2011/11/09 14:03:40  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.19  2011/10/07 18:01:23  abbott
// Corrected ERR::Bad -- replaced by Err::BadMatrixSize
//
// Revision 1.18  2011/10/07 12:22:36  bigatti
// -- added include BigRat.H
//
// Revision 1.17  2011/10/05 09:23:46  abbott
// Added pseudo ctor from vector of vector of BigRat.
//
// Revision 1.16  2011/10/04 15:36:47  abbott
// Added new pseudo ctors for DenseMatrix (from vector of vector)
//
// Revision 1.15  2011/08/24 10:24:17  bigatti
// -- renamed QQ --> BigRat
// -- sorted #include
//
// Revision 1.14  2011/08/14 15:52:17  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.13  2011/03/10 15:53:03  bigatti
// -- using len instead of size()
//
// Revision 1.12  2011/03/09 09:07:00  bigatti
// -- changed row/col args into long (instead of MachineInt)
//
// Revision 1.11  2011/03/08 17:54:09  bigatti
// -- changed: args for rows and cols are now  long  instead of  MachineInt
//
// Revision 1.10  2011/03/04 16:22:05  bigatti
// -- changed: functions args of type MachineInt instead of size_t
//             members functions args of type long instead of size_t
//
// Revision 1.9  2011/02/28 14:08:49  bigatti
// -- added det3x3
// -- using apply mapping matrix (in DetByGauss)
//
// Revision 1.8  2011/01/31 14:10:30  bigatti
// -- added mySetEntry/SetEntry with BigRat entry
//
// Revision 1.7  2009/09/24 13:41:21  abbott
// Removed useless include directives.
// Removed unnecessary "std::" prefixes.
//
// Revision 1.6  2008/12/17 12:11:52  abbott
// Changed type from long to MachineInt in operations which use a machine integer
// in place of a RingElem.  The change is "superficial" but affects many files.
//
// Revision 1.5  2008/04/18 15:35:57  abbott
// (long overdue) Major revision to matrices
//
// Revision 1.4  2008/04/16 17:24:17  abbott
// Further cleaning of the new matrix code.  Updated documentation too.
//
// Revision 1.3  2008/04/08 15:26:42  abbott
// Major revision to matrix implementation: added matrix views.
// Lots of changes.
//
// Revision 1.2  2007/10/30 17:14:08  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.6  2007/03/08 18:22:30  cocoa
// Just whitespace cleaning.
//
// Revision 1.5  2007/02/10 18:44:03  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.4  2006/11/27 16:18:33  cocoa
// -- moved classes declarations from .H to .C (DenseMatrix, DiagMatrix,
//    FieldIdeal, SpecialMatrix)
//
// Revision 1.3  2006/10/16 23:18:59  cocoa
// Corrected use of std::swap and various special swap functions.
// Improved myApply memfn for homs of RingDistrMPolyInlPP.
//
// Revision 1.2  2006/10/06 14:04:15  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.2  2006/01/25 15:58:19  cocoa
// -- fixed: size condition in mul(matrix&, ConstMatrix, ConstMatrix), by Kaspar
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.2  2005/05/06 11:57:15  cocoa
// Fixed bug in ctor for DenseMatrixImpl.
//
// Revision 1.1.1.1  2005/05/03 15:47:31  cocoa
// Imported files
//
// Revision 1.9  2005/04/20 15:40:48  cocoa
// Major change: modified the standard way errors are to be signalled
// (now via a macro which records filename and line number).  Updated
// documentation in error.txt accordingly.
//
// Improved the documentation in matrix.txt (still more work to be done).
//
// Revision 1.8  2005/04/19 15:39:55  cocoa
// Matrices now use reference counts.
//
// Revision 1.7  2005/04/19 14:06:04  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.6  2005/03/29 17:26:29  cocoa
// -- changed: myDet, myRank call MatrixArith functions
//
// Revision 1.5  2005/03/11 16:44:18  cocoa
// New abstract class structure for matrices.
// New types of special matrix.
//
// Revision 1.4  2005/03/02 18:46:41  cocoa
// Added new types ConstRefMatrix, and RefMatrix following along
// the lines of ConstRefRingElem and RefRingElem.  The semantics
// should be a bit clearer now.
//
// Revision 1.3  2005/02/11 14:15:20  cocoa
// New style ring elements and references to ring elements;
// I hope I have finally got it right!
//
// Revision 1.2  2005/02/09 11:19:30  cocoa
// -- fixed: determinant over a GCDRing (through FrF)
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.6  2004/11/29 16:22:35  cocoa
// -- added function for computing adjoint and inverse for DenseMatrix
//    (so adjoint/inverse matrix is computed by OrdvArith and is no
//    longer needed by PPOrdering)
//
// Revision 1.5  2004/11/02 14:49:03  cocoa
// -- new code for matrix orderings
//
// Revision 1.4  2004/07/20 15:04:06  cocoa
// The next step in the new "ring element" conversion process:
// handling the case of creating a "const RefRingElem" object
// (since C++ refuses to do this properly itself).
//
// Revision 1.3  2004/03/20 17:46:11  cocoa
// Check in prior to departure to RWCA
//
// Revision 1.2  2004/03/09 16:15:41  cocoa
// First version of matrices.  A simple exmaple compiles and runs,
// so I'm checking in.
//
// Revision 1.1  2004/01/28 16:29:17  cocoa
// First attempt at implementing matrices.
//
//
