//   Copyright (c)  2005,2008  John Abbott

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


// Source code for classes matrix, MatrixView and MatrixBase, MatrixViewBase

#include "CoCoA/matrix.H"
#include "CoCoA/assert.H"
#include "CoCoA/error.H"
#include "CoCoA/MatrixOperations.H"
#include "CoCoA/OpenMath.H"
#include "CoCoA/PolyRing.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/utils.H"

#include <iostream>
using std::ostream;
//#include <vector>
using std::vector;

namespace CoCoA
{


  // BUG: placeholder for Laura Torrente
  bool IsZero(const std::vector<RingElem>& v)
  {
    const long n = len(v);
    for (long i=0; i < n; ++i)
      if (!IsZero(v[i])) return false;
    return true;
  }


  RingElemAlias ConstMatrixView::operator()(long i, long j) const
  {  
    mySmartPtr->myCheckRowIndex(i, "M(i,j)");
    mySmartPtr->myCheckColIndex(j, "M(i,j)");
    return mySmartPtr->myEntry(i, j);
  }



  // Check that the row index is in range (i.e. non-neg and not too large).
  // Throws ERR::BadRowIndex if not.
  void ConstMatrixViewBase::myCheckRowIndex(long i, const char* where) const
  {
    //    if (IsNegative(i) || AsUnsignedLong(i) >= myNumRows())
    if (i < 0 || i >= myNumRows())
      CoCoA_ERROR(ERR::BadRowIndex, where);
  }

  // Check that the col index is in range (i.e. non-neg and not too large).
  // Throws ERR::BadColIndex if not.
  void ConstMatrixViewBase::myCheckColIndex(long j, const char* where) const
  {
    //    if (IsNegative(j) || AsUnsignedLong(j) >= myNumCols())
    if (j < 0 || j >= myNumCols())
      CoCoA_ERROR(ERR::BadColIndex, where);
  }


  bool ConstMatrixViewBase::IamEqual(ConstMatrixView M) const
  {
    //std::clog << "ConstMatrixViewBase::IamEqual" << std::endl;
    //std::cout << myNumRows() << " and " <<  NumRows(M) << std::endl;
    //std::cout << myNumCols() << " and " <<  NumCols(M) << std::endl;
    if (myNumRows() != NumRows(M)) return false;
    if (myNumCols() != NumCols(M)) return false;
    for (long i=0; i < myNumRows(); ++i)
      for (long j=0; j < myNumCols(); ++j)
      {
        //std::cout << myEntry(i,j) << " and " <<  M(i,j) << std::endl;
        if (myEntry(i,j) != M(i,j)) return false;
      }
    
    return true;
  }
  

  bool ConstMatrixViewBase::IamSymmetric() const
  {
    CoCoA_ASSERT(myNumRows() == myNumCols());
    for (long i=0; i < myNumRows(); ++i)
      for (long j=i+1; j < myNumCols(); ++j)
        if (myEntry(i,j) != myEntry(j,i)) return false;
    return true;
  }
  

  bool ConstMatrixViewBase::IamAntiSymmetric() const
  {
    CoCoA_ASSERT(myNumRows() == myNumCols());
    for (long i=0; i < myNumRows(); ++i)
      if (!IsZero(myEntry(i,i))) return false;
    for (long i=0; i < myNumRows(); ++i)
      for (long j=i+1; j < myNumCols(); ++j)
        if (myEntry(i,j) != -myEntry(j,i)) return false;
    return true;
  }
  

  bool ConstMatrixViewBase::IamDiagonal() const
  {
    CoCoA_ASSERT(myNumRows() == myNumCols());
    for (long i=0; i < myNumRows(); ++i)
      for (long j=0; j < myNumCols(); ++j)
        if (i!=j && !IsZero(myEntry(i,j))) return false;
    return true;
  }
  

  bool ConstMatrixViewBase::myIsZeroRow(long i) const
  {
    CoCoA_ASSERT(i < myNumRows());
    for (long j=0; j < myNumCols(); ++j)
      if (!IsZero(myEntry(i,j))) return false;
    return true;
  }


  bool ConstMatrixViewBase::myIsZeroCol(long j) const
  {
    CoCoA_ASSERT(j < myNumCols());
    for (long i=0; i < myNumRows(); ++i)
      if (!IsZero(myEntry(i,j))) return false;
    return true;
  }


  void ConstMatrixViewBase::myMulByRow(vec& lhs, const vec& v) const
  {
    CoCoA_ASSERT(len(lhs) == myNumCols() && len(v) == myNumRows());
    // Put result into ans to make code exception clean.
    vector<RingElem> ans(myNumCols(), zero(myRing()));
    for (long j=0; j < myNumCols(); ++j)
      for (long i=0; i < myNumRows(); ++i)
        ans[j] += v[i]*myEntry(i, j);
    // We have successfully computed the answer, now swap it into lhs.
//???    swap(lhs, ans); // Would this invalidate iterators on lhs???
    for (long j=0; j < myNumCols(); ++j)
      swap(lhs[j], ans[j]);
  }


  void ConstMatrixViewBase::myMulByCol(vec& lhs, const vec& v) const
  {
    CoCoA_ASSERT(len(lhs) == myNumRows() && len(v) == myNumCols());
    // Put result into ans to make code exception clean.
    vector<RingElem> ans(myNumRows(), zero(myRing()));
    for (long i=0; i < myNumRows(); ++i)
      for (long j=0; j < myNumCols(); ++j)
        ans[i] += myEntry(i, j)*v[j];

    // We have successfully computed the answer, now swap it into lhs.
//???    swap(lhs, ans); // Would this invalidate iterators on lhs???
    for (long i=0; i < myNumRows(); ++i)
      swap(lhs[i], ans[i]);
  }


  void ConstMatrixViewBase::myDet(RingElem& d) const
  {
    CoCoA_ASSERT(owner(d) == myRing());
    CoCoA_ASSERT(myNumRows() == myNumCols());
    DetByGauss(d, ConstMatrixView(this));
  }


  long ConstMatrixViewBase::myRank() const
  {
    vector<long> discard;
    return RankByGauss(discard, ConstMatrixView(this));
  }


  void ConstMatrixViewBase::myOutputSelf(std::ostream& out) const
  {
    //    out << "matrix(" << myRing() << ")(" << myNumRows() << ", " << myNumCols() << ")\n[\n";
    if (IsZZ(myRing())) out << "matrix(ZZ,";
    else
    {
      if (IsQQ(myRing())) out << "matrix(QQ,";
      else
      {
        if (IsPolyRing(myRing()) && 
            myRing() == RingQQt(NumIndets(myRing())))
          out << "matrix( RingQQt(" << NumIndets(myRing()) << "),";
        else out << "matrix( /*" << myRing() << "*/";
      }
    }
    out << "\n [[";
    for (long i=0; i < myNumRows(); ++i)
    {
      if (i >0)  out << "],\n  [";
      for (long j=0; j < myNumCols(); ++j)
      {
        if (j > 0) out << ", ";
        out << myEntry(i, j);
      }
    }
    out << "]])";
  }


  void ConstMatrixViewBase::myOutputSelf(OpenMathOutput& OMOut) const
  {
    OMOut->mySendApplyStart();
    OMOut << OpenMathSymbol("linalg2", "matrix");
    OMOut << myRing() << myNumRows() << myNumCols();
    for (long i=0; i < myNumRows(); ++i)
    {
      OMOut->mySendApplyStart();
      OMOut << OpenMathSymbol("linalg2", "matrixrow"); //??? << myNumCols; ???
      for (long j=0; j < myNumCols(); ++j)
      {
        OMOut << myEntry(i, j); //??? print only raw value???
      }
      OMOut->mySendApplyEnd();
    }
    OMOut->mySendApplyEnd();
  }


  /***************************************************************************/

  matrix::matrix(const matrix& M):
      MatrixView(M->myClone())
  {}

  /***************************************************************************/

  std::ostream& operator<<(std::ostream& out, const ConstMatrixView& M)
  {
    M->myOutputSelf(out);
    return out;
  }


  void AssignZero(MatrixView& M)
  {
    M->myAssignZero();
  }


  void SetEntry(MatrixView& M, long i, long j, ConstRefRingElem r)
  {
    if (owner(r) != RingOf(M))
      CoCoA_ERROR(ERR::MixedRings, "SetEntry(M,i,j,r)");
    M->myCheckRowIndex(i, "SetEntry(M,i,j,r)");
    M->myCheckColIndex(j, "SetEntry(M,i,j,r)");
    if (!M->myIsWritable(i, j))
      CoCoA_ERROR(ERR::BadMatrixSetEntry, "SetEntry(M,i,j,r)");
    M->mySetEntry(i, j, r);
  }


  void SetEntry(MatrixView& M, long i, long j, const MachineInt& n)
  {
    M->myCheckRowIndex(i, "SetEntry(M,i,j,n)");
    M->myCheckColIndex(j, "SetEntry(M,i,j,n)");
    if (!M->myIsWritable(i, j))
      CoCoA_ERROR(ERR::BadMatrixSetEntry, "SetEntry(M,i,j,r)");
    M->mySetEntry(i, j, n);
  }


  void SetEntry(MatrixView& M, long i, long j, const BigInt& N)
  {
    M->myCheckRowIndex(i, "SetEntry(M,i,j,N)");
    M->myCheckColIndex(j, "SetEntry(M,i,j,N)");
    if (!M->myIsWritable(i, j))
      CoCoA_ERROR(ERR::BadMatrixSetEntry, "SetEntry(M,i,j,r)");
    M->mySetEntry(i, j, N);
  }


  void SetEntry(MatrixView& M, long i, long j, const BigRat& Q)
  {
    M->myCheckRowIndex(i, "SetEntry(M,i,j,Q)");
    M->myCheckColIndex(j, "SetEntry(M,i,j,Q)");
    if (!M->myIsWritable(i, j))
      CoCoA_ERROR(ERR::BadMatrixSetEntry, "SetEntry(M,i,j,r)");
    M->mySetEntry(i, j, Q);
  }


  void SwapRows(matrix& M, long i1, long i2)
  {
    M->myCheckRowIndex(i1, "SwapRows");
    M->myCheckRowIndex(i2, "SwapRows");
    return M->mySwapRows(i1,i2);
  }


  void SwapCols(matrix& M, long i1, long i2)
  {
    M->myCheckColIndex(i1, "SwapCols");
    M->myCheckColIndex(i2, "SwapCols");
    return M->mySwapCols(i1,i2);
  }


  void swap(matrix& M1, matrix& M2)
  {
    M1.mySmartPtr.mySwap(M2.mySmartPtr);
  }


  void DeleteRow(matrix& M, long i)
  {
    M->myCheckRowIndex(i, "DeleteRow");
    for (long j=i+1; j<NumRows(M); ++j)
      M->mySwapRows(j,j-1);
    M->myResize(NumRows(M)-1, NumCols(M));
  }


  void DeleteCol(matrix& M, long i)
  {
    M->myCheckColIndex(i, "DeleteCol");
    for (long j=i+1; j<NumCols(M); ++j)
      M->mySwapCols(j,j-1);
    M->myResize(NumRows(M), NumCols(M)-1);
  }

} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/matrix.C,v 1.33 2014/08/16 13:42:26 abbott Exp $
// $Log: matrix.C,v $
// Revision 1.33  2014/08/16 13:42:26  abbott
// Summary: Added copy ctor for matrix, and myClone for MatrixBase
// Author: JAA
//
// Revision 1.32  2014/07/30 14:13:05  abbott
// Summary: Changed BaseRing into RingOf; myBaseRing --> myRing
// Author: JAA
//
// Revision 1.31  2014/07/07 13:27:29  abbott
// Summary: Removed AsPolyRing
// Author: JAA
//
// Revision 1.30  2014/04/30 16:27:30  abbott
// Summary: Replaced X.size() by len(X)
// Author: JAA
//
// Revision 1.29  2014/04/17 14:08:05  bigatti
// -- moved Is... functions from matrix.H to here
//
// Revision 1.28  2014/04/11 15:44:28  abbott
// Summary: Renamed MatrixArith to MatrixOperations (in includes)
// Author: JAA
//
// Revision 1.27  2014/03/19 15:58:16  bigatti
// -- now printing the ring explicitely (when possible)
//
// Revision 1.26  2013/07/12 15:07:18  abbott
// Corrected a string (FnName passed to index checking fn).
//
// Revision 1.25  2013/06/17 12:56:57  bigatti
// -- added DeleteRow, DeleteCol
//
// Revision 1.24  2013/03/28 15:15:34  abbott
// Modified printing of matrices (removed two newlines); consequential changes to expected output of certain tests.
//
// Revision 1.23  2012/10/24 12:24:42  abbott
// Changed return type of matrix::operator().
//
// Revision 1.22  2012/10/16 09:50:40  abbott
// Changed  RefRingElem  into  RingElem&.
// Realigned some comments.
//
// Revision 1.21  2012/04/04 13:53:54  bigatti
// -- added SwapRows, SwapCols
//
// Revision 1.20  2012/03/27 11:53:29  bigatti
// -- changed printing style (and added ring)
//
// Revision 1.19  2011/11/09 17:08:53  bigatti
// -- changed output for matrix
//
// Revision 1.18  2011/11/09 14:29:37  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.17  2011/08/24 10:32:04  bigatti
// -- renamed QQ --> BigRat
// -- sorted #include
//
// Revision 1.16  2011/08/14 15:52:16  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.15  2011/03/10 17:57:57  bigatti
// -- fixed some CoCoA_ASSERT
//
// Revision 1.14  2011/03/09 09:11:20  bigatti
// -- removed call to AsSignedLong
//
// Revision 1.13  2011/03/08 17:27:01  bigatti
// -- changed: args for rows and cols are now  long  instead of  MachineInt
//
// Revision 1.12  2011/03/04 16:26:14  bigatti
// -- changed: functions args of type MachineInt instead of size_t
//             members functions args of type long instead of size_t
//
// Revision 1.11  2011/03/03 13:50:21  abbott
// Replaced several occurrences of std::size_t by long; there's still more
// work to do though!
//
// Revision 1.10  2011/02/22 14:23:59  bigatti
// -- fixed IamAntiSymmetric
//
// Revision 1.9  2011/02/22 13:15:10  bigatti
// -- added IsSymmetric, IsAntiSymmetric, IsDiagonal, and tests
//
// Revision 1.8  2011/02/15 10:04:05  bigatti
// -- added myIsEqual (default impl)
// -- added check for myIsWritable in SetEntry (4 defn)
//
// Revision 1.7  2011/01/31 14:10:30  bigatti
// -- added mySetEntry/SetEntry with BigRat entry
//
// Revision 1.6  2009/12/03 17:40:36  abbott
// Added include directives for ZZ.H (as a consequence of removing
// the directive from ring.H).
//
// Revision 1.5  2008/12/17 12:11:52  abbott
// Changed type from long to MachineInt in operations which use a machine integer
// in place of a RingElem.  The change is "superficial" but affects many files.
//
// Revision 1.4  2008/04/16 17:24:17  abbott
// Further cleaning of the new matrix code.  Updated documentation too.
//
// Revision 1.3  2008/04/08 15:26:42  abbott
// Major revision to matrix implementation: added matrix views.
// Lots of changes.
//
// Revision 1.2  2007/10/30 17:14:06  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.4  2007/03/08 18:22:28  cocoa
// Just whitespace cleaning.
//
// Revision 1.3  2006/11/02 13:25:43  cocoa
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
// Revision 1.3  2006/05/12 16:10:58  cocoa
// Added OpenMathFwd.H, and tidied OpenMath.H.
// Many consequential but trivial changes.
//
// Revision 1.2  2006/04/27 13:45:30  cocoa
// Changed name of NewIdentityRingHom to NewIdentityHom.
// Changed name of member functions which print out their own object
// into myOutputSelf (to distinguish from "transitive" myOutput fns).
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.1.1.1  2005/05/03 15:47:31  cocoa
// Imported files
//
// Revision 1.11  2005/04/27 16:14:56  cocoa
// Cleaned up example programs -- added "free use" permit.
// Changed a couple of ErrorInfo object names, and added
// ERR::NotTrueGCDDomain.
//
// Revision 1.10  2005/04/20 15:40:47  cocoa
// Major change: modified the standard way errors are to be signalled
// (now via a macro which records filename and line number).  Updated
// documentation in error.txt accordingly.
//
// Improved the documentation in matrix.txt (still more work to be done).
//
// Revision 1.9  2005/04/19 15:39:55  cocoa
// Matrices now use reference counts.
//
// Revision 1.8  2005/04/19 14:06:03  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.7  2005/03/31 16:59:42  cocoa
// Made special matrix ctors private, so a user has to pass via the
// pseudo-ctors (which do all the arg sanity checking).
//
// Revision 1.6  2005/03/30 17:15:14  cocoa
// Cleaned the SpecialMatrix code; a simple test compiles and
// runs fine.
//
// Revision 1.5  2005/03/29 17:19:59  cocoa
// -- added: rank
//
// Revision 1.4  2005/03/11 16:44:18  cocoa
// New abstract class structure for matrices.
// New types of special matrix.
//
// Revision 1.3  2005/03/02 18:46:41  cocoa
// Added new types ConstRefMatrix, and RefMatrix following along
// the lines of ConstRefRingElem and RefRingElem.  The semantics
// should be a bit clearer now.
//
// Revision 1.2  2005/02/11 14:15:20  cocoa
// New style ring elements and references to ring elements;
// I hope I have finally got it right!
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.9  2004/11/29 16:19:02  cocoa
// -- changed syntax for function det
//
// Revision 1.8  2004/11/12 15:49:29  cocoa
// Tidying prior to 0.90 release.
// (a) documentation improved (or marked as poor)
// (b) sundry minor improvements to the code
//
// Revision 1.7  2004/11/11 13:55:10  cocoa
// -- minor changes for doxygen
//
// Revision 1.6  2004/11/02 14:49:03  cocoa
// -- new code for matrix orderings
//
// Revision 1.5  2004/07/20 15:04:06  cocoa
// The next step in the new "ring element" conversion process:
// handling the case of creating a "const RefRingElem" object
// (since C++ refuses to do this properly itself).
//
// Revision 1.4  2004/05/24 15:52:13  cocoa
// Major update:
//   new error mechanism
//   many fixes
//   RingHoms almost work now
//   RingFloat much improved
//
// Revision 1.3  2004/03/20 17:46:10  cocoa
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
