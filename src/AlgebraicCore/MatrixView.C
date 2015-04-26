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


#include "CoCoA/MatrixView.H"
#include "CoCoA/MachineInt.H"
#include "CoCoA/assert.H"
#include "CoCoA/error.H"
#include "CoCoA/matrix.H"
#include "CoCoA/ring.H"
#include "CoCoA/utils.H" // for len
#include "CoCoA/VectorOperations.H" // for HasUniqueOwner


//#include <vector>
using std::vector;
#include <algorithm>
using std::min;
//#include <limits> // included in utils.H
using std::numeric_limits;
//#include <iostream> // for debugging only

namespace CoCoA
{

  /***************************************************************************/
  bool IsZeroMatImpl(ConstMatrixView M);
  bool IsIdentityMatImpl(ConstMatrixView M);
  bool IsConstDiagMatImpl(ConstMatrixView M);
  bool IsDiagMatImpl(ConstMatrixView M);
  /***************************************************************************/
  // ZeroMatImpl -- zero matrices
  // Derived from MatrixViewBase rather than ConstMatrixViewBase so that BlockMat2x2 works well.

  class ZeroMatImpl: public MatrixViewBase
  {
  private:
    friend MatrixView ZeroMat(const ring&, long, long); // pseudo-ctor
    ZeroMatImpl(const ring& R, long nrows, long ncols);
    // default dtor works fine
  private: // disable default copy ctor and assignment
    ZeroMatImpl(const ZeroMatImpl&);            // NEVER DEFINED -- copy ctor disabled
    ZeroMatImpl& operator=(const ZeroMatImpl&); // NEVER DEFINED -- copy ctor disabled

  public:
    typedef std::vector<RingElem> vec;
    const ring& myRing() const;
    virtual long myNumRows() const;
    virtual long myNumCols() const;
    virtual RingElemAlias myEntry(long i, long j) const;
    virtual void myMulByRow(vec& lhs, const vec& v) const;
    virtual void myMulByCol(vec& lhs, const vec& v) const;
    virtual bool IamEqual(ConstMatrixView M) const;
    virtual bool IamSymmetric() const {return true;}
    virtual bool IamAntiSymmetric() const {return true;}
    virtual bool IamDiagonal() const {return true;}
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
    long myNumRowsValue;
    long myNumColsValue;
  };


  /******************************************************************/
  // Functions for ZeroMatImpl

  bool IsZeroMatImpl(ConstMatrixView M)  
  {return dynamic_cast<ZeroMatImpl*>(M.mySmartPtr.myRawPtr()) != 0;}


  ZeroMatImpl::ZeroMatImpl(const ring& R, long NumRows, long NumCols):
      MatrixViewBase(),
      myR(R),
      myNumRowsValue(NumRows),
      myNumColsValue(NumCols)
  {}


  const ring& ZeroMatImpl::myRing() const
  {
    return myR;
  }


  long ZeroMatImpl::myNumRows() const
  {
    return myNumRowsValue;
  }


  long ZeroMatImpl::myNumCols() const
  {
    return myNumColsValue;
  }


  RingElemAlias ZeroMatImpl::myEntry(long i, long j) const
  {
    (void)(i); (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    return zero(myR);
  }


  // This is exception safe only if "lhs[i]=0;" cannot throw.
  void ZeroMatImpl::myMulByRow(vec& lhs, const vec& v) const
  {
    (void)(v); // to avoid compiler warning about unused parameter
// BUG: should check all elems of v belong to the right ring!
    CoCoA_ASSERT(len(lhs) == myNumCols() && len(v) == myNumRows());
    for (long i=0; i < myNumCols(); ++i)
      lhs[i] = 0;
  }


  // This is exception safe only if "lhs[i]=0;" cannot throw.
  void ZeroMatImpl::myMulByCol(vec& lhs, const vec& v) const
  {
    (void)(v); // to avoid compiler warning about unused parameter
// BUG: should check all elems of v belong to the right ring!
    CoCoA_ASSERT(len(lhs) == myNumRows() && len(v) == myNumCols());
    for (long i=0; i < myNumRows(); ++i)
      lhs[i] = 0;
  }


  bool ZeroMatImpl::IamEqual(ConstMatrixView M) const
  {
    //    if (myRing() != RingOf(M)) return false;
    if (myRing() != RingOf(M))
      CoCoA_ERROR(ERR::MixedRings, "IamEqual");
    if (myNumRows() != NumRows(M)) return false;
    if (myNumCols() != NumCols(M)) return false;
    if (myNumCols() == 0 || myNumRows() == 0) return true;
    if (IsZeroMatImpl(M))     return true;
    if (IsIdentityMatImpl(M)) return false;
    if (IsDiagMatImpl(M) || IsConstDiagMatImpl(M))
    {
      //      std::cout << "ZeroMatImpl::IamEqual - diag" << std::endl;
      for (long i=0; i < myNumRows(); ++i) if (!IsZero(M(i,i))) return false;
      return true;
    }
    //    std::cout << "ZeroMatImpl::IamEqual" << std::endl;
    for (long i=0; i < myNumRows(); ++i)
      for (long j=0; j < myNumCols(); ++j)
        if (!IsZero(M(i,j))) return false;
    //    return ConstMatrixViewBase::IamEqual(M);
    return true;
  }


  bool ZeroMatImpl::myIsZeroRow(long i) const
  {
    (void)(i); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(0 <= 0 && i < myNumRows());
    return true;
  }


  bool ZeroMatImpl::myIsZeroCol(long j) const
  {
    (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(0 <= j && j < myNumCols());
    return true;
  }


  void ZeroMatImpl::myDet(RingElem& d) const
  {
    CoCoA_ASSERT(owner(d) == myRing());
    CoCoA_ASSERT(myNumRows() == myNumCols());
    if (myNumRows() == 0)
      d = 1;
    else
      d = 0;
  }



  long ZeroMatImpl::myRank() const
  {
    return 0;
  }


  // Do nothing - matrix is already full of zeros!
  void ZeroMatImpl::myAssignZero()
  {}


  bool ZeroMatImpl::myIsWritable(long i, long j)
  {
    (void)(i); (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    return false;
  }


  RingElemRawPtr ZeroMatImpl::myRawEntry(long /*i*/, long /*j*/)
  {
    CoCoA_ERROR(ERR::ConstMatEntry, "ZeroMat->myRawEntry(i,j)");
    return RingElemRawPtr(0); // never reached, just to keep compiler quiet
  }


  void ZeroMatImpl::mySetEntry(long /*i*/, long /*j*/, ConstRefRingElem /*r*/)
  {
    CoCoA_ERROR(ERR::ConstMatEntry, "SetEntry(ZeroMat,i,j,RingElem)");
  }

  void ZeroMatImpl::mySetEntry(long /*i*/, long /*j*/, const MachineInt& /*n*/)
  {
    CoCoA_ERROR(ERR::ConstMatEntry, "SetEntry(ZeroMat,i,j,n)");
  }

  void ZeroMatImpl::mySetEntry(long /*i*/, long /*j*/, const BigInt& /*N*/)
  {
    CoCoA_ERROR(ERR::ConstMatEntry, "SetEntry(ZeroMat,i,j,BigInt)");
  }


  void ZeroMatImpl::mySetEntry(long /*i*/, long /*j*/, const BigRat& /*Q*/)
  {
    CoCoA_ERROR(ERR::ConstMatEntry, "SetEntry(ZeroMat,i,j,BigRat)");
  }


  MatrixView ZeroMat(const ring& R, long NumRows, long NumCols)
  {
    if (NumRows < 0)  CoCoA_ERROR(ERR::BadRowIndex, "ZeroMat(R,n,m)");
    if (NumCols < 0)  CoCoA_ERROR(ERR::BadColIndex, "ZeroMat(R,n,m)");
    return MatrixView(new ZeroMatImpl(R, NumRows, NumCols));
  }


  /***************************************************************************/
  // IdentityMatImpl -- identity matrices (necessarily square)
  // Derived from MatrixViewBase rather than ConstMatrixViewBase so that BlockMat2x2 works well.

  class IdentityMatImpl: public MatrixViewBase
  {
  private:
    friend MatrixView IdentityMat(const ring& R, long dim); // pseudo-ctor
    IdentityMatImpl(const ring& R, long dim);
    // default dtor is fine
  private: // disable default copy ctor and assignment
    IdentityMatImpl(const IdentityMatImpl&);            // NEVER DEFINED -- copy ctor disabled
    IdentityMatImpl& operator=(const IdentityMatImpl&); // NEVER DEFINED -- copy ctor disabled

  public:
    typedef std::vector<RingElem> vec;
    const ring& myRing() const;
    virtual long myNumRows() const;
    virtual long myNumCols() const;
    virtual RingElemAlias myEntry(long i, long j) const;
    virtual void myMulByRow(vec& lhs, const vec& v) const;
    virtual void myMulByCol(vec& lhs, const vec& v) const;
    virtual bool IamEqual(ConstMatrixView M) const;
    virtual bool IamSymmetric() const {return true;}
    virtual bool IamAntiSymmetric() const {return myNumRows()==0;}
    virtual bool IamDiagonal() const {return true;}
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


  /***************************************************************************/
  // Functions for IdentityMatImpl

  bool IsIdentityMatImpl(ConstMatrixView M)
  {return dynamic_cast<const IdentityMatImpl*>(M.mySmartPtr.myRawPtr()) != 0;}


  IdentityMatImpl::IdentityMatImpl(const ring& R, long dim):
      MatrixViewBase(),
      myR(R),
      myDim(dim)
  {}


  const ring& IdentityMatImpl::myRing() const
  {
    return myR;
  }


  long IdentityMatImpl::myNumRows() const
  {
    return myDim;
  }


  long IdentityMatImpl::myNumCols() const
  {
    return myDim;
  }


  RingElemAlias IdentityMatImpl::myEntry(long i, long j) const
  {
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    if (i == j) return one(myR);
    return zero(myR);
  }


  //BUG: should use FreeModule elems!
  void IdentityMatImpl::myMulByRow(vec& lhs, const vec& v) const
  {
    CoCoA_ASSERT(len(lhs) == myNumCols() && len(v) == myNumRows());
    CoCoA_ASSERT(myNumRows() == 0 || myRing() == owner(v[0]));
    CoCoA_ASSERT(myNumRows() == 0 || myRing() == owner(lhs[0]));
    lhs = v;
  }


  void IdentityMatImpl::myMulByCol(vec& lhs, const vec& v) const
  {
    CoCoA_ASSERT(len(lhs) == myNumRows() && len(v) == myNumCols());
//BUG???    CoCoA_ASSERT(myRing() == RingOf(v));
    lhs = v;
  }


  bool IdentityMatImpl::IamEqual(ConstMatrixView M) const
  {
    if (myRing() != RingOf(M)) return false;
    if (myNumRows() != NumRows(M)) return false;
    if (myNumCols() != NumCols(M)) return false;
    if (IsZeroMatImpl(M)) return NumCols(M) == 0;
    if (IsIdentityMatImpl(M)) return true;
    if (IsDiagMatImpl(M) || IsConstDiagMatImpl(M))
    {
      //      std::cout << "IdentityMatImpl::IamEqual - diag" << std::endl;
      for (long i=0; i < myNumRows(); ++i) if (!IsOne(M(i,i))) return false;
      return true;
    }
    return ConstMatrixViewBase::IamEqual(M);
  }


  bool IdentityMatImpl::myIsZeroRow(long i) const
  {
    (void)(i); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(0 <= i && i < myNumRows());
    return false;
  }


  bool IdentityMatImpl::myIsZeroCol(long j) const
  {
    (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(0 <= j && j < myNumCols());
    return false;
  }


  void IdentityMatImpl::myDet(RingElem& d) const
  {
    CoCoA_ASSERT(owner(d) == myRing());
    d = 1; // NB correct even for 0x0 matrix.
  }


  long IdentityMatImpl::myRank() const
  {
    return myNumRows();
  }


  // Always gives error!
  void IdentityMatImpl::myAssignZero()
  {
//???    if (myDim == 0) return;
    CoCoA_ERROR(ERR::ConstMatEntry, "IdentityMat->myAssignZero()");
  }


  bool IdentityMatImpl::myIsWritable(long i, long j)
  {
    (void)(i); (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    return false;
  }


  RingElemRawPtr IdentityMatImpl::myRawEntry(long /*i*/, long /*j*/)
  {
    CoCoA_ERROR(ERR::ConstMatEntry, "IdentityMat->myRawEntry(i,j)");
    return RingElemRawPtr(0); // never reached, just to keep compiler quiet
  }


  void IdentityMatImpl::mySetEntry(long /*i*/, long /*j*/, ConstRefRingElem /*r*/)
  {
    CoCoA_ERROR(ERR::ConstMatEntry, "SetEntry(IdentityMat,i,j,RingElem)");
  }

  void IdentityMatImpl::mySetEntry(long /*i*/, long /*j*/, const MachineInt& /*n*/)
  {
    CoCoA_ERROR(ERR::ConstMatEntry, "SetEntry(IdentityMat,i,j,n)");
  }

  void IdentityMatImpl::mySetEntry(long /*i*/, long /*j*/, const BigInt& /*N*/)
  {
    CoCoA_ERROR(ERR::ConstMatEntry, "SetEntry(IdentityMat,i,j,BigInt)");
  }


  void IdentityMatImpl::mySetEntry(long /*i*/, long /*j*/, const BigRat& /*Q*/)
  {
    CoCoA_ERROR(ERR::ConstMatEntry, "SetEntry(IdentityMat,i,j,BigRat)");
  }


  MatrixView IdentityMat(const ring& R, long dim)
  {
    if (dim < 0)  CoCoA_ERROR(ERR::BadRowIndex, "IdentityMat(R,dim)");
    return MatrixView(new IdentityMatImpl(R, dim));
  }



  /***************************************************************************/
  // ConstTrasposeMatImpl -- a transposed view of a const matrix

  class ConstTrasposeMatImpl: public ConstMatrixViewBase
  {
  private:
    friend ConstMatrixView transpose(ConstMatrixView M); // pseudo-ctor
    ConstTrasposeMatImpl(ConstMatrixView M);
    // default dtor works fine
  private: // disable default copy ctor and assignment
    ConstTrasposeMatImpl(const ConstTrasposeMatImpl&);            // NEVER DEFINED -- copy ctor disabled
    ConstTrasposeMatImpl& operator=(const ConstTrasposeMatImpl&); // NEVER DEFINED -- copy ctor disabled

  public:
    typedef std::vector<RingElem> vec;
    const ring& myRing() const;
    virtual long myNumRows() const;
    virtual long myNumCols() const;
    virtual RingElemAlias myEntry(long i, long j) const;
    virtual void myMulByRow(vec& lhs, const vec& v) const;
    virtual void myMulByCol(vec& lhs, const vec& v) const;
    virtual bool myIsZeroRow(long i) const;
    virtual bool myIsZeroCol(long j) const;
    virtual void myDet(RingElem& d) const;
    virtual long myRank() const;

  private: // data members
    ConstMatrixView myM;
  };


  /***************************************************************************/
  // Functions for ConstTrasposeMatImpl -- all are simple one-liners

  ConstTrasposeMatImpl::ConstTrasposeMatImpl(ConstMatrixView M):
      myM(M)
  {}


  const ring& ConstTrasposeMatImpl::myRing() const
  {
    return RingOf(myM);
  }


  long ConstTrasposeMatImpl::myNumRows() const
  {
    return NumCols(myM);
  }


  long ConstTrasposeMatImpl::myNumCols() const
  {
    return NumRows(myM);
  }


  RingElemAlias ConstTrasposeMatImpl::myEntry(long i, long j) const
  {
    return myM(j, i);
  }


  void ConstTrasposeMatImpl::myMulByRow(vec& lhs, const vec& v) const
  {
    CoCoA_ASSERT(len(lhs) == myNumCols() && len(v) == myNumRows());
    myM->myMulByCol(lhs, v);
  }


  void ConstTrasposeMatImpl::myMulByCol(vec& lhs, const vec& v) const
  {
    CoCoA_ASSERT(len(lhs) == myNumRows() && len(v) == myNumCols());
    myM->myMulByRow(lhs, v);
  }



  bool ConstTrasposeMatImpl::myIsZeroRow(long i) const
  {
    CoCoA_ASSERT(i >= 0 && i < myNumRows());
    return myM->myIsZeroCol(i);
  }


  bool ConstTrasposeMatImpl::myIsZeroCol(long j) const
  {
    CoCoA_ASSERT(j >= 0 && j < myNumCols());
    return myM->myIsZeroRow(j);
  }


  void ConstTrasposeMatImpl::myDet(RingElem& d) const
  {
    CoCoA_ASSERT(owner(d) == myRing());
    CoCoA_ASSERT(myNumRows() == myNumCols());
    myM->myDet(d);
  }


  long ConstTrasposeMatImpl::myRank() const
  {
    return myM->myRank();
  }


  ConstMatrixView transpose(ConstMatrixView M)
  {
    return ConstMatrixView(new ConstTrasposeMatImpl(M));
  }


  /***************************************************************************/
  // TrasposeMatImpl -- a transposed view of a (non-const) matrix

  class TrasposeMatImpl: public MatrixViewBase
  {
  private:
    friend MatrixView transpose(MatrixView M); // pseudo-ctor
    TrasposeMatImpl(MatrixView M);
    // default dtor works fine
  private: // disable default copy ctor and assignment
    TrasposeMatImpl(const TrasposeMatImpl&);            // NEVER DEFINED -- copy ctor disabled
    TrasposeMatImpl& operator=(const TrasposeMatImpl&); // NEVER DEFINED -- copy ctor disabled

  public:
    typedef std::vector<RingElem> vec;
    const ring& myRing() const;
    virtual long myNumRows() const;
    virtual long myNumCols() const;
    virtual RingElemAlias myEntry(long i, long j) const;
    virtual void myMulByRow(vec& lhs, const vec& v) const;
    virtual void myMulByCol(vec& lhs, const vec& v) const;
    virtual bool myIsZeroRow(long i) const;
    virtual bool myIsZeroCol(long j) const;
    virtual void myDet(RingElem& d) const;
    virtual long myRank() const;

    // non-const member fns
    virtual void myAssignZero();
    virtual bool myIsWritable(long i, long j);
    virtual RingElemRawPtr myRawEntry(long i, long j);
    virtual void mySetEntry(long i, long j, ConstRefRingElem r);
    virtual void mySetEntry(long i, long j, const MachineInt& n);
    virtual void mySetEntry(long i, long j, const BigInt& N);
    virtual void mySetEntry(long i, long j, const BigRat& Q);

  private: // data members
    MatrixView myM;
  };


  /***************************************************************************/
  // Functions for TrasposeMatImpl -- all are simple one-liners

  TrasposeMatImpl::TrasposeMatImpl(MatrixView M):
      myM(M)
  {}


  const ring& TrasposeMatImpl::myRing() const
  {
    return RingOf(myM);
  }


  long TrasposeMatImpl::myNumRows() const
  {
    return NumCols(myM);
  }


  long TrasposeMatImpl::myNumCols() const
  {
    return NumRows(myM);
  }


  RingElemAlias TrasposeMatImpl::myEntry(long i, long j) const
  {
    return myM(j, i);
  }


  void TrasposeMatImpl::myMulByRow(vec& lhs, const vec& v) const
  {
    CoCoA_ASSERT(len(lhs) == myNumCols() && len(v) == myNumRows());
    myM->myMulByCol(lhs, v);
  }


  void TrasposeMatImpl::myMulByCol(vec& lhs, const vec& v) const
  {
    CoCoA_ASSERT(len(lhs) == myNumRows() && len(v) == myNumCols());
    myM->myMulByRow(lhs, v);
  }



  bool TrasposeMatImpl::myIsZeroRow(long i) const
  {
    CoCoA_ASSERT(i >= 0 && i < myNumRows());
    return myM->myIsZeroCol(i);
  }


  bool TrasposeMatImpl::myIsZeroCol(long j) const
  {
    CoCoA_ASSERT(j >= 0 && j < myNumCols());
    return myM->myIsZeroRow(j);
  }


  void TrasposeMatImpl::myDet(RingElem& d) const
  {
    CoCoA_ASSERT(owner(d) == myRing());
    CoCoA_ASSERT(myNumRows() == myNumCols());
    myM->myDet(d);
  }


  long TrasposeMatImpl::myRank() const
  {
    return myM->myRank();
  }


  void TrasposeMatImpl::myAssignZero()
  {
    myM->myAssignZero();
  }


  bool TrasposeMatImpl::myIsWritable(long i, long j)
  {
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    return myM->myIsWritable(j,i);
  }


  RingElemRawPtr TrasposeMatImpl::myRawEntry(long i, long j)
  {
    CoCoA_ASSERT(myIsWritable(i,j));
    return myM->myRawEntry(j,i);
  }


  void TrasposeMatImpl::mySetEntry(long i, long j, ConstRefRingElem r)
  {
    CoCoA_ASSERT(myIsWritable(i,j));
    myM->mySetEntry(j,i,r);
  }


  void TrasposeMatImpl::mySetEntry(long i, long j, const MachineInt& n)
  {
    CoCoA_ASSERT(myIsWritable(i,j));
    myM->mySetEntry(j,i,n);
  }


  void TrasposeMatImpl::mySetEntry(long i, long j, const BigInt& N)
  {
    CoCoA_ASSERT(myIsWritable(i,j));
    myM->mySetEntry(j,i,N);
  }


  void TrasposeMatImpl::mySetEntry(long i, long j, const BigRat& Q)
  {
    CoCoA_ASSERT(myIsWritable(i,j));
    myM->mySetEntry(j,i,Q);
  }


  MatrixView transpose(MatrixView M)
  {
    return MatrixView(new TrasposeMatImpl(M));
  }


  /***************************************************************************/
  // SubmatImpl -- see only certain rows/cols of a matrix

  class ConstSubmatImpl: public ConstMatrixViewBase
  {
  private:
    friend ConstMatrixView submat(ConstMatrixView M, const std::vector<long>& RowIndices, const std::vector<long>& ColIndices); // pseudo-ctor
    ConstSubmatImpl(ConstMatrixView M, const std::vector<long>& rows, const std::vector<long>& cols);
    // default dtor works fine
  private: // disable default copy ctor and assignment
    ConstSubmatImpl(const ConstSubmatImpl&);            // NEVER DEFINED -- copy ctor disabled
    ConstSubmatImpl& operator=(const ConstSubmatImpl&); // NEVER DEFINED -- assignment disabled

  public:
    typedef std::vector<RingElem> vec;
    const ring& myRing() const;
    virtual long myNumRows() const;
    virtual long myNumCols() const;
    virtual RingElemAlias myEntry(long i, long j) const;
// Use DEFAULT definitions of the following six functions.
//     virtual void myMulByRow(vec& lhs, const vec& v) const;
//     virtual void myMulByCol(vec& lhs, const vec& v) const;
//     virtual bool myIsZeroRow(long i) const;
//     virtual bool myIsZeroCol(long j) const;
//     virtual void myDet(RingElem& d) const;
//     virtual long myRank() const;

  private: // data members
    ConstMatrixView myM;
    std::vector<long> myRowTable;
    std::vector<long> myColTable;
  };


  /***************************************************************************/
  // Functions for ConstSubmatImpl

  ConstSubmatImpl::ConstSubmatImpl(ConstMatrixView M, const std::vector<long>& rows, const std::vector<long>& cols):
      myM(M),
      myRowTable(rows),
      myColTable(cols)
  {
    // Args already checked by pseudo-ctor submat, the only
    // function which calls this ctor
  }


  const ring& ConstSubmatImpl::myRing() const
  {
    return RingOf(myM);
  }


  long ConstSubmatImpl::myNumRows() const
  {
    return len(myRowTable);
  }


  long ConstSubmatImpl::myNumCols() const
  {
    return len(myColTable);
  }


  RingElemAlias ConstSubmatImpl::myEntry(long i, long j) const
  {
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    return myM(myRowTable[i], myColTable[j]);
  }


  ConstMatrixView submat(ConstMatrixView M, const std::vector<long>& rows, const vector<long>& cols)
  {
    const long LenRows = len(rows);
    const long LenCols = len(cols);
    for (long i=0; i < LenRows; ++i)
      M->myCheckRowIndex(rows[i], "submat(M,rows,cols)");
    for (long j=0; j < LenCols; ++j)
      M->myCheckColIndex(cols[j], "submat(M,rows,cols)");

    // Check all row indices are distinct.
    vector<bool> RowUsed(NumRows(M));
    for (long i=0; i < LenRows; ++i)
    {
      if (RowUsed[rows[i]])
        CoCoA_ERROR(ERR::BadRowIndex, "submat(M,rows,cols):  repeated row indices");
      RowUsed[rows[i]] = true;
    }

    // Check all col indices are distinct.
    vector<bool> ColUsed(NumCols(M));
    for (long i=0; i < LenCols; ++i)
    {
      if (ColUsed[cols[i]])
        CoCoA_ERROR(ERR::BadColIndex, "submat(M,rows,cols):  repeated col indices");
      ColUsed[cols[i]] = true;
    }

    return ConstMatrixView(new ConstSubmatImpl(M, rows, cols));
  }


  /***************************************************************************/
  // SubmatImpl -- see only certain rows/cols of a matrix

  class SubmatImpl: public MatrixViewBase
  {
  private:
    friend MatrixView submat(MatrixView M, const std::vector<long>& RowIndices, const std::vector<long>& ColIndices); // pseudo-ctor
    SubmatImpl(MatrixView M, const std::vector<long>& rows, const std::vector<long>& cols);
    // default dtor works fine
  private: // disable default copy ctor and assignment
    SubmatImpl(const SubmatImpl&);            // NEVER DEFINED -- copy ctor disabled
    SubmatImpl& operator=(const SubmatImpl&); // NEVER DEFINED -- assignment disabled

  public:
    typedef std::vector<RingElem> vec;
    const ring& myRing() const;
    virtual long myNumRows() const;
    virtual long myNumCols() const;
    virtual RingElemAlias myEntry(long i, long j) const;
// Use DEFAULT definitions of the following six functions.
//     virtual void myMulByRow(vec& lhs, const vec& v) const;
//     virtual void myMulByCol(vec& lhs, const vec& v) const;
//     virtual bool myIsZeroRow(long i) const;
//     virtual bool myIsZeroCol(long j) const;
//     virtual void myDet(RingElem& d) const;
//     virtual long myRank() const;

    // non-const member fns
    virtual void myAssignZero();
    virtual bool myIsWritable(long i, long j);
    virtual RingElemRawPtr myRawEntry(long i, long j);
    virtual void mySetEntry(long i, long j, ConstRefRingElem r);
    virtual void mySetEntry(long i, long j, const MachineInt& n);
    virtual void mySetEntry(long i, long j, const BigInt& N);
    virtual void mySetEntry(long i, long j, const BigRat& Q);

  private: // data members
    MatrixView myM;
    std::vector<long> myRowTable;
    std::vector<long> myColTable;
  };


  /***************************************************************************/
  // Functions for SubmatImpl


  SubmatImpl::SubmatImpl(MatrixView M, const std::vector<long>& rows, const std::vector<long>& cols):
      myM(M),
      myRowTable(rows),
      myColTable(cols)
  {
    // Args already checked by pseudo-ctor submat, the only
    // function which calls this ctor
  }


  const ring& SubmatImpl::myRing() const
  {
    return RingOf(myM);
  }


  long SubmatImpl::myNumRows() const
  {
    return len(myRowTable);
  }


  long SubmatImpl::myNumCols() const
  {
    return len(myColTable);
  }


  RingElemAlias SubmatImpl::myEntry(long i, long j) const
  {
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    return myM(myRowTable[i], myColTable[j]);
  }


  void SubmatImpl::myAssignZero()
  {
    for (long i=0; i < myNumRows(); ++i)
      for (long j=0; j < myNumCols(); ++j)
        mySetEntry(i,j,0);
  }


  bool SubmatImpl::myIsWritable(long i, long j)
  {
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    return myM->myIsWritable(myRowTable[i], myColTable[j]);
  }


  RingElemRawPtr SubmatImpl::myRawEntry(long i, long j)
  {
    CoCoA_ASSERT(myIsWritable(i, j));
    return myM->myRawEntry(myRowTable[i], myColTable[j]);
  }


  void SubmatImpl::mySetEntry(long i, long j, ConstRefRingElem r)
  {
    CoCoA_ASSERT(myIsWritable(i, j));
    myM->mySetEntry(myRowTable[i], myColTable[j], r);
  }

  void SubmatImpl::mySetEntry(long i, long j, const MachineInt& n)
  {
    CoCoA_ASSERT(myIsWritable(i, j));
    myM->mySetEntry(myRowTable[i], myColTable[j], n);
  }

  void SubmatImpl::mySetEntry(long i, long j, const BigInt& N)
  {
    CoCoA_ASSERT(myIsWritable(i, j));
    myM->mySetEntry(myRowTable[i], myColTable[j], N);
  }


  void SubmatImpl::mySetEntry(long i, long j, const BigRat& Q)
  {
    CoCoA_ASSERT(myIsWritable(i, j));
    myM->mySetEntry(myRowTable[i], myColTable[j], Q);
  }


  MatrixView submat(MatrixView M, const std::vector<long>& rows, const vector<long>& cols)
  {
    const long LenRows = len(rows);
    const long LenCols = len(cols);
    for (long i=0; i < LenRows; ++i)
      M->myCheckRowIndex(rows[i], "submat(M,rows,cols)");
    for (long j=0; j < LenCols; ++j)
      M->myCheckColIndex(cols[j], "submat(M,rows,cols)");

    // Check all row indices are distinct.
    for (long i=0; i < LenRows; ++i)
      for (long j=i+1; j < LenRows; ++j)
        if (rows[i] == rows[j])
          CoCoA_ERROR(ERR::BadRowIndex, "submat(M,rows,cols):  repeated row indices");

    // Check all col indices are distinct.
    for (long i=0; i < LenCols; ++i)
      for (long j=i+1; j < LenCols; ++j)
        if (cols[i] == cols[j])
          CoCoA_ERROR(ERR::BadColIndex, "submat(M,rows,cols):  repeated col indices");

    return MatrixView(new SubmatImpl(M, rows, cols));
  }


  /****************************************************************************/

  class ConstColMatImpl: public ConstMatrixViewBase
  {
  private:
    friend ConstMatrixView ColMat(const std::vector<RingElem>& v); // pseudo-ctor
    ConstColMatImpl(const std::vector<RingElem>& v);
    // default dtor works fine
  private: // disable default copy ctor and assignment
    ConstColMatImpl(const ConstColMatImpl&);            // NEVER DEFINED -- disable copy ctor
    ConstColMatImpl& operator=(const ConstColMatImpl&); // NEVER DEFINED -- disable assignment

  public:
    typedef std::vector<RingElem> vec;
    const ring& myRing() const;
    virtual long myNumRows() const;
    virtual long myNumCols() const;
    virtual RingElemAlias myEntry(long i, long j) const;
// Use default defns of the following two functions
//    virtual void myMulByRow(vec& lhs, const vec& v) const;
//    virtual void myMulByCol(vec& lhs, const vec& v) const;
    virtual bool myIsZeroRow(long i) const;
    virtual bool myIsZeroCol(long j) const;
    virtual void myDet(RingElem& d) const;
    virtual long myRank() const;

  private: // data members
    const ring myR;
    long myNumRowsValue;
    const std::vector<RingElem>& myV; // NB this is a reference!!!
  };


  ConstColMatImpl::ConstColMatImpl(const std::vector<RingElem>& v):
      myR(owner(v[0])),
      myNumRowsValue(len(v)),
      myV(v)
  {
    CoCoA_ASSERT(!v.empty());
    CoCoA_ASSERT(HasUniqueOwner(v));
  }


  const ring& ConstColMatImpl::myRing() const
  {
    return myR;
  }


  long ConstColMatImpl::myNumRows() const
  {
    return myNumRowsValue;
  }


  long ConstColMatImpl::myNumCols() const
  {
    return 1;
  }


  RingElemAlias ConstColMatImpl::myEntry(long i, long j) const
  {
    (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    return myV[i];
  }


  bool ConstColMatImpl::myIsZeroRow(long i) const
  {
    CoCoA_ASSERT(i >= 0 && i < myNumRows());
    return IsZero(myV[i]);
  }

  bool ConstColMatImpl::myIsZeroCol(long j) const
  {
    (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(j >= 0 && j < myNumCols());
    for (long i=0; i < myNumRows(); ++i)
      if (!IsZero(myV[i])) return false;
    return true;
  }


  void ConstColMatImpl::myDet(RingElem& d) const
  {
    CoCoA_ASSERT(myNumRows() == 1);
    d = myV[0];
  }


  long ConstColMatImpl::myRank() const
  {
    if (myIsZeroCol(0)) return 0;
    return 1;
  }


  ConstMatrixView ColMat(const std::vector<RingElem>& v)
  {
    if (v.empty())
      CoCoA_ERROR("empty vector of ring elems", "ColMat(v)");
    if (!HasUniqueOwner(v))
      CoCoA_ERROR(ERR::MixedRings, "ColMat(v)");
    return ConstMatrixView(new ConstColMatImpl(v));
  }

  /****************************************************************************/

  class ColMatImpl: public MatrixViewBase
  {
  private:
    friend MatrixView ColMat(std::vector<RingElem>& v); // pseudo-ctor
    ColMatImpl(std::vector<RingElem>& v);
    // default dtor works fine
  private: // disable default copy ctor and assignment
    ColMatImpl(const ColMatImpl&);            // NEVER DEFINED -- disable copy ctor
    ColMatImpl& operator=(const ColMatImpl&); // NEVER DEFINED -- disable assignment

  public:
    typedef std::vector<RingElem> vec;
    const ring& myRing() const;
    virtual long myNumRows() const;
    virtual long myNumCols() const;
    virtual RingElemAlias myEntry(long i, long j) const;
// Use default defns of the following two functions
//    virtual void myMulByRow(vec& lhs, const vec& v) const;
//    virtual void myMulByCol(vec& lhs, const vec& v) const;
    virtual bool myIsZeroRow(long i) const;
    virtual bool myIsZeroCol(long j) const;
    virtual void myDet(RingElem& d) const;
    virtual long myRank() const;

    // non-const member fns
    virtual void myAssignZero();
    virtual bool myIsWritable(long i, long j);
    virtual RingElemRawPtr myRawEntry(long i, long j);
    virtual void mySetEntry(long i, long j, ConstRefRingElem r);
    virtual void mySetEntry(long i, long j, const MachineInt& n);
    virtual void mySetEntry(long i, long j, const BigInt& N);
    virtual void mySetEntry(long i, long j, const BigRat& Q);

  private: // data members
    const ring myR;
    long myNumRowsValue;
    std::vector<RingElem>& myV; // NB this is a reference!!!
  };


  ColMatImpl::ColMatImpl(std::vector<RingElem>& v):
      myR(owner(v[0])),
      myNumRowsValue(len(v)),
      myV(v)
  {
    CoCoA_ASSERT(!v.empty());
    CoCoA_ASSERT(HasUniqueOwner(v));
  }


  const ring& ColMatImpl::myRing() const
  {
    return myR;
  }


  long ColMatImpl::myNumRows() const
  {
    return myNumRowsValue;
  }


  long ColMatImpl::myNumCols() const
  {
    return 1;
  }


  RingElemAlias ColMatImpl::myEntry(long i, long j) const
  {
    (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    return myV[i];
  }


  bool ColMatImpl::myIsZeroRow(long i) const
  {
    CoCoA_ASSERT(i >= 0 && i < myNumRows());
    return IsZero(myV[i]);
  }

  bool ColMatImpl::myIsZeroCol(long j) const
  {
    (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(j >= 0 && j < myNumCols());
    for (long i=0; i < myNumRows(); ++i)
      if (!IsZero(myV[i])) return false;
    return true;
  }


  void ColMatImpl::myDet(RingElem& d) const
  {
    CoCoA_ASSERT(myNumRows() == 1);
    d = myV[0];
  }


  long ColMatImpl::myRank() const
  {
    if (myIsZeroCol(0)) return 0;
    return 1;
  }


  void ColMatImpl::myAssignZero()
  {
    for (long i=0; i < myNumRows(); ++i)
      myV[i] = 0;
  }


  bool ColMatImpl::myIsWritable(long i, long j)
  {
    (void)(i); (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    return true;
  }

  RingElemRawPtr ColMatImpl::myRawEntry(long i, long j)
  {
    (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    return raw(myV[i]);
  }


  void ColMatImpl::mySetEntry(long i, long j, ConstRefRingElem r)
  {
    (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(owner(r) == myR);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    myV[i] = r;
  }

  void ColMatImpl::mySetEntry(long i, long j, const MachineInt& n)
  {
    (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    myV[i] = n;
  }

  void ColMatImpl::mySetEntry(long i, long j, const BigInt& N)
  {
    (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    myV[i] = N;
  }


  void ColMatImpl::mySetEntry(long i, long j, const BigRat& Q)
  {
    (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    myV[i] = Q;
  }


  MatrixView ColMat(std::vector<RingElem>& v)
  {
    if (v.empty())
      CoCoA_ERROR("empty vector of ring elems", "ColMat(v)");
    if (!HasUniqueOwner(v))
      CoCoA_ERROR(ERR::MixedRings, "ColMat(v)");
    return MatrixView(new ColMatImpl(v));
  }


  /***************************************************************************/
  // Very simple definitions... but may not be very efficient.

  MatrixView RowMat(std::vector<RingElem>& v)
  {
    return transpose(ColMat(v));
  }

  ConstMatrixView RowMat(const std::vector<RingElem>& v)
  {
    return transpose(ColMat(v));
  }



  /***************************************************************************/
  // ConstDiagMatImpl -- diagonal matrices implemented in a pretty obvious manner

  class ConstDiagMatImpl: public ConstMatrixViewBase
  {
  public:
    ConstDiagMatImpl(const ring& R, const std::vector<RingElem>& DiagEntries);
    // Default dtor is OK

  public:
    typedef std::vector<RingElem> vec;
    const ring& myRing() const;
    virtual long myNumRows() const;
    virtual long myNumCols() const;
    virtual RingElemAlias myEntry(long i, long j) const;
    virtual void myMulByRow(vec& lhs, const vec& v) const;
    virtual void myMulByCol(vec& lhs, const vec& v) const;
    virtual bool IamEqual(ConstMatrixView M) const;
    virtual bool IamSymmetric() const {return true;}
    virtual bool IamAntiSymmetric() const;
    virtual bool IamDiagonal() const {return true;}
    virtual bool myIsZeroRow(long i) const;
    virtual bool myIsZeroCol(long j) const;
    virtual void myDet(RingElem& d) const;
    virtual long myRank() const;

  private: // data members
    const ring myR;
    const long myDim;
    const std::vector<RingElem>& myDiagEntries; // NB reference to an external vector!!!
  };

  //----------------------------------------------------------------------

  bool IsConstDiagMatImpl(ConstMatrixView M)  
  {return dynamic_cast<ConstDiagMatImpl*>(M.mySmartPtr.myRawPtr()) != 0;}


  ConstDiagMatImpl::ConstDiagMatImpl(const ring& R, const std::vector<RingElem>& DiagEntries):
      myR(R),
      myDim(len(DiagEntries)),
      myDiagEntries(DiagEntries)
  {
#ifdef CoCoA_DEBUG
    for (long i=0; i < myDim; ++i)
      CoCoA_ASSERT(owner(myDiagEntries[i]) == myR);
#endif
  }


  const ring& ConstDiagMatImpl::myRing() const
  {
    return myR;
  }


  long ConstDiagMatImpl::myNumRows() const
  {
    return myDim;
  }


  long ConstDiagMatImpl::myNumCols() const
  {
    return myDim;
  }


  RingElemAlias ConstDiagMatImpl::myEntry(long i, long j) const
  {
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    if (i == j)
      return myDiagEntries[i];
    else
      return zero(myR);
  }


  void ConstDiagMatImpl::myMulByRow(vec& lhs, const vec& v) const
  {
    CoCoA_ASSERT(len(lhs) == myNumCols() && len(v) == myNumRows());
    //BUG: Should also check owner of elements of v and lhs.
    vector<RingElem> ans(myDim, zero(myR));
    for (long i=0; i < myDim; ++i)
      ans[i] = v[i]*myDiagEntries[i];
    for (long i=0; i < myDim; ++i)
      swap(lhs[i], ans[i]);
//???    std::swap(lhs, ans); //??? does this clobber iterators on lhs???
  }


  void ConstDiagMatImpl::myMulByCol(vec& lhs, const vec& v) const
  {
    CoCoA_ASSERT(len(lhs) == myNumRows() && len(v) == myNumCols());
    //BUG: Should also check owner of elements of v and lhs.
    vector<RingElem> ans(myDim, zero(myR));
    for (long i=0; i < myDim; ++i)
      ans[i] = v[i]*myDiagEntries[i];
    for (long i=0; i < myDim; ++i)
      swap(lhs[i], ans[i]);
//???    std::swap(lhs, ans); //??? does this clobber iterators on lhs???
  }


  bool ConstDiagMatImpl::IamEqual(ConstMatrixView M) const
  {
    if (myRing() != RingOf(M)) return false;
    if (myNumRows() != NumRows(M)) return false;
    if (myNumCols() != NumCols(M)) return false;
    if (IsZeroMatImpl(M) || IsIdentityMatImpl(M) || IsDiagMatImpl(M) || IsConstDiagMatImpl(M))
    {
      //      std::cout << "ConstDiagMatImpl::IamEqual" << std::endl;
      for (long i=0; i < myNumRows(); ++i)
        if (myDiagEntries[i] != M(i,i)) return false;
      return true;
    }
    return ConstMatrixViewBase::IamEqual(M);
  }


  bool ConstDiagMatImpl::IamAntiSymmetric() const
  {
    for (long i=0; i < myNumRows(); ++i)
      if (!IsZero(myDiagEntries[i])) return false;
    return true;
  }


  bool ConstDiagMatImpl::myIsZeroRow(long i) const
  {
    CoCoA_ASSERT(i < myNumRows());
    return IsZero(myDiagEntries[i]);
  }


  bool ConstDiagMatImpl::myIsZeroCol(long j) const
  {
    CoCoA_ASSERT(j >= 0 && j < myNumCols());
    return IsZero(myDiagEntries[j]);
  }


  void ConstDiagMatImpl::myDet(RingElem& d) const
  {
    CoCoA_ASSERT(owner(d) == myR);
    RingElem ans = one(myR);
    for (long i=0; i < myDim; ++i)
      ans *= myDiagEntries[i];
    swap(d, ans);
  }


  // ??? better to use "count"?
  long ConstDiagMatImpl::myRank() const
  {
    CoCoA_ASSERT(IsIntegralDomain(myR));
    long rk=0;
    for (long i = 0; i < myDim; ++i)
      if (!IsZero(myDiagEntries[i]))
        ++rk;
    return rk;
  }


  ConstMatrixView DiagMat(const std::vector<RingElem>& DiagEntries)
  {
    if (DiagEntries.empty())
      CoCoA_ERROR("Empty vector", "DiagMat(v)");
    if (!HasUniqueOwner(DiagEntries))
      CoCoA_ERROR(ERR::MixedRings, "DiagMat(v)");
    return ConstMatrixView(new ConstDiagMatImpl(owner(DiagEntries[0]), DiagEntries));
  }


  /***************************************************************************/
  // DiagMatImpl -- diagonal matrices implemented in a pretty obvious manner

  class DiagMatImpl: public MatrixViewBase
  {
  public:
    DiagMatImpl(const ring& R, std::vector<RingElem>& DiagEntries);
    // Default dtor is OK

  public:
    typedef std::vector<RingElem> vec;
    const ring& myRing() const;
    virtual long myNumRows() const;
    virtual long myNumCols() const;
    virtual RingElemAlias myEntry(long i, long j) const;
    virtual bool IamEqual(ConstMatrixView M) const;
    virtual bool IamSymmetric() const {return true;}
    virtual bool IamAntiSymmetric() const;
    virtual bool IamDiagonal() const {return true;}
    virtual void myMulByRow(vec& lhs, const vec& v) const;
    virtual void myMulByCol(vec& lhs, const vec& v) const;
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
    const long myDim;
    std::vector<RingElem>& myDiagEntries; // NB reference to an external vector!!!
  };

  //----------------------------------------------------------------------

  bool IsDiagMatImpl(ConstMatrixView M)  
  {return dynamic_cast<DiagMatImpl*>(M.mySmartPtr.myRawPtr()) != 0;}


  DiagMatImpl::DiagMatImpl(const ring& R, std::vector<RingElem>& DiagEntries):
      myR(R),
      myDim(len(DiagEntries)),
      myDiagEntries(DiagEntries)
  {
#ifdef CoCoA_DEBUG
    for (long i=0; i < myDim; ++i)
      CoCoA_ASSERT(owner(myDiagEntries[i]) == myR);
#endif
  }


  const ring& DiagMatImpl::myRing() const
  {
    return myR;
  }


  long DiagMatImpl::myNumRows() const
  {
    return myDim;
  }


  long DiagMatImpl::myNumCols() const
  {
    return myDim;
  }


  RingElemAlias DiagMatImpl::myEntry(long i, long j) const
  {
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    if (i == j)
      return myDiagEntries[i];
    else
      return zero(myR);
  }


  void DiagMatImpl::myMulByRow(vec& lhs, const vec& v) const
  {
    CoCoA_ASSERT(len(lhs) == myNumCols() && len(v) == myNumRows());
    // Should also check owner of elements of v and lhs.
    vector<RingElem> ans(myDim, zero(myR));
    for (long i=0; i < myDim; ++i)
      ans[i] = v[i]*myDiagEntries[i];
    for (long i=0; i < myDim; ++i)
      swap(lhs[i], ans[i]);
//???    std::swap(lhs, ans); //??? does this clobber iterators on lhs???
  }


  void DiagMatImpl::myMulByCol(vec& lhs, const vec& v) const
  {
    CoCoA_ASSERT(len(lhs) == myNumRows() && len(v) == myNumCols());
    // Should also check owner of elements of v and lhs.
    vector<RingElem> ans(myDim, zero(myR));
    for (long i=0; i < myDim; ++i)
      ans[i] = v[i]*myDiagEntries[i];
    for (long i=0; i < myDim; ++i)
      swap(lhs[i], ans[i]);
//???    std::swap(lhs, ans); //??? would this clobber iterators on lhs???
  }


  bool DiagMatImpl::IamEqual(ConstMatrixView M) const
  {
    if (myRing() != RingOf(M)) return false;
    if (myNumRows() != NumRows(M)) return false;
    if (myNumCols() != NumCols(M)) return false;
    if (IsZeroMatImpl(M) || IsIdentityMatImpl(M) || IsDiagMatImpl(M) || IsConstDiagMatImpl(M))
    {
      //      std::cout << "DiagMatImpl::IamEqual" << std::endl;
      for (long i=0; i < myNumRows(); ++i)
        if (myDiagEntries[i] != M(i,i)) return false;
      return true;
    }
    return ConstMatrixViewBase::IamEqual(M);
  }


  bool DiagMatImpl::IamAntiSymmetric() const
  {
    for (long i=0; i < myNumRows(); ++i)
      if (!IsZero(myDiagEntries[i])) return false;
    return true;
  }


  bool DiagMatImpl::myIsZeroRow(long i) const
  {
    CoCoA_ASSERT(i < myNumRows());
    return IsZero(myDiagEntries[i]);
  }


  bool DiagMatImpl::myIsZeroCol(long j) const
  {
    CoCoA_ASSERT(j >= 0 && j < myNumCols());
    return IsZero(myDiagEntries[j]);
  }


  void DiagMatImpl::myDet(RingElem& d) const
  {
    CoCoA_ASSERT(owner(d) == myR);
    RingElem ans = one(myR);
    for (long i=0; i < myDim; ++i)
      ans *= myDiagEntries[i];
    swap(d, ans);
  }


  // ??? better to use "std::count"???
  long DiagMatImpl::myRank() const
  {
    CoCoA_ASSERT(IsIntegralDomain(myR));
    long rk=0;
    for (long i = 0; i < myDim; ++i)
      if (!IsZero(myDiagEntries[i]))
        ++rk;
    return rk;
  }


  // Exception clean only if the call to fill cannot throw
  void DiagMatImpl::myAssignZero()
  {
    fill(myDiagEntries.begin(), myDiagEntries.end(), zero(myR));
  }


  bool DiagMatImpl::myIsWritable(long i, long j)
  {
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    return i==j;
  }

  RingElemRawPtr DiagMatImpl::myRawEntry(long i, long j)
  {
    (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(myIsWritable(i,j));
    return raw(myDiagEntries[i]);
  }


  void DiagMatImpl::mySetEntry(long i, long j, ConstRefRingElem r)
  {
    (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(myIsWritable(i,j));
    myDiagEntries[i] = r;
  }


  void DiagMatImpl::mySetEntry(long i, long j, const MachineInt& n)
  {
    (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(myIsWritable(i,j));
    myDiagEntries[i] = n;
  }


  void DiagMatImpl::mySetEntry(long i, long j, const BigInt& N)
  {
    (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(myIsWritable(i,j));
    myDiagEntries[i] = N;
  }


  void DiagMatImpl::mySetEntry(long i, long j, const BigRat& Q)
  {
    (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(myIsWritable(i,j));
    myDiagEntries[i] = Q;
  }


  MatrixView DiagMat(std::vector<RingElem>& DiagEntries)
  {
    if (DiagEntries.empty())
      CoCoA_ERROR("Empty vector", "DiagMat(v)");
    if (!HasUniqueOwner(DiagEntries))
      CoCoA_ERROR(ERR::MixedRings, "DiagMat(v)");
    return MatrixView(new DiagMatImpl(owner(DiagEntries[0]), DiagEntries));
  }

  /****************************************************************************/

  class MatByRowsImpl: public MatrixViewBase
  {
  private:
    friend MatrixView MatByRows(long nrows, long ncols, std::vector<RingElem>& v); // pseudo-ctor
    MatByRowsImpl(long nrows, long ncols, std::vector<RingElem>& v);
    static std::vector<RingElem>& CtorArgCheck(long nrows, long ncols, std::vector<RingElem>& v);
    // default dtor works fine
  private: // disable default copy ctor and assignment
    MatByRowsImpl(const MatByRowsImpl&);            // NEVER DEFINED -- disable copy ctor
    MatByRowsImpl& operator=(const MatByRowsImpl&); // NEVER DEFINED -- disable assignment

  public:
    typedef std::vector<RingElem> vec;
    const ring& myRing() const;
    virtual long myNumRows() const;
    virtual long myNumCols() const;
    virtual RingElemAlias myEntry(long i, long j) const;
// Use default defns of the following two functions
//    virtual void myMulByRow(vec& lhs, const vec& v) const;
//    virtual void myMulByCol(vec& lhs, const vec& v) const;
    virtual bool myIsZeroRow(long i) const;
    virtual bool myIsZeroCol(long j) const;
/// USE DEFAULT?    virtual void myDet(RingElem& d) const;
/// USE DEFAULT?    virtual long myRank() const;

    // non-const member fns
    virtual void myAssignZero();
    virtual bool myIsWritable(long i, long j);
    virtual RingElemRawPtr myRawEntry(long i, long j);
    virtual void mySetEntry(long i, long j, ConstRefRingElem r);
    virtual void mySetEntry(long i, long j, const MachineInt& n);
    virtual void mySetEntry(long i, long j, const BigInt& N);
    virtual void mySetEntry(long i, long j, const BigRat& Q);

  private: // data members
    const long myNumRowsValue;
    const long myNumColsValue;
    const ring myR;
    std::vector<RingElem>& myV; // NB this is a reference!!!
  private: // impl detail
    long myIndex(long i, long j) const { return i*myNumColsValue+j; }

  };


  std::vector<RingElem>& MatByRowsImpl::CtorArgCheck(long nrows, long ncols, std::vector<RingElem>& v)
  {
    if (nrows < 1 || ncols < 1) CoCoA_ERROR(ERR::NotPositive, "MatByRows ctor");
    const long n = len(v);
    if (n/nrows != ncols || n != nrows*ncols) CoCoA_ERROR(ERR::BadArraySize, "MatByRows ctor");
//???    if (v.empty()) CoCoA_ERROR(ERR::Empty, "MatByRows ctor");
    if (!HasUniqueOwner(v))
      CoCoA_ERROR(ERR::MixedRings, "MatByRows ctor");
    return v;
  }


  MatByRowsImpl::MatByRowsImpl(long nrows, long ncols, std::vector<RingElem>& v):
      myNumRowsValue(nrows),
      myNumColsValue(ncols),
      myR(owner(CtorArgCheck(nrows, ncols,v)[0])),
      myV(v)
  {}


  const ring& MatByRowsImpl::myRing() const
  {
    return myR; // return owner(myV[0]);
  }


  long MatByRowsImpl::myNumRows() const
  {
    return myNumRowsValue;
  }


  long MatByRowsImpl::myNumCols() const
  {
    return myNumColsValue;
  }


  RingElemAlias MatByRowsImpl::myEntry(long i, long j) const
  {
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    return myV[myIndex(i,j)];
  }


  bool MatByRowsImpl::myIsZeroRow(long i) const
  {
    CoCoA_ASSERT(i >= 0 && i < myNumRows());
    const long start = myIndex(i,0);
    const long end = myIndex(i,myNumColsValue-1);
    for (long j=start; j <= end; ++j)
      if (!IsZero(myV[j])) return false;
    return true;
  }

  bool MatByRowsImpl::myIsZeroCol(long j) const
  {
    CoCoA_ASSERT(j >= 0 && j < myNumCols());
    for (long i=0; i < myNumRows(); ++i)
      if (!IsZero(myV[myIndex(i,j)])) return false;
    return true;
  }


//   void MatByRowsImpl::myDet(RingElem& d) const
//   {
//     CoCoA_ASSERT(myNumRows() == myNumCols());
//     d = myV[0];
//   }


//   long MatByRowsImpl::myRank() const
//   {
//     if (myIsZeroCol(0)) return 0;
//     return 1;
//   }


  void MatByRowsImpl::myAssignZero()
  {
    const long n = len(myV);
    for (long i=0; i < n; ++i)
      myV[i] = 0;
  }


  bool MatByRowsImpl::myIsWritable(long i, long j)
  {
    (void)(i); (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    return true;
  }

  RingElemRawPtr MatByRowsImpl::myRawEntry(long i, long j)
  {
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    return raw(myV[myIndex(i,j)]);
  }


  void MatByRowsImpl::mySetEntry(long i, long j, ConstRefRingElem r)
  {
    CoCoA_ASSERT(owner(r) == myR);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    myV[myIndex(i,j)] = r;
  }

  void MatByRowsImpl::mySetEntry(long i, long j, const MachineInt& n)
  {
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    myV[myIndex(i,j)] = n;
  }

  void MatByRowsImpl::mySetEntry(long i, long j, const BigInt& N)
  {
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    myV[myIndex(i,j)] = N;
  }


  void MatByRowsImpl::mySetEntry(long i, long j, const BigRat& Q)
  {
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    myV[myIndex(i,j)] = Q;
  }


  MatrixView MatByRows(long nrows, long ncols, std::vector<RingElem>& v)
  {
    if (nrows <= 0 || ncols <= 0)
      CoCoA_ERROR(ERR::NotPositive, "MatByRows pseudo-ctor");
    // We let ctor for MatByRowsImpl do the arg checking...
//     if (vv.empty())
//       CoCoA_ERROR("empty vector of ring elems", "MatByRows(vv)");
//     const long nrows = len(vv);
//     for (long i=0; i < nrows; ++i)
//       if (!HasUniqueOwner(v[i]) || owner(vv[i][0]) != owner(vv[0][0]))
//       CoCoA_ERROR(ERR::MixedRings, "MatByRows(v)");
    return MatrixView(new MatByRowsImpl(nrows,ncols,v));
  }


  /****************************************************************************/

  class MatByRowsImpl2: public MatrixViewBase
  {
  private:
    friend MatrixView MatByRows(std::vector< std::vector<RingElem> >& vv); // pseudo-ctor
    MatByRowsImpl2(std::vector< std::vector<RingElem> >& vv);
    static std::vector< std::vector<RingElem> >& CtorArgCheck(std::vector< std::vector<RingElem> >& vv);
    // default dtor works fine
  private: // disable default copy ctor and assignment
    MatByRowsImpl2(const MatByRowsImpl&);            // NEVER DEFINED -- disable copy ctor
    MatByRowsImpl2& operator=(const MatByRowsImpl2&); // NEVER DEFINED -- disable assignment

  public:
    typedef std::vector<RingElem> vec;
    const ring& myRing() const;
    virtual long myNumRows() const;
    virtual long myNumCols() const;
    virtual RingElemAlias myEntry(long i, long j) const;
// Use default defns of the following two functions
//    virtual void myMulByRow(vec& lhs, const vec& v) const;
//    virtual void myMulByCol(vec& lhs, const vec& v) const;
    virtual bool myIsZeroRow(long i) const;
    virtual bool myIsZeroCol(long j) const;
/// USE DEFAULT?    virtual void myDet(RingElem& d) const;
/// USE DEFAULT?    virtual long myRank() const;

    // non-const member fns
    virtual void myAssignZero();
    virtual bool myIsWritable(long i, long j);
    virtual RingElemRawPtr myRawEntry(long i, long j);
    virtual void mySetEntry(long i, long j, ConstRefRingElem r);
    virtual void mySetEntry(long i, long j, const MachineInt& n);
    virtual void mySetEntry(long i, long j, const BigInt& N);
    virtual void mySetEntry(long i, long j, const BigRat& Q);

  private: // data members
    const ring myR;
    const long myNumRowsValue;
    const long myNumColsValue;
    std::vector< std::vector<RingElem> >& myVV; // NB this is a reference!!!
  };


  std::vector< std::vector<RingElem> >& MatByRowsImpl2::CtorArgCheck(std::vector< std::vector<RingElem> >& vv)
  {
    if (vv.empty()) CoCoA_ERROR(ERR::Empty, "MatByRows ctor");
    const long nrows = len(vv);
    for (long i=0; i < nrows; ++i)
      if (!HasUniqueOwner(vv[i]) ||
          len(vv[i]) != len(vv[0]) ||
          owner(vv[i][0]) != owner(vv[0][0]))
      CoCoA_ERROR(ERR::MixedRings, "MatByRows ctor");
    return vv;
  }


  MatByRowsImpl2::MatByRowsImpl2(std::vector< std::vector<RingElem> >& vv):
      myR(owner(CtorArgCheck(vv)[0][0])),
      myNumRowsValue(len(vv)),
      myNumColsValue(len(vv[0])),
      myVV(vv)
  {}


  const ring& MatByRowsImpl2::myRing() const
  {
    return myR;
  }


  long MatByRowsImpl2::myNumRows() const
  {
    return myNumRowsValue;
  }


  long MatByRowsImpl2::myNumCols() const
  {
    return myNumColsValue;
  }


  RingElemAlias MatByRowsImpl2::myEntry(long i, long j) const
  {
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    return myVV[i][j];
  }


  bool MatByRowsImpl2::myIsZeroRow(long i) const
  {
    CoCoA_ASSERT(i >= 0 && i < myNumRows());
    for (long j=0; j < myNumColsValue; ++j)
      if (!IsZero(myVV[i][j])) return false;
    return true;
  }

  bool MatByRowsImpl2::myIsZeroCol(long j) const
  {
    CoCoA_ASSERT(j >= 0 && j < myNumCols());
    for (long i=0; i < myNumRows(); ++i)
      if (!IsZero(myVV[i][j])) return false;
    return true;
  }


//   void MatByRowsImpl2::myDet(RingElem& d) const
//   {
//     CoCoA_ASSERT(myNumRows() == myNumCols());
//     d = myV[0];
//   }


//   long MatByRowsImpl2::myRank() const
//   {
//     if (myIsZeroCol(0)) return 0;
//     return 1;
//   }


  void MatByRowsImpl2::myAssignZero()
  {
    for (long i=0; i < myNumRows(); ++i)
      for (long j=0; j < myNumCols(); ++j)
        myVV[i][j] = 0;
  }


  bool MatByRowsImpl2::myIsWritable(long i, long j)
  {
    (void)(i); (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    return true;
  }

  RingElemRawPtr MatByRowsImpl2::myRawEntry(long i, long j)
  {
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    return raw(myVV[i][j]);
  }


  void MatByRowsImpl2::mySetEntry(long i, long j, ConstRefRingElem r)
  {
    CoCoA_ASSERT(owner(r) == myR);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    myVV[i][j] = r;
  }

  void MatByRowsImpl2::mySetEntry(long i, long j, const MachineInt& n)
  {
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    myVV[i][j] = n;
  }

  void MatByRowsImpl2::mySetEntry(long i, long j, const BigInt& N)
  {
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    myVV[i][j] = N;
  }


  void MatByRowsImpl2::mySetEntry(long i, long j, const BigRat& Q)
  {
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    myVV[i][j] = Q;
  }


  MatrixView MatByRows(std::vector< std::vector<RingElem> >& vv)
  {
    // We let ctor for MatByRowsImpl2 do the arg checking...
//     if (vv.empty())
//       CoCoA_ERROR("empty vector of ring elems", "MatByRows(vv)");
//     const long nrows = len(vv);
//     for (long i=0; i < nrows; ++i)
//       if (!HasUniqueOwner(v[i]) || owner(vv[i][0]) != owner(vv[0][0]))
//       CoCoA_ERROR(ERR::MixedRings, "MatByRows(v)");
    return MatrixView(new MatByRowsImpl2(vv));
  }


  /****************************************************************************/

  class MatByColsImpl: public MatrixViewBase
  {
  private:
    friend MatrixView MatByCols(long nrows, long ncols, std::vector<RingElem>& v); // pseudo-ctor
    MatByColsImpl(long nrows, long ncols, std::vector<RingElem>& v);
    static std::vector<RingElem>& CtorArgCheck(long nrows, long ncols, std::vector<RingElem>& v);
    // default dtor works fine
  private: // disable default copy ctor and assignment
    MatByColsImpl(const MatByRowsImpl&);            // NEVER DEFINED -- disable copy ctor
    MatByColsImpl& operator=(const MatByRowsImpl&); // NEVER DEFINED -- disable assignment

  public:
    typedef std::vector<RingElem> vec;
    const ring& myRing() const;
    virtual long myNumRows() const;
    virtual long myNumCols() const;
    virtual RingElemAlias myEntry(long i, long j) const;
// Use default defns of the following two functions
//    virtual void myMulByRow(vec& lhs, const vec& v) const;
//    virtual void myMulByCol(vec& lhs, const vec& v) const;
    virtual bool myIsZeroRow(long i) const;
    virtual bool myIsZeroCol(long j) const;
/// USE DEFAULT?    virtual void myDet(RingElem& d) const;
/// USE DEFAULT?    virtual long myRank() const;

    // non-const member fns
    virtual void myAssignZero();
    virtual bool myIsWritable(long i, long j);
    virtual RingElemRawPtr myRawEntry(long i, long j);
    virtual void mySetEntry(long i, long j, ConstRefRingElem r);
    virtual void mySetEntry(long i, long j, const MachineInt& n);
    virtual void mySetEntry(long i, long j, const BigInt& N);
    virtual void mySetEntry(long i, long j, const BigRat& Q);

  private: // data members
    const long myNumRowsValue;
    const long myNumColsValue;
    const ring myR;
    std::vector<RingElem>& myV; // NB this is a reference!!!
  private: // impl detail
    long myIndex(long i, long j) const { return j*myNumRowsValue+i; }

  };


  std::vector<RingElem>& MatByColsImpl::CtorArgCheck(long nrows, long ncols, std::vector<RingElem>& v)
  {
    if (nrows < 1 || ncols < 1) CoCoA_ERROR(ERR::NotPositive, "MatByCols ctor");
    const long n = len(v);
    // Check that n = nrows*ncols (avoiding overflow)
    if (n/nrows != ncols || n%nrows != 0) CoCoA_ERROR(ERR::BadArraySize, "MatByCols ctor");
//???    if (v.empty()) CoCoA_ERROR(ERR::Empty, "MatByCols ctor");
    if (!HasUniqueOwner(v))
      CoCoA_ERROR(ERR::MixedRings, "MatByCols ctor");
    return v;
  }


  MatByColsImpl::MatByColsImpl(long nrows, long ncols, std::vector<RingElem>& v):
      myNumRowsValue(nrows),
      myNumColsValue(ncols),
      myR(owner(CtorArgCheck(nrows, ncols,v)[0])),
      myV(v)
  {}


  const ring& MatByColsImpl::myRing() const
  {
    return myR; // return owner(myV[0]);
  }


  long MatByColsImpl::myNumRows() const
  {
    return myNumRowsValue;
  }


  long MatByColsImpl::myNumCols() const
  {
    return myNumColsValue;
  }


  RingElemAlias MatByColsImpl::myEntry(long i, long j) const
  {
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    return myV[myIndex(i,j)];
  }


  bool MatByColsImpl::myIsZeroRow(long i) const
  {
    CoCoA_ASSERT(i >= 0 && i < myNumRows());
    for (long j=0; j < myNumCols(); ++j)
      if (!IsZero(myV[myIndex(i,j)])) return false;
    return true;
  }

  bool MatByColsImpl::myIsZeroCol(long j) const
  {
    CoCoA_ASSERT(j >= 0 && j < myNumCols());
    const long start = myIndex(0,j);
    const long end = myIndex(myNumRowsValue-1,j);
    for (long i=start; i <= end; ++i)
      if (!IsZero(myV[i])) return false;
    return true;
  }


//   void MatByColsImpl::myDet(RingElem& d) const
//   {
//     CoCoA_ASSERT(myNumRows() == myNumCols());
//     d = myV[0];
//   }


//   long MatByColsImpl::myRank() const
//   {
//     if (myIsZeroCol(0)) return 0;
//     return 1;
//   }


  void MatByColsImpl::myAssignZero()
  {
    const long n = len(myV);
    for (long i=0; i < n; ++i)
      myV[i] = 0;
  }


  bool MatByColsImpl::myIsWritable(long i, long j)
  {
    (void)(i); (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    return true;
  }

  RingElemRawPtr MatByColsImpl::myRawEntry(long i, long j)
  {
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    return raw(myV[myIndex(i,j)]);
  }


  void MatByColsImpl::mySetEntry(long i, long j, ConstRefRingElem r)
  {
    CoCoA_ASSERT(owner(r) == myR);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    myV[myIndex(i,j)] = r;
  }

  void MatByColsImpl::mySetEntry(long i, long j, const MachineInt& n)
  {
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    myV[myIndex(i,j)] = n;
  }

  void MatByColsImpl::mySetEntry(long i, long j, const BigInt& N)
  {
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    myV[myIndex(i,j)] = N;
  }


  void MatByColsImpl::mySetEntry(long i, long j, const BigRat& Q)
  {
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    myV[myIndex(i,j)] = Q;
  }


  MatrixView MatByCols(long nrows, long ncols, std::vector<RingElem>& v)
  {
    if (nrows <= 0 || ncols <= 0)
      CoCoA_ERROR(ERR::NotPositive, "MatByCols pseudo-ctor");
    // We let ctor for MatByColsImpl do the arg checking...
//     if (vv.empty())
//       CoCoA_ERROR("empty vector of ring elems", "MatByCols(vv)");
//     const long nrows = len(vv);
//     for (long i=0; i < nrows; ++i)
//       if (!HasUniqueOwner(v[i]) || owner(vv[i][0]) != owner(vv[0][0]))
//       CoCoA_ERROR(ERR::MixedRings, "MatByCols(v)");
    return MatrixView(new MatByColsImpl(nrows,ncols,v));
  }


  /***************************************************************************/

  class ConstBlockMat2x2Impl: public ConstMatrixViewBase
  {
  private:
    friend ConstMatrixView BlockMat2x2(ConstMatrixView A, ConstMatrixView B, ConstMatrixView C, ConstMatrixView D); // pseudo-ctor
    ConstBlockMat2x2Impl(ConstMatrixView A, ConstMatrixView B, ConstMatrixView C, ConstMatrixView D);
    // default dtor works fine
  private: // disable default copy ctor and assignment
    ConstBlockMat2x2Impl(const ConstBlockMat2x2Impl&);            // NEVER DEFINED -- disable copy ctor
    ConstBlockMat2x2Impl& operator=(const ConstBlockMat2x2Impl&); // NEVER DEFINED -- disable assignment

  public:
    typedef std::vector<RingElem> vec;
    const ring& myRing() const;
    virtual long myNumRows() const;
    virtual long myNumCols() const;
    virtual RingElemAlias myEntry(long i, long j) const;
// Use default defns of the following two functions
//    virtual void myMulByRow(vec& lhs, const vec& v) const;
//    virtual void myMulByCol(vec& lhs, const vec& v) const;
    virtual bool myIsZeroRow(long i) const;
    virtual bool myIsZeroCol(long j) const;
// Use default defns of the following two functions
//    virtual void myDet(RingElem& d) const;
//    virtual long myRank() const;

  private: // data members
    const ring myR;
    long myNumRowsValue;
    long myNumColsValue;
    ConstMatrixView myA, myB, myC, myD;
  };


  /***************************************************************************/
  // Functions for ConstBlockMat2x2Impl -- a matrix of the form  (A B)
  //                                                             (C D)
  // This implementation requires that the submatrices correspond
  // to two straight cuts -- several changes would be needed to handle
  // a more general assembly of submatrices.

  ConstBlockMat2x2Impl::ConstBlockMat2x2Impl(ConstMatrixView A, ConstMatrixView B, ConstMatrixView C, ConstMatrixView D):
      myR(RingOf(A)),
      myNumRowsValue(NumRows(A)+NumRows(C)),
      myNumColsValue(NumCols(A)+NumCols(B)),
      myA(A),
      myB(B),
      myC(C),
      myD(D)
  {
    // Args are checked for compatibility by pseudo-ctor BlockMat2x2, the
    // only function which calls this ctor.
  }


  const ring& ConstBlockMat2x2Impl::myRing() const
  {
    return myR;
  }


  long ConstBlockMat2x2Impl::myNumRows() const
  {
    return myNumRowsValue;
  }


  long ConstBlockMat2x2Impl::myNumCols() const
  {
    return myNumColsValue;
  }


  RingElemAlias ConstBlockMat2x2Impl::myEntry(long i, long j) const
  {
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    if (i < NumRows(myA) && j < NumCols(myA))
      return myA(i, j);
    if (i < NumRows(myB))
      return myB(i, j-NumCols(myA));
    if (j < NumCols(myC))
      return myC(i-NumRows(myA), j);
    return myD(i-NumRows(myA), j-NumCols(myA));
  }


  bool ConstBlockMat2x2Impl::myIsZeroRow(long i) const
  {
    CoCoA_ASSERT(i < myNumRows());
    if (i < NumRows(myA))
      return myA->myIsZeroRow(i) && myB->myIsZeroRow(i);
    return myC->myIsZeroRow(i-NumRows(myA)) && myD->myIsZeroRow(i-NumRows(myA));
  }


  bool ConstBlockMat2x2Impl::myIsZeroCol(long j) const
  {
    CoCoA_ASSERT(j >= 0 && j < myNumCols());
    if (j < NumCols(myA))
      return myA->myIsZeroCol(j) && myC->myIsZeroCol(j);
    return myB->myIsZeroCol(j-NumCols(myA)) && myD->myIsZeroCol(j-NumCols(myA));
  }


  ConstMatrixView BlockMat2x2(ConstMatrixView A, ConstMatrixView B, ConstMatrixView C, ConstMatrixView D)
  {
    {
      // Check that the rings are all the same
      const ring R = RingOf(A);
      if (RingOf(B) != R || RingOf(C) != R || RingOf(D) != R)
        CoCoA_ERROR(ERR::MixedRings, "BlockMat2x2(A,B,C,D)");
    }
    // Check that the matrix sizes are suitably matched.
    if (NumRows(A) != NumRows(B) ||
        NumRows(C) != NumRows(D) ||
        NumCols(A) != NumCols(C) ||
        NumCols(B) != NumCols(D))
      CoCoA_ERROR(ERR::BadMatrixSize, "BlockMat2x2(A,B,C,D)");
    if (NumCols(A) > numeric_limits<long>::max() -NumCols(B)) //avoid overflow  
      CoCoA_ERROR(ERR::BadColIndex, "BlockMat2x2(A,B,C,D)");
    if (NumRows(A) > numeric_limits<long>::max() -NumRows(C)) //avoid overflow  
      CoCoA_ERROR(ERR::BadRowIndex, "BlockMat2x2(A,B,C,D)");
    return ConstMatrixView(new ConstBlockMat2x2Impl(A, B, C, D));
  }


/***************************************************************************/

  class BlockMat2x2Impl: public MatrixViewBase
  {
  private:
    friend MatrixView BlockMat2x2(MatrixView A, MatrixView B, MatrixView C, MatrixView D); // pseudo-ctor
    BlockMat2x2Impl(MatrixView A, MatrixView B, MatrixView C, MatrixView D);
    // default dtor works fine
  private: // disable default copy ctor and assignment
    BlockMat2x2Impl(const BlockMat2x2Impl&);            // NEVER DEFINED -- disable copy ctor
    BlockMat2x2Impl& operator=(const BlockMat2x2Impl&); // NEVER DEFINED -- disable assignment

  public:
    typedef std::vector<RingElem> vec;
    const ring& myRing() const;
    virtual long myNumRows() const;
    virtual long myNumCols() const;
    virtual RingElemAlias myEntry(long i, long j) const;
// Use default defns of the following two functions
//    virtual void myMulByRow(vec& lhs, const vec& v) const;
//    virtual void myMulByCol(vec& lhs, const vec& v) const;
    virtual bool myIsZeroRow(long i) const;
    virtual bool myIsZeroCol(long j) const;
// Use default defns of the following two functions
//    virtual void myDet(RingElem& d) const;
//    virtual long myRank() const;

    // Non-const member functions
    virtual void myAssignZero();
    virtual bool myIsWritable(long i, long j);
    virtual RingElemRawPtr myRawEntry(long i, long j);
    virtual void mySetEntry(long i, long j, ConstRefRingElem r);
    virtual void mySetEntry(long i, long j, const MachineInt& n);
    virtual void mySetEntry(long i, long j, const BigInt& N);
    virtual void mySetEntry(long i, long j, const BigRat& Q);

  private: // data members
    const ring myR;
    long myNumRowsValue;
    long myNumColsValue;
    MatrixView myA, myB, myC, myD;
  };


  /***************************************************************************/
  // Functions for BlockMat2x2Impl -- a matrix of the form  (A B)
  //                                                        (C D)
  // This implementation requires that the submatrices correspond
  // to two straight cuts -- several changes would be needed to handle
  // a more general assembly of submatrices.

  BlockMat2x2Impl::BlockMat2x2Impl(MatrixView A, MatrixView B, MatrixView C, MatrixView D):
      myR(RingOf(A)),
      myNumRowsValue(NumRows(A)+NumRows(C)),
      myNumColsValue(NumCols(A)+NumCols(B)),
      myA(A),
      myB(B),
      myC(C),
      myD(D)
  {
    // Args are checked for compatibility by pseudo-ctor BlockMat2x2, the
    // only function which calls this ctor.
  }


  const ring& BlockMat2x2Impl::myRing() const
  {
    return myR;
  }


  long BlockMat2x2Impl::myNumRows() const
  {
    return myNumRowsValue;
  }


  long BlockMat2x2Impl::myNumCols() const
  {
    return myNumColsValue;
  }


  RingElemAlias BlockMat2x2Impl::myEntry(long i, long j) const
  {
// This one-liner doesn't work:    return ConstRefRingElem(myR, myRawEntry(i,j));
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    if (i < NumRows(myA) && j < NumCols(myA))
      return myA(i, j);
    if (i < NumRows(myB))
      return myB(i, j-NumCols(myA));
    if (j < NumCols(myC))
      return myC(i-NumRows(myA), j);
    return myD(i-NumRows(myA), j-NumCols(myA));
  }


  bool BlockMat2x2Impl::myIsZeroRow(long i) const
  {
    CoCoA_ASSERT(i < myNumRows());
    if (i < NumRows(myA))
      return myA->myIsZeroRow(i) && myB->myIsZeroRow(i);
    return myC->myIsZeroRow(i-NumRows(myA)) && myD->myIsZeroRow(i-NumRows(myA));
  }


  bool BlockMat2x2Impl::myIsZeroCol(long j) const
  {
    CoCoA_ASSERT(j >= 0 && j < myNumCols());
    if (j < NumCols(myA))
      return myA->myIsZeroCol(j) && myC->myIsZeroCol(j);
    return myB->myIsZeroCol(j-NumCols(myA)) && myD->myIsZeroCol(j-NumCols(myA));
  }


  bool BlockMat2x2Impl::myIsWritable(long i, long j)
  {
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    if (j >= NumCols(myA))
    {
      if (i >= NumRows(myB))
        return myD->myIsWritable(i-NumRows(myA), j-NumCols(myB));
      else
        return myB->myIsWritable(i, j-NumCols(myA));
    }
    if (i >= NumRows(myA))
      return myC->myIsWritable(i-NumRows(myA), j);
    else
      return myA->myIsWritable(i, j);
  }


  void BlockMat2x2Impl::myAssignZero()
  {
    myA->myAssignZero();
    myB->myAssignZero();
    myC->myAssignZero();
    myD->myAssignZero();
  }


  RingElemRawPtr BlockMat2x2Impl::myRawEntry(long i, long j)
  {
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    if (i < NumRows(myA) && j < NumCols(myA))
      return myA->myRawEntry(i, j);
    if (i < NumRows(myB))
      return myB->myRawEntry(i, j-NumCols(myA));
    if (j < NumCols(myC))
      return myC->myRawEntry(i-NumRows(myA), j);
    return myD->myRawEntry(i-NumRows(myA), j-NumCols(myA));
  }


  void BlockMat2x2Impl::mySetEntry(long i, long j, ConstRefRingElem r)
  {
    CoCoA_ASSERT(myR == owner(r));
    CoCoA_ASSERT(myIsWritable(i,j));
    myR->myAssign(myRawEntry(i,j), raw(r)); // myRefEntry(i,j) = r;
  }

  void BlockMat2x2Impl::mySetEntry(long i, long j, const MachineInt& n)
  {
    CoCoA_ASSERT(myIsWritable(i,j));
    myR->myAssign(myRawEntry(i,j), n); //    myRefEntry(i,j) = n;
  }

  void BlockMat2x2Impl::mySetEntry(long i, long j, const BigInt& N)
  {
    CoCoA_ASSERT(myIsWritable(i,j));
    myR->myAssign(myRawEntry(i,j), N); //    myRefEntry(i,j) = N;
  }


  void BlockMat2x2Impl::mySetEntry(long i, long j, const BigRat& Q)
  {
    CoCoA_ASSERT(myIsWritable(i,j));
    myR->myAssign(myRawEntry(i,j), Q); //    myRefEntry(i,j) = Q;
  }


  MatrixView BlockMat2x2(MatrixView A, MatrixView B, MatrixView C, MatrixView D)
  {
    {
      // Check that the rings are all the same
      const ring R = RingOf(A);
      if (RingOf(B) != R || RingOf(C) != R || RingOf(D) != R)
        CoCoA_ERROR(ERR::MixedRings, "BlockMat2x2(A,B,C,D)");
    }
    // Check that the matrix sizes are suitably matched.
    if (NumRows(A) != NumRows(B) ||
        NumRows(C) != NumRows(D) ||
        NumCols(A) != NumCols(C) ||
        NumCols(B) != NumCols(D))
      CoCoA_ERROR(ERR::BadMatrixSize, "BlockMat2x2(A,B,C,D)");
    if (NumCols(A) > numeric_limits<long>::max() -NumCols(B)) //avoid overflow  
      CoCoA_ERROR(ERR::BadColIndex, "BlockMat2x2(A,B,C,D)");
    if (NumRows(A) > numeric_limits<long>::max() -NumRows(C)) //avoid overflow  
      CoCoA_ERROR(ERR::BadRowIndex, "BlockMat2x2(A,B,C,D)");
    return MatrixView(new BlockMat2x2Impl(A, B, C, D));
  }


  /****************************************************************************/

  class ConstConcatVerImpl: public ConstMatrixViewBase
  {
  private:
    friend ConstMatrixView ConcatVer(ConstMatrixView A, ConstMatrixView B); // pseudo-ctor
    ConstConcatVerImpl(ConstMatrixView A, ConstMatrixView B);
    // default dtor works fine
  private: // disable default copy ctor and assignment
    ConstConcatVerImpl(const ConstConcatVerImpl&);            // NEVER DEFINED -- disable copy ctor
    ConstConcatVerImpl& operator=(const ConstConcatVerImpl&); // NEVER DEFINED -- disable assignment

  public:
    typedef std::vector<RingElem> vec;
    const ring& myRing() const;
    virtual long myNumRows() const;
    virtual long myNumCols() const;
    virtual RingElemAlias myEntry(long i, long j) const;
// Use default defns of the following two functions
//    virtual void myMulByRow(vec& lhs, const vec& v) const;
//    virtual void myMulByCol(vec& lhs, const vec& v) const;
    virtual bool myIsZeroRow(long i) const;
    virtual bool myIsZeroCol(long j) const;
// Use default defns of the following two functions
//    virtual void myDet(RingElem& d) const;
//    virtual long myRank() const;

  private: // data members
    const ring myR;
    long myNumRowsValue;
    long myNumColsValue;
    ConstMatrixView myA, myB;
  };



  ConstConcatVerImpl::ConstConcatVerImpl(ConstMatrixView A, ConstMatrixView B):
      myR(RingOf(A)),
      myNumRowsValue(NumRows(A)+NumRows(B)),
      myNumColsValue(NumCols(A)),
      myA(A),
      myB(B)
  {
    // Args are checked for compatibility by ConcatVer, the
    // only function which calls this ctor.
  }


  const ring& ConstConcatVerImpl::myRing() const
  {
    return myR;
  }


  long ConstConcatVerImpl::myNumRows() const
  {
    return myNumRowsValue;
  }


  long ConstConcatVerImpl::myNumCols() const
  {
    return myNumColsValue;
  }


  RingElemAlias ConstConcatVerImpl::myEntry(long i, long j) const
  {
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    if (i < NumRows(myA))
      return myA(i, j);
    return myB(i-NumRows(myA), j);
  }


  bool ConstConcatVerImpl::myIsZeroRow(long i) const
  {
    CoCoA_ASSERT(i < myNumRows());
    if (i < NumRows(myA))
      return myA->myIsZeroRow(i);
    return myB->myIsZeroRow(i-NumRows(myA));
  }


  bool ConstConcatVerImpl::myIsZeroCol(long j) const
  {
    CoCoA_ASSERT(j >= 0 && j < myNumCols());
    return myA->myIsZeroCol(j) && myB->myIsZeroCol(j);
  }


  ConstMatrixView ConcatVer(ConstMatrixView A, ConstMatrixView B)
  {
    if (RingOf(A) != RingOf(B))
      CoCoA_ERROR(ERR::MixedRings, "ConcatVer(A,B)");
    if (NumCols(A) != NumCols(B))
      CoCoA_ERROR(ERR::BadMatrixSize, "ConcatVer(A,B)");
    if (NumRows(A) > numeric_limits<long>::max() -NumRows(B)) //avoid overflow  
      CoCoA_ERROR(ERR::BadRowIndex, "ConcatVer(A,B)");

    return ConstMatrixView(new ConstConcatVerImpl(A, B));
  }


  /***************************************************************************/

  class ConcatVerImpl: public MatrixViewBase
  {
  private:
    friend MatrixView ConcatVer(MatrixView A, MatrixView B); // pseudo-ctor
    ConcatVerImpl(MatrixView A, MatrixView B);
    // default dtor works fine
  private: // disable default copy ctor and assignment
    ConcatVerImpl(const ConcatVerImpl&);            // NEVER DEFINED -- disable copy ctor
    ConcatVerImpl& operator=(const ConcatVerImpl&); // NEVER DEFINED -- disable assignment

  public:
    typedef std::vector<RingElem> vec;
    const ring& myRing() const;
    virtual long myNumRows() const;
    virtual long myNumCols() const;
    virtual RingElemAlias myEntry(long i, long j) const;
// Use default defns of the following two functions
//    virtual void myMulByRow(vec& lhs, const vec& v) const;
//    virtual void myMulByCol(vec& lhs, const vec& v) const;
    virtual bool myIsZeroRow(long i) const;
    virtual bool myIsZeroCol(long j) const;
// Use default defns of the following two functions
//    virtual void myDet(RingElem& d) const;
//    virtual long myRank() const;

    // Non-const member functions
    virtual void myAssignZero();
    virtual bool myIsWritable(long i, long j);
    virtual RingElemRawPtr myRawEntry(long i, long j);
    virtual void mySetEntry(long i, long j, ConstRefRingElem r);
    virtual void mySetEntry(long i, long j, const MachineInt& n);
    virtual void mySetEntry(long i, long j, const BigInt& N);
    virtual void mySetEntry(long i, long j, const BigRat& Q);

  private: // data members
    const ring myR;
    long myNumRowsValue;
    long myNumColsValue;
    MatrixView myA, myB;
  };


  /***************************************************************************/
  // Functions for ConcatVerImpl -- a matrix of the form  (A)
  //                                                      (B)
  // This implementation requires that the submatrices correspond
  // to two straight cuts -- several changes would be needed to handle
  // a more general assembly of submatrices.

  ConcatVerImpl::ConcatVerImpl(MatrixView A, MatrixView B):
      myR(RingOf(A)),
      myNumRowsValue(NumRows(A)+NumRows(B)),
      myNumColsValue(NumCols(A)),
      myA(A),
      myB(B)
  {
    // Args are checked for compatibility by ConcatVer, the
    // only function which calls this ctor.
  }


  const ring& ConcatVerImpl::myRing() const
  {
    return myR;
  }


  long ConcatVerImpl::myNumRows() const
  {
    return myNumRowsValue;
  }


  long ConcatVerImpl::myNumCols() const
  {
    return myNumColsValue;
  }


  RingElemAlias ConcatVerImpl::myEntry(long i, long j) const
  {
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    if (i < NumRows(myA))
      return myA(i, j);
    return myB(i-NumRows(myA), j);
  }


  bool ConcatVerImpl::myIsZeroRow(long i) const
  {
    CoCoA_ASSERT(i < myNumRows());
    if (i < NumRows(myA))
      return myA->myIsZeroRow(i);
    return myB->myIsZeroRow(i-NumRows(myA));
  }


  bool ConcatVerImpl::myIsZeroCol(long j) const
  {
    CoCoA_ASSERT(j >= 0 && j < myNumCols());
    return myA->myIsZeroCol(j) && myB->myIsZeroCol(j);
  }


  bool ConcatVerImpl::myIsWritable(long i, long j)
  {
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    if (i >= NumRows(myA))
      return myB->myIsWritable(i-NumRows(myA), j);
    else
      return myA->myIsWritable(i, j);
  }


  void ConcatVerImpl::myAssignZero()
  {
    myA->myAssignZero();
    myB->myAssignZero();
  }


  RingElemRawPtr ConcatVerImpl::myRawEntry(long i, long j)
  {
    CoCoA_ASSERT(myIsWritable(i,j));
    if (i >= NumRows(myA))
      return myB->myRawEntry(i-NumRows(myA), j);
    else
      return myA->myRawEntry(i, j);
  }


  void ConcatVerImpl::mySetEntry(long i, long j, ConstRefRingElem r)
  {
    CoCoA_ASSERT(myR == owner(r));
    CoCoA_ASSERT(myIsWritable(i,j));
    myR->myAssign(myRawEntry(i,j), raw(r)); //    myRefEntry(i,j) = r;
  }

  void ConcatVerImpl::mySetEntry(long i, long j, const MachineInt& n)
  {
    CoCoA_ASSERT(myIsWritable(i,j));
    myR->myAssign(myRawEntry(i,j), n); //    myRefEntry(i,j) = n;
  }

  void ConcatVerImpl::mySetEntry(long i, long j, const BigInt& N)
  {
    CoCoA_ASSERT(myIsWritable(i,j));
    myR->myAssign(myRawEntry(i,j), N); //    myRefEntry(i,j) = N;
  }


  void ConcatVerImpl::mySetEntry(long i, long j, const BigRat& Q)
  {
    CoCoA_ASSERT(myIsWritable(i,j));
    myR->myAssign(myRawEntry(i,j), Q); //    myRefEntry(i,j) = Q;
  }


  MatrixView ConcatVer(MatrixView A, MatrixView B)
  {
    {
      // Check that the base rings are the same.
      if (RingOf(A) != RingOf(B))
        CoCoA_ERROR(ERR::MixedRings, "ConcatVer(A,B)");
    }
    if (NumCols(A) != NumCols(B))
      CoCoA_ERROR(ERR::BadMatrixSize, "ConcatVer(A,B)");
    return MatrixView(new ConcatVerImpl(A, B));
  }


  /***************************************************************************/
  // ConcatHor implemented in terms of ConcatVer.
  // Very simple implementations... but may not be so efficient.

  ConstMatrixView ConcatHor(ConstMatrixView M1, ConstMatrixView M2)
  {
    if (RingOf(M1) != RingOf(M2))
      CoCoA_ERROR(ERR::MixedRings, "ConcatHor(A,B)");
    if (NumRows(M1) != NumRows(M2))
      CoCoA_ERROR(ERR::BadMatrixSize, "ConcatHor(A,B)");

    return transpose(ConcatVer(transpose(M1), transpose(M2)));
  }

  MatrixView ConcatHor(MatrixView M1, MatrixView M2)
  {
    if (RingOf(M1) != RingOf(M2))
      CoCoA_ERROR(ERR::MixedRings, "ConcatHor(A,B)");
    if (NumRows(M1) != NumRows(M2))
      CoCoA_ERROR(ERR::BadMatrixSize, "ConcatHor(A,B)");

    return transpose(ConcatVer(transpose(M1), transpose(M2)));
  }


  /***************************************************************************/

  MatrixView ConcatDiag(MatrixView M1, MatrixView M2)
  {
    ring R = RingOf(M1);
    if (RingOf(M2) != R)
      CoCoA_ERROR(ERR::MixedRings, "ConcatDiag(A,B)");

    return BlockMat2x2(M1,
                    ZeroMat(R, NumRows(M1), NumCols(M2)),
                    ZeroMat(R, NumRows(M2), NumCols(M1)),
                    M2);
  }

  ConstMatrixView ConcatDiag(ConstMatrixView M1, ConstMatrixView M2)
  {
    ring R = RingOf(M1);
    if (RingOf(M2) != R)
      CoCoA_ERROR(ERR::MixedRings, "ConcatDiag(A,B)");

    return BlockMat2x2(M1,
                    ZeroMat(R, NumRows(M1), NumCols(M2)),
                    ZeroMat(R, NumRows(M2), NumCols(M1)),
                    M2);
  }


  MatrixView ConcatAntiDiag(MatrixView M1, MatrixView M2)
  {
    ring R = RingOf(M1);
    if (RingOf(M2) != R)
      CoCoA_ERROR(ERR::MixedRings, "ConcatAntiDiag(A,B)");

    return BlockMat2x2(ZeroMat(R, NumRows(M1), NumCols(M2)),
                    M1,
                    M2,
                    ZeroMat(R, NumRows(M2), NumCols(M1)));
  }

  ConstMatrixView ConcatAntiDiag(ConstMatrixView M1, ConstMatrixView M2)
  {
    ring R = RingOf(M1);
    if (RingOf(M2) != R)
      CoCoA_ERROR(ERR::MixedRings, "ConcatAntiDiag(A,B)");

    return BlockMat2x2(ZeroMat(R, NumRows(M1), NumCols(M2)),
                    M1,
                    M2,
                    ZeroMat(R, NumRows(M2), NumCols(M1)));
  }



} // End of namespace CoCoA

// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/MatrixView.C,v 1.9 2014/07/31 14:45:17 abbott Exp $
// $Log: MatrixView.C,v $
// Revision 1.9  2014/07/31 14:45:17  abbott
// Summary: Merged io.H and UtilsTemplate.H into new header VectorOperations.H
// Author: JAA
//
// Revision 1.8  2014/07/30 14:07:13  abbott
// Summary: Changed BaseRing into RingOf; myBaseRing --> myRing
// Author: JAA
//
// Revision 1.7  2014/07/17 14:51:02  abbott
// Summary: Added MatByColsImpl; minor cleaning to MatByRows
// Author: JAA
//
// Revision 1.6  2014/07/17 13:19:06  abbott
// Summary: Corrected (lots of) silly typos -- it's too hot & brain is not functioning right
// Author: JAA
//
// Revision 1.5  2014/07/17 13:17:07  abbott
// Summary: Revised code that checks for duplicate indices in submat(...)
// Author: JAA
//
// Revision 1.4  2014/07/14 15:07:27  abbott
// Summary: Changed include of tmp.H into UtilsTemplate.H
// Author: JAA
//
// Revision 1.3  2014/06/17 10:10:00  abbott
// Summary: Added many (void)(var_name) to avoid compiler warning about unused param
// Author: JAA
//
// Revision 1.2  2014/04/30 16:09:19  abbott
// Summary: Replaced X.size() by len(X)
// Author: JAA
//
// Revision 1.1  2014/04/17 13:40:56  abbott
// Summary: renamed from MatrixViews.C
// Author: JAA
//
// Revision 1.25  2014/04/10 15:41:12  abbott
// Summary: Removed FilledMat (matrix view)
// Author: JAA
//
// Revision 1.24  2013/07/30 15:00:47  bigatti
// -- equality of ZeroMat now throws if mixed rings
//
// Revision 1.23  2013/05/31 14:52:03  abbott
// Added new views "MatByRows" and "MatByCols" (impl. is untested &  may be incomplete).
//
// Revision 1.22  2013/05/31 12:42:02  bigatti
// changed BlockMat into BlockMat2x2
// changed MultiBlockMat into BlockMat
//
// Revision 1.21  2012/10/24 12:14:43  abbott
// Changed return type of <matrix>::myEntry  (many times).
//
// Revision 1.20  2012/10/16 09:57:20  abbott
// Replaced  RefRingElem  by RingElem&  as param for <matrix>::myDet (many times)
//
// Revision 1.19  2012/10/11 14:57:20  abbott
// Replaced mem fn  myRefEntry  by equivalent new mem fn  myRawEntry;
// this way we remove the need for RefRingElem.
// Also the new name/semantics should discourage casual use.
//
// Revision 1.18  2011/12/23 14:55:56  bigatti
// -- changed name of parameter (x-->val)
//
// Revision 1.17  2011/11/09 14:09:53  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.16  2011/08/24 10:25:53  bigatti
// -- renamed QQ --> BigRat
// -- sorted #include
//
// Revision 1.15  2011/08/14 15:52:17  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.14  2011/03/21 13:21:28  bigatti
// -- added FilledMat
// -- added some CoCoA_ASSERT
//
// Revision 1.13  2011/03/10 18:10:12  bigatti
// -- using len instead of size() in CoCoA_ASSERT
//
// Revision 1.12  2011/03/09 17:20:28  bigatti
// -- size checks are now better organised
// -- using len(v) instead of v.size()
//
// Revision 1.11  2011/03/09 08:00:36  bigatti
// -- changed row/col args into long (instead of MachineInt)
// -- added #include "MachineInt.H"
//
// Revision 1.10  2011/03/04 16:35:10  bigatti
// -- changed: functions args of type MachineInt instead of size_t
//             members functions args of type long instead of size_t
//
// Revision 1.9  2011/03/03 13:50:22  abbott
// Replaced several occurrences of std::size_t by long; there's still more
// work to do though!
//
// Revision 1.8  2011/02/22 14:23:59  bigatti
// -- fixed IamAntiSymmetric
//
// Revision 1.7  2011/02/22 13:15:10  bigatti
// -- added IsSymmetric, IsAntiSymmetric, IsDiagonal, and tests
//
// Revision 1.6  2011/02/16 17:37:10  bigatti
// -- optimized IamEqual for diagonal MatrixViews
//
// Revision 1.5  2011/01/31 14:10:30  bigatti
// -- added mySetEntry/SetEntry with BigRat entry
//
// Revision 1.4  2009/12/03 17:40:36  abbott
// Added include directives for ZZ.H (as a consequence of removing
// the directive from ring.H).
//
// Revision 1.3  2009/06/25 16:59:42  abbott
// Minor improvement to some error messages (better coherence & comprehensibility).
//
// Revision 1.2  2008/12/17 12:11:52  abbott
// Changed type from long to MachineInt in operations which use a machine integer
// in place of a RingElem.  The change is "superficial" but affects many files.
//
// Revision 1.1  2008/04/18 15:35:57  abbott
// (long overdue) Major revision to matrices
//
//
