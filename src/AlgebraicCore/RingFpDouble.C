//   Copyright (c)  2005-2012,2014  John Abbott

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


#include "CoCoA/RingFpDouble.H"

#include "CoCoA/BigRat.H"
#include "CoCoA/FieldIdeal.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/IntOperations.H"
#include "CoCoA/MemPool.H"
#include "CoCoA/OpenMath.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/SmallFpDoubleImpl.H"
#include "CoCoA/assert.H"
#include "CoCoA/convert.H"
#include "CoCoA/error.H"
#include "CoCoA/ideal.H"


#include <algorithm>
//using std::swap;         // only in mySwap
#include <iostream>
using std::ostream;        // only in myOutput
// #include <limits>  ---  included in MachineInt.H (included via BigRat.H)
using std::numeric_limits; // only in ctor
// #include <memory>  ---  included in MemPool.H
// #include <vector>  ---  included in ideal.H
using std::vector;


namespace CoCoA
{

  class RingFpDoubleImpl: public QuotientRingBase
  {
  private: // data members
    typedef SmallFpDoubleImpl::value_t value_t;
    const unsigned long myModulus;
    const SmallFpDoubleImpl myImpl;
    mutable MemPool myMemMgr;       // MemPool must come *BEFORE* myZeroPtr and myOnePtr
    std::unique_ptr<RingElem> myZeroPtr;  ///< Every ring stores its own zero.
    std::unique_ptr<RingElem> myOnePtr;   ///< Every ring stores its own one.

  private: // auxiliary functions
    static unsigned long PrincipalGen(const ideal& I); // used for arg checking in ctor
    static value_t& import(RingElemRawPtr rawx);
    static const value_t& import(RingElemConstRawPtr rawx);

  private:
    RingFpDoubleImpl(const ideal& I, GlobalSettings::ResidueSetting ResidueChoice = DefaultResidueSetting()); // called only by NewRingFpDouble
    ~RingFpDoubleImpl();
    friend QuotientRing NewRingFpDouble(const MachineInt& p, GlobalSettings::ResidueSetting ResidueChoice);
    friend QuotientRing NewRingFpDouble(const BigInt& P);
    friend QuotientRing NewRingFpDouble(const ideal& I);
  private:
    RingFpDoubleImpl(const RingFpDoubleImpl&);            ///< NEVER DEFINED -- disallow copy construction
    RingFpDoubleImpl& operator=(const RingFpDoubleImpl&); ///< NEVER DEFINED -- disallow assignment
  public:

    // functions inherited by RingBase (ring.H)
    //-- ring
    virtual void myCharacteristic(BigInt& p) const;
    virtual long myLogCardinality() const;
    virtual bool IamCommutative() const;
    virtual bool3 IamIntegralDomain3(bool) const;
    virtual bool IamField() const;
    virtual bool IamFiniteField() const;
    virtual bool IamExact() const;
    virtual void myOutputSelf(std::ostream& out) const;
    virtual void myOutputSelf(OpenMathOutput& OMOut) const;
    //-- RingElem: memory
    virtual ConstRefRingElem myZero() const;
    virtual ConstRefRingElem myOne() const;
    using RingBase::myNew;    // disable warnings of overloading
    using RingBase::myAssign; // disable warnings of overloading
    virtual RingElemRawPtr myNew() const;
    virtual RingElemRawPtr myNew(const MachineInt& n) const;
    virtual RingElemRawPtr myNew(const BigInt& N) const;
    virtual RingElemRawPtr myNew(ConstRawPtr rawx) const;
    virtual void myDelete(RawPtr rawx) const;
    //-- RingElem: assignments
    virtual void myAssign(RawPtr rawlhs, ConstRawPtr rawx) const;
    virtual void myAssign(RawPtr rawlhs, const MachineInt& n) const;
    virtual void myAssign(RawPtr rawlhs, const BigInt& N) const;
    virtual void myAssignZero(RawPtr rawlhs) const;
    virtual void myRecvTwinFloat(RawPtr rawlhs, ConstRawPtr rawx) const;
    virtual void mySwap(RawPtr rawx, RawPtr rawy) const;
    //-- RingElem: arithmetic
    virtual void myNegate(RawPtr rawlhs, ConstRawPtr rawx) const;
    virtual void myAdd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const;
    virtual void mySub(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const;
    virtual void myMul(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const;
    virtual void myDiv(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const;
    virtual void myGcd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const; // lhs = gcd(x,y) in a field
    //-- RingElem: queries
    virtual bool myIsZero(ConstRawPtr rawx) const;
    virtual bool myIsOne(ConstRawPtr rawx) const;
    virtual bool myIsMinusOne(ConstRawPtr rawx) const;
    virtual bool myIsEqual(ConstRawPtr rawx, ConstRawPtr rawy) const;
    virtual bool myIsInteger(BigInt& N, ConstRawPtr rawx) const; ///< always true
    virtual bool myIsRational(BigRat& Q, ConstRawPtr rawx) const;
    virtual bool myIsDouble(double& d, ConstRawPtr rawx) const;
    virtual bool myIsDivisible(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const;
    virtual bool myIsInvertible(ConstRawPtr rawx) const;
    virtual bool myIsZeroAddMul(RawPtr rawlhs, ConstRawPtr rawy, ConstRawPtr rawz) const;
    virtual bool myIsZeroAddMul(RawPtr rawlhs, RawPtr /*rawtmp*/, ConstRawPtr rawy, ConstRawPtr rawz) const;
    //-- RingElem: printing
    virtual void myOutput(std::ostream& out, ConstRawPtr rawx) const;
    virtual void myOutput(OpenMathOutput& OMOut, ConstRawPtr rawx) const;
    virtual bool myIsPrintAtom(ConstRawPtr rawx) const;
    virtual bool myIsPrintedWithMinus(ConstRawPtr rawx) const;
    //-- ideal
    virtual ideal myIdealCtor(const std::vector<RingElem>& gens) const;
    //-- RingHom
    virtual RingHom myCompose(const RingHom& phi, const RingHom& theta) const; // phi(theta(...))

    virtual bool myImageLiesInSubfield(const RingHom& phi) const;

    //-- protected in RingBase
    virtual void myPowerSmallExp(RawPtr rawlhs, ConstRawPtr rawx, long n) const;
    virtual void myPowerBigExp(RawPtr rawlhs, ConstRawPtr rawx, const BigInt& N) const; // lhs = x^N, N big, x not -1,0,1

    // functions inherited by QuotientRingBase (QuotientRing.H)
    virtual RingElem myCanonicalRepr(ConstRawPtr rawx) const;
    virtual void myReduction(RawPtr rawimage, ConstRawPtr rawarg) const;
    virtual RingHom myInducedHomCtor(const RingHom& InducingHom) const;


  private:
    class InducedHomImpl: public RingHomBase
    {
    public:
      InducedHomImpl(const QuotientRing& domain, const RingHom& InducingHom);
      virtual void myApply(RingElemRawPtr rawimage, RingElemConstRawPtr rawarg) const;
      virtual bool IamPartial() const { return false; }
    };

  };



  inline RingFpDoubleImpl::value_t& RingFpDoubleImpl::import(RingElemRawPtr rawx)
  {
    return *static_cast<value_t*>(rawx.myRawPtr());
  }

  inline const RingFpDoubleImpl::value_t& RingFpDoubleImpl::import(RingElemConstRawPtr rawx)
  {
    return *static_cast<const value_t*>(rawx.myRawPtr());
  }


  // Returns generator of I as a value_t; returns 0 if value is too large to fit.
  unsigned long RingFpDoubleImpl::PrincipalGen(const ideal& I)
  {
    if (IsZero(I)) return 0;
    const BigInt GenI = ConvertTo<BigInt>(TidyGens(I)[0]);
    unsigned long p;
    if (!IsConvertible(p, GenI))  // check that the value of the principal generator will fit
      return 0;
    return p;
  }


  RingFpDoubleImpl::RingFpDoubleImpl(const ideal& I, GlobalSettings::ResidueSetting ResidueChoice):
      QuotientRingBase(RingZZ(), I),  // confirms that I is an ideal of Z
      myModulus(PrincipalGen(I)),
      myImpl(myModulus, ResidueChoice),   // also checks that myModulus is a small prime
      myMemMgr(SmallFpDoubleImpl::ourDatumSize, "RingFpDoubleImpl.myMemMgr")
  {
    myRefCountInc();  // this is needed for exception cleanliness, in case one of the lines below throws
    myZeroPtr.reset(new RingElem(ring(this)));
    myOnePtr.reset(new RingElem(ring(this), 1));
    myRefCountZero();
  }


  RingFpDoubleImpl::~RingFpDoubleImpl()
  {}


  void RingFpDoubleImpl::myCharacteristic(BigInt& p) const
  {
    p = myModulus;
  }


  long RingFpDoubleImpl::myLogCardinality() const
  {
    return 1;
  }


  bool RingFpDoubleImpl::IamCommutative() const
  {
    return true;
  }


  bool3 RingFpDoubleImpl::IamIntegralDomain3(bool) const
  {
    return true3;
  }


  bool RingFpDoubleImpl::IamField() const
  {
    return true;
  }


  bool RingFpDoubleImpl::IamFiniteField() const
  {
    return true;
  }


  bool RingFpDoubleImpl::IamExact() const
  {
    return true;
  }


  ConstRefRingElem RingFpDoubleImpl::myZero() const
  {
    return *myZeroPtr;
  }


  ConstRefRingElem RingFpDoubleImpl::myOne() const
  {
    return *myOnePtr;
  }


  RingElemRawPtr RingFpDoubleImpl::myNew() const
  {
    value_t* ans = static_cast<value_t*>(myMemMgr.alloc());
    *ans = 0;
    return RingElemRawPtr(ans);
  }


  RingElemRawPtr RingFpDoubleImpl::myNew(const MachineInt& n) const
  {
    value_t* ans = static_cast<value_t*>(myMemMgr.alloc());
    *ans = myImpl.myReduce(n);
    return RingElemRawPtr(ans);
  }


  RingElemRawPtr RingFpDoubleImpl::myNew(const BigInt& N) const
  {
    value_t* ans = static_cast<value_t*>(myMemMgr.alloc());
    *ans = myImpl.myReduce(N);
    return RingElemRawPtr(ans);
  }


  RingElemRawPtr RingFpDoubleImpl::myNew(ConstRawPtr rawy) const
  {
    value_t* ans = static_cast<value_t*>(myMemMgr.alloc());
    *ans  = import(rawy);
    return RingElemRawPtr(ans);
  }


  void RingFpDoubleImpl::myDelete(RawPtr rawx) const
  {
    myMemMgr.free(rawx.myRawPtr());
  }


  void RingFpDoubleImpl::mySwap(RawPtr rawx, RawPtr rawy) const
  {
    std::swap(import(rawx), import(rawy));
  }


  void RingFpDoubleImpl::myAssign(RawPtr rawlhs, ConstRawPtr rawx) const
  {
    import(rawlhs) = import(rawx);
  }


  void RingFpDoubleImpl::myAssign(RawPtr rawlhs, const MachineInt& n) const
  {
    import(rawlhs) = myImpl.myReduce(n);
  }


  void RingFpDoubleImpl::myAssign(RawPtr rawlhs, const BigInt& N) const
  {
    import(rawlhs) = myImpl.myReduce(N);
  }


  void RingFpDoubleImpl::myAssignZero(RawPtr rawlhs) const
  {
    import(rawlhs) = 0;
  }


  void RingFpDoubleImpl::myRecvTwinFloat(RawPtr /*rawlhs*/, ConstRawPtr /*rawx*/) const
  {
    CoCoA_ERROR(ERR::SERIOUS, "RingFpDoubleImpl::myRecvTwinFloat");
  }


  void RingFpDoubleImpl::myNegate(RawPtr rawlhs, ConstRawPtr rawx) const
  {
    import(rawlhs) = myImpl.myNegate(import(rawx));
  }


  void RingFpDoubleImpl::myAdd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    import(rawlhs) = myImpl.myAdd(import(rawx), import(rawy));
  }


  void RingFpDoubleImpl::mySub(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    import(rawlhs) = myImpl.mySub(import(rawx), import(rawy));
  }


  void RingFpDoubleImpl::myMul(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    import(rawlhs) = myImpl.myMul(import(rawx), import(rawy));
  }


  void RingFpDoubleImpl::myDiv(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    CoCoA_ASSERT(!myIsZero(rawy));
    import(rawlhs) = myImpl.myDiv(import(rawx), import(rawy));
  }


  bool RingFpDoubleImpl::myIsDivisible(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    if (myIsZero(rawy)) return false;
    myDiv(rawlhs, rawx, rawy);
    return true;
  }


  bool RingFpDoubleImpl::myIsInvertible(ConstRawPtr rawx) const
  {
    return !myIsZero(rawx);
  }


  void RingFpDoubleImpl::myGcd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    myGcdInField(rawlhs, rawx, rawy);
  }


  void RingFpDoubleImpl::myPowerSmallExp(RawPtr rawlhs, ConstRawPtr rawx, long n) const  // assumes n > 1
  {
    // Assert that we have a genuinely non-trivial case.
    CoCoA_ASSERT(n > 1);
    CoCoA_ASSERT(!myIsZero(rawx) && !myIsOne(rawx) && !myIsMinusOne(rawx));
    import(rawlhs) = myImpl.myPower(import(rawx), n);
  }

  void RingFpDoubleImpl::myPowerBigExp(RawPtr rawlhs, ConstRawPtr rawx, const BigInt& N) const
  {
    // Assert that we have a genuinely non-trivial case.
    CoCoA_ASSERT(N > 1);
    CoCoA_ASSERT(!myIsZero(rawx) && !myIsOne(rawx) && !myIsMinusOne(rawx));
    // Use Fermat's Little Theorem to reduce exponent...
    import(rawlhs) = myImpl.myPower(import(rawx), N%(myModulus-1));
  }


  void RingFpDoubleImpl::myOutput(ostream& out, ConstRawPtr rawx) const
  {
    out << myImpl.myExport(import(rawx));
  }


  bool RingFpDoubleImpl::myIsPrintAtom(ConstRawPtr rawx) const
  {
    return myImpl.myExport(import(rawx)) >= 0;
  }


  bool RingFpDoubleImpl::myIsPrintedWithMinus(ConstRawPtr rawx) const
  {
    return myImpl.myExport(import(rawx)) < 0;
  }


  void RingFpDoubleImpl::myOutputSelf(ostream& out) const
  {
    //    out << "FFp(" << myModulus << ")";
    out << "RingWithID(" << myID << ",\"FFp(" << myModulus << ")\")";
  }


  void RingFpDoubleImpl::myOutputSelf(OpenMathOutput& OMOut) const
  {
    OMOut->mySendApplyStart();
    OMOut << OpenMathSymbol("setname2", "GFp");
    OMOut << myModulus;
    OMOut->mySendApplyEnd();
  }


  void RingFpDoubleImpl::myOutput(OpenMathOutput& OMOut, ConstRawPtr rawx) const
  {
    OMOut << myImpl.myExport(import(rawx));
  }


  bool RingFpDoubleImpl::myIsZero(ConstRawPtr rawx) const
  {
    return (import(rawx) == 0);
  }


  bool RingFpDoubleImpl::myIsOne(ConstRawPtr rawx) const
  {
    return (import(rawx) == 1);
  }


  bool RingFpDoubleImpl::myIsMinusOne(ConstRawPtr rawx) const
  {
    return (import(rawx) == myImpl.myReduce(-1));
  }


  bool RingFpDoubleImpl::myIsInteger(BigInt& N, ConstRawPtr rawx) const
  {
    N = myImpl.myExport(import(rawx));
    return true;
  }


  bool RingFpDoubleImpl::myIsRational(BigRat& Q, ConstRawPtr rawx) const
  {
    Q = myImpl.myExport(import(rawx));
    return true;
  }


  bool RingFpDoubleImpl::myIsDouble(double& d, ConstRawPtr rawx) const
  {
    d = myImpl.myExport(import(rawx));
    return true;
  }


  bool RingFpDoubleImpl::myIsZeroAddMul(RawPtr rawlhs, ConstRawPtr rawfact1, ConstRawPtr rawfact2) const
  {
    return myImpl.myIsZeroAddMul(import(rawlhs), import(rawfact1), import(rawfact2));
  }


  bool RingFpDoubleImpl::myIsZeroAddMul(RawPtr rawlhs, RawPtr /*rawtmp*/, ConstRawPtr rawfact1, ConstRawPtr rawfact2) const
  {
    return myImpl.myIsZeroAddMul(import(rawlhs), import(rawfact1), import(rawfact2));
  }


  bool RingFpDoubleImpl::myIsEqual(ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    return (import(rawx) == import(rawy));
  }




  ideal RingFpDoubleImpl::myIdealCtor(const std::vector<RingElem>& gens) const
  {
    return NewFieldIdeal(ring(this), gens);
  }


  RingHom RingFpDoubleImpl::myCompose(const RingHom& phi, const RingHom& theta) const
  {
    // No need to check compatibility -- it was checked when theta and phi were built
    return RingHom(new InducedHomImpl(QuotientRing(this), phi(theta(myQuotientingHomCtor()))));
  }


  bool RingFpDoubleImpl::myImageLiesInSubfield(const RingHom& phi) const
  {
    (void)(phi); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(codomain(phi) == ring(this));
    return true;
  }


  RingElem RingFpDoubleImpl::myCanonicalRepr(ConstRawPtr rawx) const
  {
    return RingElem(myReprRing, myImpl.myExport(import(rawx)));
  }


  void RingFpDoubleImpl::myReduction(RawPtr rawimage, ConstRawPtr rawarg) const
  {
    BigInt tmp;
    CoCoA_ASSERT(myReprRing->myIsInteger(tmp, rawarg));
    myReprRing->myIsInteger(tmp, rawarg);
    import(rawimage) = myImpl.myReduce(tmp);
  }


  RingHom RingFpDoubleImpl::myInducedHomCtor(const RingHom& InducingHom) const
  {
    // Compatibility has already been checked (see InducedHom in QuotientRing.C)
    CoCoA_ASSERT(IsZero(InducingHom(myModulus)));
    return RingHom(new InducedHomImpl(QuotientRing(this), InducingHom));
  }


  //---------------------------------------------------------------------------
  // Functions to do with ring homomorphisms


  RingFpDoubleImpl::InducedHomImpl::InducedHomImpl(const QuotientRing& domain, const RingHom& InducingHom):
      RingHomBase(domain, codomain(InducingHom))
  { /* Compatibility already checked in InducedHom in QuotientRing.C */  }


  void RingFpDoubleImpl::InducedHomImpl::myApply(RingElemRawPtr rawimage, RingElemConstRawPtr rawarg) const
  {
    BigInt tmp;  //??? wasteful new/delete
    CoCoA_ASSERT(myDomain->myIsInteger(tmp, rawarg));
    myDomain->myIsInteger(tmp, rawarg);  // must necessarily succeed
    myCodomain->myAssign(rawimage, tmp);
  }



  QuotientRing NewRingFpDouble(const MachineInt& p, GlobalSettings::ResidueSetting ResidueChoice)
  {
    return QuotientRing(new RingFpDoubleImpl(ideal(RingElem(RingZZ(), p)), ResidueChoice));
  }

  QuotientRing NewRingFpDouble(const BigInt& P)
  {
    return QuotientRing(new RingFpDoubleImpl(ideal(RingElem(RingZZ(), P))));
  }

  QuotientRing NewRingFpDouble(const ideal& I)
  {
    if (!IsZZ(RingOf(I))) CoCoA_ERROR(ERR::IdealNotInRing, "NewRingFpDouble(I)");
    return QuotientRing(new RingFpDoubleImpl(I));
  }


  bool IsGoodForRingFpDouble(const MachineInt& p)
  {
    if (IsNegative(p) || !IsSignedLong(p)) return false;
    const long n = AsSignedLong(p);
    return SmallFpDoubleImpl::IsGoodCtorArg(n);
  }

  bool IsGoodForRingFpDouble(const BigInt& P)
  {
    if (P <= 0) return false;
    long p;
    if (!IsConvertible(p, P)) return false;
    return IsGoodForRingFpDouble(p);
  }

  bool IsGoodForRingFpDouble(const ideal& I)
  {
    if (!IsZZ(RingOf(I))) return false;
    if (IsZero(I)) return false;
    return IsGoodForRingFpDouble(ConvertTo<BigInt>(TidyGens(I)[0]));
  }


} // end of namespace CoCoA


// RCS header/log
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/RingFpDouble.C,v 1.41 2014/07/30 14:08:57 abbott Exp $
// $Log: RingFpDouble.C,v $
// Revision 1.41  2014/07/30 14:08:57  abbott
// Summary: Changed name AmbientRing --> RingOf
// Author: JAA
//
// Revision 1.40  2014/07/28 16:04:49  abbott
// Summary: Renamed myQuotientingHom to myQuotientingHomCtor
// Author: JAA
//
// Revision 1.39  2014/07/28 15:49:35  abbott
// Summary: Redesign: ringhoms no longer cached in rings (caused ref count trouble)
// Author: JAA
//
// Revision 1.38  2014/07/04 13:08:08  bigatti
// -- RingID into RingWithID
//
// Revision 1.37  2014/07/02 16:35:04  bigatti
// -- new way of printing ring with ID
//
// Revision 1.36  2014/06/17 10:13:10  abbott
// Summary: Added (void)(phi) to avoid compiler warning about unused param
// Author: JAA
//
// Revision 1.35  2014/05/14 15:57:15  bigatti
// -- added "using" for clang with superpedantic flag
//
// Revision 1.34  2014/04/02 10:57:46  abbott
// Summary: Revised design of IamIntegralDomain3
// Author: JAA
//
// Revision 1.33  2014/03/27 17:17:31  abbott
// Summary: Added new fn IsIntegralDomain3 (and mem fn IamIntegralDomain3)
// Author: JAA
//
// Revision 1.32  2014/01/28 10:02:48  abbott
// Replaced some calls to IsInteger by calls to ConvertTo<BigInt>.
//
// Revision 1.31  2013/03/25 17:04:19  abbott
// Major clean-up of interface to SmallFpImpl/SmallFpLogImpl/SmallFpDoubleImpl
// (underlying impl remains much the same).  Removed lots of cruft.
// Consequential changes to RingFp* classes; small change to SparsePolyRing.
//
// Revision 1.30  2012/09/07 15:21:13  abbott
// First stage of revision of SmallFpImpl interface (and SmallFpLog, SmallFpDouble).
//
// Revision 1.29  2012/05/28 09:18:21  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.28  2012/04/27 15:04:10  abbott
// Added mem fns IamFiniteField & myLogCardinality.
//
// Revision 1.27  2012/03/16 14:43:48  bigatti
// -- re-sorted
//
// Revision 1.26  2012/02/10 10:28:08  bigatti
// -- changed RingZ.H, RingQ.H --> RingZZ.H, RingQQ.H
//
// Revision 1.25  2012/02/08 17:06:40  bigatti
// -- changed: Z,Q -> ZZ,QQ
// -- reordered member functions (to be completed)
//
// Revision 1.24  2012/01/30 23:20:48  abbott
// Realigned comments (I hope, CVS did something funny though).
//
// Revision 1.23  2012/01/25 13:33:51  bigatti
// -- added myIsZeroAddMul with 4 args
// -- some tidying
//
// Revision 1.22  2011/11/09 14:11:58  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.21  2011/08/24 10:28:49  bigatti
// -- renamed QQ --> BigRat
// -- sorted #include
//
// Revision 1.20  2011/08/14 15:52:16  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.19  2011/06/23 16:04:47  abbott
// Added IamExact mem fn for rings.
// Added myRecvTwinFloat mem fn for rings.
// Added first imple of RingHom from RingTwinFloat to other rings.
//
// Revision 1.18  2011/05/24 14:54:29  abbott
// Consquential changes from removing several ctors for principal ideals.
//
// Revision 1.17  2011/05/20 19:26:05  abbott
// Updated SmallFp*Impl: removed all output-related fns (must use myExport instead).
//
// Revision 1.16  2011/05/20 09:45:04  abbott
// Harmonized RingFp, RingFpLog, RingFpDouble -- they are now almost ready to be merged into a single class!
//
// Revision 1.15  2011/05/19 14:38:27  abbott
// Updated small prime finite field impls to allow user to specify
// separately for each whether to use symmetric or non-negative
// residues for export operations (myExport and printing).
//
// Revision 1.14  2011/03/22 20:00:37  abbott
// Added IsGoodForXXX fns to test whether a given arg is suitable as
// characteristic for the given type of small prime finite field.
//
// Revision 1.13  2011/03/10 16:39:34  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.12  2010/10/06 14:10:24  abbott
// Added increments to the ref count in ring and PPMonoid ctors to make
// them exception safe.
//
// Revision 1.11  2009/12/11 11:46:32  abbott
// Changed fn  convert  into  IsConvertible.
// Added template procedure  convert.
// New version because change is not backward compatible.
//
// Revision 1.10  2009/09/25 14:01:11  abbott
// Cleaned up include directives.
//
// Revision 1.9  2009/09/24 16:22:29  abbott
// Tidied up include directives.  Removed some unnecessary "std::" prefixes.
//
// Revision 1.8  2009/07/02 16:32:11  abbott
// Consequential changes stemming from new class QQ, and modified interface to the member
// function RingBase::myIsRational.  Also some new conversion functions.
//
// Revision 1.7  2009/06/05 12:08:28  abbott
// Changed return type of operator%(ZZ,MachineInteger); it is now unsigned long
// instead of ZZ.
//
// Revision 1.6  2008/12/17 12:11:52  abbott
// Changed type from long to MachineInteger in operations which use a machine integer
// in place of a RingElem.  The change is "superficial" but affects many files.
//
// Revision 1.5  2007/10/30 17:14:07  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.4  2007/05/22 22:45:14  abbott
// Changed fn name IsUnit to IsInvertible.
//
// Revision 1.3  2007/05/21 12:44:25  abbott
// Changed impl of powering routine as a consequence of changed signature
// to the modulus operator for ZZ modulo machine integer.
//
// Revision 1.2  2007/03/23 18:38:42  abbott
// Separated the "convert" functions (and CheckedCast) into their own files.
// Many consequential changes.  Also corrected conversion to doubles.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.14  2007/03/08 18:22:29  cocoa
// Just whitespace cleaning.
//
// Revision 1.13  2007/03/08 10:23:29  bigatti
// -- CanonHom --> CanonicalHom
//
// Revision 1.12  2007/03/05 21:06:07  cocoa
// New names for homomorphism pseudo-ctors: removed the "New" prefix.
//
// Revision 1.11  2007/03/03 14:07:23  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.10  2007/03/02 10:47:53  cocoa
// First stage of RingZ modifications -- tests do not compile currently, Anna will fix this.
//
// Revision 1.9  2007/01/17 12:32:39  cocoa
// Changed a few more "raw" variable names so that the code compiles fine
// also when CoCoA_DEBUG is set.
//
// Revision 1.8  2007/01/15 16:00:10  cocoa
// -- added prefix "raw" to RawPtr arguments names
// -- changed rhs into rawx, n, or N
//
// Revision 1.7  2007/01/13 14:14:34  cocoa
// Overhaul of RingHom code: it nows uses SmartPtrIRC, and printing is more logical.
// Have not yet updated the documentation.
//
// Revision 1.6  2006/11/08 16:21:59  cocoa
// Structural cleaning of RingHom; many consequential changes.
//
// Revision 1.5  2006/11/02 13:25:44  cocoa
// Simplification of header files: the OpenMath classes have been renamed.
// Many minor consequential changes.
//
// Revision 1.4  2006/10/16 23:18:59  cocoa
// Corrected use of std::swap and various special swap functions.
// Improved myApply memfn for homs of RingDistrMPolyInlPP.
//
// Revision 1.3  2006/10/06 14:04:14  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.2  2006/10/06 10:15:52  cocoa
// In response to Susan's bug: a fiasco when compiling with CoCoA_MEMPOOL_DEBUG
// set wrongly.  Moved several implementation classes out of their header files
// into the implementation files.  Several functions had to be uninlined.
// Also corrected position of #include, etc.
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.7  2006/05/29 16:22:37  cocoa
// Third time lucky???
// Added myIsInteger member function to all rings (NYI for RingFloat).
//
// Revision 1.6  2006/05/12 16:10:58  cocoa
// Added OpenMathFwd.H, and tidied OpenMath.H.
// Many consequential but trivial changes.
//
// Revision 1.5  2006/04/21 15:01:36  cocoa
// Changed default implementation of RingBase::myGcd -- it now gives a SERIOUS
// error.  All fields must now handle a call to gcd explicitly: they can use
// the new myGcdInField function.  It's now cleaner than it was.
//
// Revision 1.4  2006/03/15 18:09:31  cocoa
// Changed names of member functions which print out their object
// into myOutputSelf -- hope this will appease the Intel C++ compiler.
//
// Revision 1.3  2006/03/14 15:01:49  cocoa
// Improved the implementation of ring member fns for computing powers.
// Should keep Intel C++ compiler quieter too.
//
// Revision 1.2  2006/03/12 21:28:33  cocoa
// Major check in after many changes
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.3  2005/10/14 15:25:07  cocoa
// Major tidying and cleaning to small prime finite fields.
// Several consequential changes.  Improved their documentation.
//
// Added Makefile and script to include/CoCoA/ directory to
// keep library.H up to date.
//
// Revision 1.2  2005/10/12 15:52:09  cocoa
// Completed test-RingFp1 and corrected/cleaned the SmallFp*
// and RingFp* files.
//
// Some minor tidying elsewhere.
//
// Revision 1.1  2005/10/11 16:37:30  cocoa
// Added new small prime finite field class (see RingFpDouble).
//
// Cleaned makefiles and configuration script.
//
// Tidied PPMonoid code (to eliminate compiler warnings).
//
// Fixed bug in RingFloat::myIsInteger.
//
