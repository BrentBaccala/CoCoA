//   Copyright (c)  2003-2009,2011,2014  John Abbott

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


// Source code for class GeneralQuotientRingImpl

#include "CoCoA/QuotientRing.H"

#include "CoCoA/GlobalManager.H"
#include "CoCoA/BigInt.H"
#include "CoCoA/IntOperations.H"
#include "CoCoA/NumTheory.H"
#include "CoCoA/OpenMath.H"
#include "CoCoA/PolyRing.H"
#include "CoCoA/RingFp.H"
#include "CoCoA/RingFpDouble.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/assert.H"
#include "CoCoA/convert.H"
#include "CoCoA/error.H"
#include "CoCoA/utils.H"


#include <iostream>
using std::ostream;
#include <memory>
using std::auto_ptr;
#include <vector>
using std::vector;


namespace CoCoA
{

  const QuotientRingBase* QuotientRingPtr(const ring& R, const char* const FnName)
  {
    const QuotientRingBase* ptr = QuotientRingPtr(R);
    if (ptr == 0/*nullptr*/) CoCoA_ERROR(ERR::NotQuotientRing, FnName);
    return ptr;
  }


  class QuotientingHomImpl: public RingHomBase
  {
  public:
    QuotientingHomImpl(const QuotientRing& RmodI);
    virtual void myApply(RingElemRawPtr rawimage, RingElemConstRawPtr rawarg) const;
    virtual bool IamPartial() const { return false; }
  private:
    virtual void myOutputSelfDetails(std::ostream& out) const;
  };


  RingHom QuotientRingBase::myQuotientingHomCtor() const
  {
    return RingHom(new QuotientingHomImpl(QuotientRing(this)));
  }


  RingElem QuotientRingBase::mySymbolValue(const symbol& sym) const
  {
    return myQuotientingHomCtor()(myReprRing->mySymbolValue(sym));
  }


  class GeneralQuotientRingImpl: public QuotientRingBase
  {
  private: // data members
    auto_ptr<RingElemAlias> myZeroPtr;
    auto_ptr<RingElemAlias> myOnePtr;

  public:
    GeneralQuotientRingImpl(const ring& R, const ideal& I);
    // Default copy ctor works fine.
    // Assignment disabled
    ~GeneralQuotientRingImpl() {};
  private:
    GeneralQuotientRingImpl& operator=(const GeneralQuotientRingImpl&); // NEVER DEFINED
  public: // functions that every ring must have
    virtual void myCharacteristic(BigInt& p) const;
    virtual long myLogCardinality() const;
    virtual bool IamCommutative() const {return IsCommutative(myBaseRingValue);} // BUG??? this could fail to recognize commutativity
    virtual bool3 IamIntegralDomain3(bool QuickMode) const { if (QuickMode) return uncertain3; return bool3(IsPrime(myReducingIdeal)); /*return IsPrime3(myReducingIdeal);*/ }
    virtual bool IamTrueGCDDomain() const {return false;} // we don't know how to compute GCDs
    virtual bool IamField() const {return IsMaximal(myReducingIdeal);}
    virtual bool IamFiniteField() const;
    virtual bool IamExact() const {return IsExact(myReprRing);}
    virtual ConstRefRingElem myZero() const {return *myZeroPtr;}
    virtual ConstRefRingElem myOne() const  {return *myOnePtr;}
    using RingBase::myNew;    // disable warnings of overloading
    using RingBase::myAssign; // disable warnings of overloading
    virtual RingElemRawPtr myNew() const;
    virtual RingElemRawPtr myNew(const MachineInt& n) const;
    virtual RingElemRawPtr myNew(const BigInt& N) const;
    virtual RingElemRawPtr myNew(ConstRawPtr rawt) const;
    virtual void myDelete(RawPtr rawx) const;                      // destroys x (incl all resources)
    virtual void mySwap(RawPtr rawx, RawPtr rawy) const;                           // swap(x, y)
    virtual void myAssign(RawPtr rawlhs, ConstRawPtr rawx) const;                  // lhs = x
    virtual void myAssign(RawPtr rawlhs, const MachineInt& n) const;           // lhs = n
    virtual void myAssign(RawPtr rawlhs, const BigInt& N) const;                   // lhs = N
    virtual void myAssignZero(RawPtr rawlhs) const;                                // lhs = 0
    virtual void myRecvTwinFloat(RawPtr rawlhs, ConstRawPtr rawx) const;
    virtual void myNegate(RawPtr rawlhs, ConstRawPtr rawx) const;                  // lhs = -x
    virtual void myAdd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const;   // lhs = x+y
    virtual void mySub(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const;   // lhs = x-y
    virtual void myMul(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const;   // lhs = x*y
    virtual void myDiv(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const;   // lhs = x/y
    virtual bool myIsDivisible(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const;// lhs = x/y, if divisible
    virtual bool myIsInvertible(ConstRawPtr rawx) const;                           // true iff x is invertible
    virtual void myGcd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const;   // lhs = gcd(x,y) if TrueGCDDomain; o/w error.
    virtual void myPowerSmallExp(RawPtr rawlhs, ConstRawPtr rawx, long n) const;   // lhs = x^n, n>1, x not -1,0,1
    virtual void mySymbols(std::vector<symbol>& SymList) const {myReprRing->mySymbols(SymList);}
    virtual void myOutput(std::ostream& out, ConstRawPtr rawx) const;              // out << x
    virtual bool myIsPrintAtom(ConstRawPtr rawx) const;                            ///< x^n may be printed without parentheses
    virtual bool myIsPrintedWithMinus(ConstRawPtr rawx) const;                     ///< first character of x printed is a minus sign
    virtual void myOutputSelf(std::ostream& out) const;                            // out << R
    virtual void myOutputSelfLong(std::ostream& out) const; // out << R (descr)
    virtual void myOutputSelf(OpenMathOutput& OMOut) const;                        // OMOut << R
    virtual void myOutput(OpenMathOutput& OMOut, ConstRawPtr rawx) const;          // OMOut << x
    virtual bool myIsZero(ConstRawPtr rawx) const;                                 // x == 0
    virtual bool myIsOne(ConstRawPtr rawx) const;                                  // x == 1
    virtual bool myIsMinusOne(ConstRawPtr rawx) const;                             // x == -1
    virtual bool myIsInteger(BigInt& N, ConstRawPtr rawx) const;                   // true iff x is integer
    virtual bool myIsRational(BigRat& Q, ConstRawPtr rawx) const;                      // true iff x is rational
    virtual bool myIsDouble(double& d, ConstRawPtr rawx) const;                    // false iff x overflows
    virtual bool myIsEqual(ConstRawPtr rawx, ConstRawPtr rawy) const;              // x == y
  //???  virtual void convert(string&, ConstRawPtr) const;

    virtual ideal myIdealCtor(const std::vector<RingElem>& gens) const;

    virtual RingHom myCompose(const RingHom& phi, const RingHom& theta) const; // phi(theta(...))
    virtual bool myImageLiesInSubfield(const RingHom& phi) const;

    virtual void myPowerBigExp(RawPtr rawlhs, ConstRawPtr rawx, const BigInt& N) const;  ///< lhs = x^N (non-triv, N large); default gives error

  public: // functions that every QuotientRing must have
    virtual RingElem myCanonicalRepr(ConstRawPtr rawx) const;
    virtual void myReduction(RawPtr rawimage, ConstRawPtr rawarg) const;
    virtual RingHom myInducedHomCtor(const RingHom& InducingHom) const; // NB domain(InducingHom) == myReprRing, and ker(InducingHom) contains myReducingIdeal.

  private:
    void myReduction(RawPtr rawx) const; // implementation detail


  private:
    class InducedHomImpl: public RingHomInducedBase
    {
    public:
      InducedHomImpl(const QuotientRing& RmodI, const RingHom& InducingHom);
      // Default copy ctor & assignment disabled in RingHomBase.
      // Default dtor works fine.
      virtual void myApply(RingElemRawPtr rawimage, RingElemConstRawPtr rawarg) const;
      virtual bool IamPartial() const { return IsPartial(myInducingHom); }
    };

  private:
    class IdealImpl: public IdealBase
    {
    public:
      friend class GeneralQuotientRingImpl;
      IdealImpl(const QuotientRing& R, const vector<RingElem>& gens);
      ~IdealImpl();
      virtual IdealImpl* myClone() const;
//???    virtual void swap(ideal& other);

      virtual const QuotientRing& myRing() const; // covariant return type!
      virtual bool IamZero() const;
      virtual bool IamOne() const;
      virtual void mySetMaximalFlag(bool b) const;    // virtual + default impl
      virtual void mySetPrimeFlag(bool b) const;      // virtual + default impl

      virtual void myReduceMod(RingElemRawPtr rawx) const; // r elem of R, where I is ideal in R
      virtual bool IhaveElem(RingElemConstRawPtr rawx) const;
      virtual void myAdd(const ideal&);
      virtual void myMul(const ideal&);
      virtual void myIntersect(const ideal&);
      virtual void myColon(const ideal&);
      virtual bool myDivMod(RingElemRawPtr rawlhs, RingElemConstRawPtr rawnum, RingElemConstRawPtr rawden) const; // return true iff quot exists & is unique; if true then lhs = num/den modulo I

      virtual const std::vector<RingElem>& myGens() const; // gens as specified by user
      virtual const std::vector<RingElem>& myTidyGens() const; // tidier set of gens
//??? use default??      virtual void myOutputSelf(OpenMathOutput&) const;
    protected:
      virtual void myMaximalTest() const;
      virtual void myPrimeTest() const;

    private:
      QuotientRing myR;
      vector<RingElem> myGensValue;
      vector<RingElem> myTidyGensValue;
      ideal myPreimageIdeal;
      void SetTidyGensFromPreimageIdeal(); // implementation detail
      static const IdealImpl* ourGetPtr(const ideal& I);
    };
  };


  QuotientRingBase::QuotientRingBase(const ring& R, const ideal& I):
      myBaseRingValue(R),
      myDefiningIdeal(I),
      myReprRing(IsQuotientRing(R)?ReprRing(R):R),
      myReducingIdeal(IsQuotientRing(R)?(I+ReducingIdeal(R)):I)
  {
    CoCoA_ASSERT(RingOf(I) == R);
  }


  QuotientRing::QuotientRing(const QuotientRingBase* RingPtr):
      ring(RingPtr)
  {}


  QuotientRing::QuotientRing(const ring& R):
      ring(QuotientRingPtr(R, "QuotientRing ctor"))
  {}


  QuotientingHomImpl::QuotientingHomImpl(const QuotientRing& RmodI):
      RingHomBase(BaseRing(RmodI), RmodI)
  {}


  void QuotientingHomImpl::myApply(RingElemRawPtr rawimage, RingElemConstRawPtr rawarg) const
  {
    QuotientRingPtr(myCodomain)->myReduction(rawimage, rawarg);
  }


  void QuotientingHomImpl::myOutputSelfDetails(std::ostream& out) const
  {
    out << " canonical quotient";
  }



  GeneralQuotientRingImpl::GeneralQuotientRingImpl(const ring& R, const ideal& I):
      QuotientRingBase(R, I)
  {
    myRefCountInc();  // this is needed for exception cleanliness, in case one of the lines below throws
    // NOTE: currently this impl does not make its own versions of 0 and 1:
    //       it simply aliases those of myReprRing, and relabels them as belonging to this ring.
    myZeroPtr.reset(new RingElemAlias(ring(this), raw(zero(myReprRing))));
    myOnePtr.reset(new RingElemAlias(ring(this), raw(one(myReprRing))));
    myRefCountZero();
  }


  void GeneralQuotientRingImpl::myCharacteristic(BigInt& p) const
  {
    // Handle the simple case of Z/I
    if (IsZZ(myReprRing))
    {
      if (IsZero(myReducingIdeal)) p = 0;
      else IsInteger(p, TidyGens(myReducingIdeal)[0]);
      return;
    }
    if (IsOne(myReducingIdeal)) { p = 1; return; } // trivial ring has char=1
    if (IsPolyRing(myBaseRingValue) && IsField(CoeffRing(myBaseRingValue)))
    {
      myBaseRingValue->myCharacteristic(p);
      return;
    }
    CoCoA_ERROR(ERR::NYI, "GeneralQuotientRingImpl::characteristic");
    myBaseRingValue->myCharacteristic(p); // just to keep compiler quiet
  }


  long GeneralQuotientRingImpl::myLogCardinality() const
  {
    if (!IamFiniteField()) return 0;
    if (IsZZ(myReprRing)) return 1;
    if (IsPolyRing(myBaseRingValue))
      return CoeffRing(myBaseRingValue)->myLogCardinality() * len(QuotientBasis(myReducingIdeal)); // SLUG SLUG SLUG!!!
    CoCoA_ERROR(ERR::SERIOUS, "GeneralQuotientRingImpl::myLogCardinality");
    return -1;
  }


  bool GeneralQuotientRingImpl::IamFiniteField() const
  {
    if (IsZZ(myBaseRingValue)) //return IsProbPrime(myCharacteristic());
      /* STOPGAP IMPL*/ { BigInt P; myCharacteristic(P); return IsProbPrime(P); }

    // Currently we recognise only k[x]/I where k is finite field
    if (!IsPolyRing(myReprRing)) return false;
    if (!IsField(CoeffRing(myReprRing)))
      CoCoA_ERROR(ERR::NYI, "GeneralQuotientRingImpl::IamFiniteField");
    return IsFiniteField(CoeffRing(myReprRing)) &&
           IsZeroDim(myReducingIdeal) &&
           IsMaximal(myReducingIdeal);
  }


  RingElem GeneralQuotientRingImpl::myCanonicalRepr(ConstRawPtr rawx) const
  {
    RingElem ans(myReprRing);
    myReprRing->myAssign(raw(ans), rawx);
    return ans;
  }


  void GeneralQuotientRingImpl::myReduction(RawPtr rawimage, ConstRawPtr rawarg) const
  {
    myReprRing->myAssign(rawimage, rawarg);
    myReducingIdeal->myReduceMod(rawimage);
  }


  RingElemRawPtr GeneralQuotientRingImpl::myNew() const
  {
    return myReprRing->myNew();
  }


  RingElemRawPtr GeneralQuotientRingImpl::myNew(const MachineInt& n) const
  {
    RingElemRawPtr rawans = myReprRing->myNew(n);
    myReduction(rawans);
    return rawans;
  }


  RingElemRawPtr GeneralQuotientRingImpl::myNew(const BigInt& N) const
  {
    RingElemRawPtr rawans = myReprRing->myNew(N);
    myReduction(rawans);
    return rawans;
  }


  RingElemRawPtr GeneralQuotientRingImpl::myNew(ConstRawPtr rawCopyMe) const
  {
    return myReprRing->myNew(rawCopyMe);
  }


  void GeneralQuotientRingImpl::myDelete(RawPtr rawx) const
  {
    myReprRing->myDelete(rawx);
  }


  void GeneralQuotientRingImpl::mySwap(RawPtr rawx, RawPtr rawy) const
  {
    myReprRing->mySwap(rawx, rawy);
  }


  void GeneralQuotientRingImpl::myAssign(RawPtr rawlhs, ConstRawPtr rawx) const
  {
    myReprRing->myAssign(rawlhs, rawx);
  }


  void GeneralQuotientRingImpl::myAssign(RawPtr rawlhs, const MachineInt& n) const
  {
    myReprRing->myAssign(rawlhs, n);
    myReduction(rawlhs);
  }


  void GeneralQuotientRingImpl::myAssign(RawPtr rawlhs, const BigInt& N) const
  {
    myReprRing->myAssign(rawlhs, N);
    myReduction(rawlhs);
  }


  void GeneralQuotientRingImpl::myAssignZero(RawPtr rawlhs) const
  {
    myReprRing->myAssignZero(rawlhs); // assume zero is always reduced
  }


  void GeneralQuotientRingImpl::myRecvTwinFloat(RawPtr rawlhs, ConstRawPtr rawx) const
  {
    CoCoA_ASSERT(!IsExact(myReprRing));
    myReprRing->myRecvTwinFloat(rawlhs, rawx);
  }


  void GeneralQuotientRingImpl::myNegate(RawPtr rawlhs, ConstRawPtr rawx) const
  {
    myReprRing->myNegate(rawlhs, rawx);
    myReduction(rawlhs);
  }


  void GeneralQuotientRingImpl::myAdd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    myReprRing->myAdd(rawlhs, rawx, rawy);
    myReduction(rawlhs);
  }


  void GeneralQuotientRingImpl::mySub(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    myReprRing->mySub(rawlhs, rawx, rawy);
    myReduction(rawlhs);
  }


  void GeneralQuotientRingImpl::myMul(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    myReprRing->myMul(rawlhs, rawx, rawy);
    myReduction(rawlhs);
  }


  void GeneralQuotientRingImpl::myDiv(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    RingElem X = myCanonicalRepr(rawx);
    RingElem Y = myCanonicalRepr(rawy);
    RingElem ans(myReprRing);
    if (!myReducingIdeal->myDivMod(raw(ans), raw(X), raw(Y)))
      CoCoA_ERROR(ERR::BadQuot, "GQR::myDiv");
    myReduction(rawlhs, raw(ans));
  }


  bool GeneralQuotientRingImpl::myIsDivisible(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    if (myIsOne(rawy)) { myAssign(rawlhs, rawx); return true; }
    const RingElem X = myCanonicalRepr(rawx);
    const RingElem Y = myCanonicalRepr(rawy);
    RingElem ans(myReprRing);
    if (!myReducingIdeal->myDivMod(raw(ans), raw(X), raw(Y))) return false;
    myReduction(rawlhs, raw(ans));
    return true;
  }


  // We could simply use the default definition of this function.
  bool GeneralQuotientRingImpl::myIsInvertible(ConstRawPtr rawx) const
  {
    ring R(myReprRing);
    RingElem xcopy(R);
    R->myAssign(raw(xcopy), rawx);
    return IsOne(myReducingIdeal + ideal(xcopy));
  }


  //  This seems hard to define in full generality; e.g. there are some GCD domains which are not euclidean.
  void GeneralQuotientRingImpl::myGcd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
 {
   if (!IamTrueGCDDomain())
     CoCoA_ERROR(ERR::NotTrueGCDDomain, "gcd(lhs, x, y) in a quotient ring");
   if (IamField()) { myGcdInField(rawlhs, rawx, rawy); return; }

   //??? BUG INCOMPLETE BUG ???
   CoCoA_ERROR(ERR::NYI, "gcd in generic quotient ring which is not a field");
 }


  void GeneralQuotientRingImpl::myPowerSmallExp(RawPtr rawlhs, ConstRawPtr rawx, long n) const  // assumes n > 1
  {
    // Assert that we have a genuinely non-trivial case.
    CoCoA_ASSERT(n > 1);
    CoCoA_ASSERT(!myIsZero(rawx) && !myIsOne(rawx) && !myIsMinusOne(rawx));
    myBinaryPower(rawlhs, rawx, n); //??? BUG/SLUG this is not always the best choice ???
  }


  void GeneralQuotientRingImpl::myOutput(std::ostream& out, ConstRawPtr rawx) const
  {
    // THIS PRINT FORMAT IS UGLY!  BUG???
    RingElemAlias x(myReprRing, rawx); // alias for value to be printed
    if (!IsZZ(myReprRing))
    {
      out << "(" << x << ")";
      return;
    }
    if (DefaultResidueSetting()==GlobalSettings::NonNegResidues)
    {
      out << x;
      return;
    }
    // Special case of Z/(nZ) and using symmetric residues...
    BigInt n; myCharacteristic(n);
    if (x > n/2)
      out << x-n;
    else
      out << x;
  }


  bool GeneralQuotientRingImpl::myIsPrintAtom(ConstRawPtr rawx) const
  {
    if (myIsPrintedWithMinus(rawx))  return false;
    return true;  // either positive int or between parentheses
  }


  bool GeneralQuotientRingImpl::myIsPrintedWithMinus(ConstRawPtr rawx) const
  {
    if (IsZZ(myReprRing) && DefaultResidueSetting()==GlobalSettings::SymmResidues )
    {
      RingElemAlias x(myReprRing, rawx); // alias for value to be printed
      BigInt n; myCharacteristic(n);
      return x > n/2;
    }
    return false;    
  }
  

  void GeneralQuotientRingImpl::myOutputSelf(std::ostream& out) const
  {
    //    out << "QuotientRing(" << myBaseRingValue << ", " << myDefiningIdeal << ")";
    out << "RingWithID(" << myID 
        << ",\"RingWithID(" << ID(myBaseRingValue) << ")/"<<myDefiningIdeal<<"\")";
  }


  void GeneralQuotientRingImpl::myOutputSelfLong(std::ostream& out) const
  {
    out << "RingWithID(" << myID 
        << ",\"RingWithID(" << ID(myBaseRingValue) << ")/"<<myDefiningIdeal<<"\")";
    out <<"\n  with BaseRing  ";
    myBaseRingValue->myOutputSelfLong(out);
    out << std::endl;
  }


  void GeneralQuotientRingImpl::myOutputSelf(OpenMathOutput& OMOut) const
  {
    OMOut->mySendApplyStart();
    OMOut << OpenMathSymbol("cocoa", "QuotientRing");
    OMOut << myBaseRingValue;     // superfluous since the ideal already contains the ring info...
    OMOut << myDefiningIdeal;
    OMOut->mySendApplyEnd();
  }


  void GeneralQuotientRingImpl::myOutput(OpenMathOutput& OMOut, ConstRawPtr rawx) const
  {
    myReprRing->myOutput(OMOut, rawx);
  }


  bool GeneralQuotientRingImpl::myIsZero(ConstRawPtr rawx) const
  {
    return myReprRing->myIsZero(rawx);
  }


  bool GeneralQuotientRingImpl::myIsOne(ConstRawPtr rawx) const
  {
    return myReprRing->myIsOne(rawx) || IsOne(myReducingIdeal);
  }


  bool GeneralQuotientRingImpl::myIsMinusOne(ConstRawPtr rawx) const
  {
    if (myReprRing->myIsMinusOne(rawx) || IsOne(myReducingIdeal)) return true;
    RingElem tmp(myReprRing);
    myReprRing->myAdd(raw(tmp), raw(myOne()), rawx);
    return IsElem(tmp, myReducingIdeal);
  }


  bool GeneralQuotientRingImpl::myIsInteger(BigInt& N, ConstRawPtr rawx) const
  {
    return myReprRing->myIsInteger(N, rawx);
  }


  bool GeneralQuotientRingImpl::myIsRational(BigRat& Q, ConstRawPtr rawx) const
  {
    return myReprRing->myIsRational(Q, rawx);
  }


  bool GeneralQuotientRingImpl::myIsDouble(double& d, ConstRawPtr rawx) const
  {
    return myReprRing->myIsDouble(d, rawx);
  }


  bool GeneralQuotientRingImpl::myIsEqual(ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    return myReprRing->myIsEqual(rawx, rawy);
  }


  ideal GeneralQuotientRingImpl::myIdealCtor(const std::vector<RingElem>& gens) const
  {
    return ideal(new IdealImpl(QuotientRing(this), gens));
  }


  void GeneralQuotientRingImpl::myReduction(RawPtr rawx) const
  {
    myReducingIdeal->myReduceMod(rawx);
  }


  //---------------------------------------------------------------------------
  // Functions to do with homomorphisms


  RingHom GeneralQuotientRingImpl::myInducedHomCtor(const RingHom& InducingHom) const
  {
    CoCoA_ASSERT(IsInKer(myDefiningIdeal, InducingHom));
    return RingHom(new InducedHomImpl(QuotientRing(this), InducingHom));
  }


  RingHom GeneralQuotientRingImpl::myCompose(const RingHom& phi, const RingHom& theta) const
  {
    return myInducedHomCtor(phi(theta(myQuotientingHomCtor())));
  }


  bool GeneralQuotientRingImpl::myImageLiesInSubfield(const RingHom& phi) const
  {
    CoCoA_ASSERT(codomain(phi) == ring(this));
    return IamField() || myReprRing->myImageLiesInSubfield(phi); // ***BUG BUG BUG*** wrong codomain!!!!
  }


  void GeneralQuotientRingImpl::myPowerBigExp(RawPtr rawlhs, ConstRawPtr rawx, const BigInt& N) const
  {
    myBinaryPower(rawlhs, rawx, N);
  }


  GeneralQuotientRingImpl::InducedHomImpl::InducedHomImpl(const QuotientRing& RmodI, const RingHom& InducingHom):
      RingHomInducedBase(RmodI, InducingHom)
  { /* Compatibility already checked in InducedHom in QuotientRing.C */ }


  void GeneralQuotientRingImpl::InducedHomImpl::myApply(RawPtr rawimage, ConstRawPtr rawarg) const
  {
    myInducingHom->myApply(rawimage, rawarg);
  }


  //---------------------------------------------------------------------------
  // Functions to do with GQRIdeals

  inline const GeneralQuotientRingImpl::IdealImpl* GeneralQuotientRingImpl::IdealImpl::ourGetPtr(const ideal& I)
  {
    return dynamic_cast<const IdealImpl*>(I.myIdealPtr());
  }


  GeneralQuotientRingImpl::IdealImpl::IdealImpl(const QuotientRing& RmodI, const vector<RingElem>& gens):
      myR(RmodI),
      myGensValue(gens),
      myPreimageIdeal(ReducingIdeal(RmodI))
  {
    vector<RingElem> OverGens;
    const long ngens = len(gens);
    for (long i=0; i < ngens; ++i)
      OverGens.push_back(RmodI->myCanonicalRepr(raw(gens[i])));
    myPreimageIdeal += ideal(ReprRing(RmodI), OverGens);
    SetTidyGensFromPreimageIdeal();
    //    IamPrime3Flag = uncertain3;   // default value
    //    IamMaximal3Flag = uncertain3; // default value
  }


  GeneralQuotientRingImpl::IdealImpl::~IdealImpl()
  {}


  GeneralQuotientRingImpl::IdealImpl* GeneralQuotientRingImpl::IdealImpl::myClone() const
  {
    return new IdealImpl(*this);
  }


  // NB the fowarded call to myPreimageIdeal does all the consistency checking
  void GeneralQuotientRingImpl::IdealImpl::mySetPrimeFlag(bool b) const
  {
    myPreimageIdeal->mySetPrimeFlag(b);
    IamPrime3Flag = b;
  }


  // NB the fowarded call to myPreimageIdeal does all the consistency checking
  void GeneralQuotientRingImpl::IdealImpl::mySetMaximalFlag(bool b) const
  {
    myPreimageIdeal->mySetMaximalFlag(b);
    IamMaximal3Flag = b;
  }


  const QuotientRing& GeneralQuotientRingImpl::IdealImpl::myRing() const
  {
    return myR;
  }


  bool GeneralQuotientRingImpl::IdealImpl::IamZero() const
  {
    return myPreimageIdeal == ReducingIdeal(myR);
  }


  bool GeneralQuotientRingImpl::IdealImpl::IamOne() const
  {
    return IsOne(myPreimageIdeal);
  }


  void GeneralQuotientRingImpl::IdealImpl::myMaximalTest() const
  {
    CoCoA_ERROR(ERR::NYI, "GeneralQuotientRingImpl::IdealImpl::myMaximalTest()");
    //???    mySetMaximalFlag(IsMaximal(myPreimageIdeal));
  }


  void GeneralQuotientRingImpl::IdealImpl::myPrimeTest() const
  {
    CoCoA_ERROR(ERR::NYI, "GeneralQuotientRingImpl::IdealImpl::myPrimeTest()");
    //???    mySetPrimeFlag(IsPrime(myPreimageIdeal));
  }


  void GeneralQuotientRingImpl::IdealImpl::myReduceMod(RawPtr rawx) const
  {
    if (IamZero()) return;
    myPreimageIdeal->myReduceMod(rawx);
  }


  bool GeneralQuotientRingImpl::IdealImpl::IhaveElem(RingElemConstRawPtr rawx) const
  {
    return myPreimageIdeal->IhaveElem(rawx);
  }


  void GeneralQuotientRingImpl::IdealImpl::myAdd(const ideal& J)
  {
    myGensValue.insert(myGensValue.end(), gens(J).begin(), gens(J).end());
    myPreimageIdeal += ourGetPtr(J)->myPreimageIdeal;
    SetTidyGensFromPreimageIdeal();
    IamPrime3Flag = ourGetPtr(J)->IamPrime3Flag;
    IamMaximal3Flag = ourGetPtr(J)->IamMaximal3Flag;
  }


  void GeneralQuotientRingImpl::IdealImpl::myMul(const ideal& J)
  {
    if (IamZero()) return;
    //    CoCoA_ERROR(ERR::NYI, "GeneralQuotientRingImpl::IdealImpl::mul");
    vector<RingElem> tmpV;
    for (vector<RingElem>::const_iterator itI=myGens().begin(); itI!=myGens().end(); ++itI)
      for (vector<RingElem>::const_iterator itJ=gens(J).begin(); itJ!=gens(J).end(); ++itJ)
        tmpV.push_back((*itI)*(*itJ));
    swap(tmpV, myGensValue);
  }


  void GeneralQuotientRingImpl::IdealImpl::myIntersect(const ideal& J)
  {
    MakeUnique(myPreimageIdeal)->myIntersect(ourGetPtr(J)->myPreimageIdeal);
    SetTidyGensFromPreimageIdeal();
    myGensValue = myTidyGensValue;
    IamPrime3Flag = ourGetPtr(J)->IamPrime3Flag;
    IamMaximal3Flag = ourGetPtr(J)->IamMaximal3Flag;
  }


  void GeneralQuotientRingImpl::IdealImpl::myColon(const ideal& J)
  {
    MakeUnique(myPreimageIdeal)->myColon(ourGetPtr(J)->myPreimageIdeal);
    SetTidyGensFromPreimageIdeal();
    myGensValue = myTidyGensValue;
    IamPrime3Flag = ourGetPtr(J)->IamPrime3Flag;
    IamMaximal3Flag = ourGetPtr(J)->IamMaximal3Flag;
    //??? CANNOT WORK OUT HOW TO USE THE ALGORITHM transform OR for_each :-(
    //transform(myPreimageIdeal->gens().begin(), myPreimageIdeal->gens().end(), back_inserter(myGens), myR->myReduction);
  }


  bool GeneralQuotientRingImpl::IdealImpl::myDivMod(RingElemRawPtr rawlhs, RingElemConstRawPtr rawnum, RingElemConstRawPtr rawden) const
  {
    const QuotientRing& RmodI = myR;
    const ring& R = ReprRing(RmodI);
    const RingElem N = RmodI->myCanonicalRepr(rawnum);
    const RingElem D = RmodI->myCanonicalRepr(rawden);
    RingElem ans(R);
    if (!myPreimageIdeal->myDivMod(raw(ans), raw(N), raw(D))) return false;
    RmodI->myReduction(rawlhs, raw(ans));
    return true;
  }


  const std::vector<RingElem>& GeneralQuotientRingImpl::IdealImpl::myGens() const
  {
    return myGensValue;
  }


  const std::vector<RingElem>& GeneralQuotientRingImpl::IdealImpl::myTidyGens() const
  {
    return myTidyGensValue;
  }


  void GeneralQuotientRingImpl::IdealImpl::SetTidyGensFromPreimageIdeal()
  {
    vector<RingElem> NewTidyGens;
    const vector<RingElem>& OverGens = TidyGens(myPreimageIdeal);
    const long ngens = len(OverGens);
    for (long i=0; i < ngens; ++i)
    {
      RingElem g(myR);
      myR->myReduction(raw(g), raw(OverGens[i]));
      if (!IsZero(g))
        NewTidyGens.push_back(g);
    }
    myTidyGensValue = NewTidyGens;
  }


  //---------------------------------------------------------------------------

  QuotientRing NewQuotientRing(const ring& R, const ideal& I)
  {
    if (R != RingOf(I)) CoCoA_ERROR(ERR::IdealNotInRing, "NewQuotientRing");
    if (IsOne(I)) CoCoA_ERROR(ERR::BadQuotRing, "NewQuotientRing");
//???    if (IsZero(I)) CoCoA_ERROR(ERR::BadQuotRing, "NewQuotientRing"); ???
    if (!IsZero(I) && IsZZ(R))
    {
      // Try to make a RingFp first as it is fastest.
      if (IsGoodForRingFp(I)) return NewRingFp(I);
      // Couldn't make a RingFp, so try to make a RingFpDouble.
      if (IsGoodForRingFpDouble(I)) return NewRingFpDouble(I);
      // Couldn't make RingFp or RingFpDouble (char not prime or too big),
      // so fall through to (relatively slow) generic implementation.
    }
    return QuotientRing(new GeneralQuotientRingImpl(R, I));
  }


  QuotientRing NewZZmod(const MachineInt& n)
  {
    return NewZZmod(BigInt(n));
  }


  QuotientRing NewZZmod(const BigInt& N)
  {
    if (N == 1 || N == -1) CoCoA_ERROR(ERR::BadQuotRing, "NewZZmod");
    if (IsZero(N)) CoCoA_ERROR(ERR::BadQuotRing, "NewZZmod"); //??? disallow N==0???
    return NewQuotientRing(RingZZ(),
                           ideal(RingElem(RingZZ(), N)));
  }


  RingHom QuotientingHom(const QuotientRing& RmodI)
  {
    return RmodI->myQuotientingHomCtor();
  }


  RingHom InducedHom(const QuotientRing& RmodI, const RingHom& phi)
  {
    if (domain(phi) != BaseRing(RmodI))
      CoCoA_ERROR(ERR::BadInducingHom, "InducedHom(QuotientRing,hom)");

    if (!IsInKer(DefiningIdeal(RmodI), phi))
      CoCoA_ERROR(ERR::BadInducingHomKer, "InducedHom(QuotientRing,hom)");

    return RmodI->myInducedHomCtor(phi);
  }


  RingElem CanonicalRepr(ConstRefRingElem r)
  {
    if (!IsQuotientRing(owner(r)))
      CoCoA_ERROR(ERR::NotElemQuotientRing, "CanonicalRepr(r)");
    return QuotientRingPtr(owner(r))->myCanonicalRepr(raw(r));
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/QuotientRing.C,v 1.54 2014/07/30 14:08:37 abbott Exp $
// $Log: QuotientRing.C,v $
// Revision 1.54  2014/07/30 14:08:37  abbott
// Summary: Changed name AmbientRing --> RingOf; myAmbientRingValue --> myR
// Author: JAA
//
// Revision 1.53  2014/07/28 16:04:28  abbott
// Summary: Renamed myQuotientingHom to myQuotientingHomCtor
// Author: JAA
//
// Revision 1.52  2014/07/28 15:46:25  abbott
// Summary: Redesign: ringhoms no longer cached in rings (caused ref count trouble); myQuotientingHomCtor
// Author: JAA
//
// Revision 1.51  2014/07/11 15:44:13  bigatti
// -- added  myOutputSelfLong
//
// Revision 1.50  2014/07/09 13:00:19  abbott
// Summary: Removed some cruft
// Author: JAA
//
// Revision 1.49  2014/07/08 15:22:35  abbott
// Summary: Added new ctor (from general ring)
// Author: JAA
//
// Revision 1.48  2014/07/08 13:14:40  abbott
// Summary: Removed AsQuotientRing; added new defn of BaseRing
// Author: JAA
//
// Revision 1.47  2014/07/07 12:43:25  abbott
// Summary: Removed AsPolyRing
// Author: JAA
//
// Revision 1.46  2014/07/04 13:08:08  bigatti
// -- RingID into RingWithID
//
// Revision 1.45  2014/07/02 16:53:16  bigatti
// -- new way of printing ring with ID
//
// Revision 1.44  2014/05/14 15:57:15  bigatti
// -- added "using" for clang with superpedantic flag
//
// Revision 1.43  2014/04/15 13:29:13  abbott
// Summary: Corrected grammar in a comment
// Author: JAA
//
// Revision 1.42  2014/04/02 10:57:46  abbott
// Summary: Revised design of IamIntegralDomain3
// Author: JAA
//
// Revision 1.41  2014/03/28 18:23:02  bigatti
// -- added myMul for ideals
//
// Revision 1.40  2014/03/27 17:17:31  abbott
// Summary: Added new fn IsIntegralDomain3 (and mem fn IamIntegralDomain3)
// Author: JAA
//
// Revision 1.39  2013/06/28 17:03:51  abbott
// Modified semantics of IdealBase::myDivMod;
// it now returns a boolean.
// Several consequential changes.
//
// Revision 1.38  2013/05/29 17:11:40  bigatti
// -- fixed IdealImpl::myMaximalTest IdealImpl::myPrimeTest (NYI)
// -- inlined some one-liners
//
// Revision 1.37  2013/02/21 14:14:42  abbott
// First attempt at implementing PartialRingHom -- some problems remain!!
//
// Revision 1.36  2012/10/24 13:36:36  abbott
// IMPORTANT CHANGE:
//  * zero/one are now held inside std::auto_ptr objects;
//  * consequent changes to ring ctor & to myZero/myOne fns.
// Replaced ConstRefRingElem by RingElemAlias in 2 ctor calls.
//
// Revision 1.35  2012/05/30 16:04:55  bigatti
// -- applied "3" convention on bool3 functions and member fields
//
// Revision 1.34  2012/05/28 10:35:54  abbott
// Corrected defn of IsTrueGCDDomain.
//
// Revision 1.33  2012/05/28 09:18:21  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.32  2012/05/22 10:02:37  abbott
// Removed IsGCDDomain; substituted by IsTrueGCDDomain.
// Added IsFractionFieldOfGCDDomain.
//
// Revision 1.31  2012/05/07 14:32:45  abbott
// Improved defn of  GeneralQuotientRingImpl::IamGCDDomain.
//
// Revision 1.30  2012/04/27 15:05:24  abbott
// Added mem fns IamFiniteField &  myLogCardinality
//
// Revision 1.29  2012/02/10 10:28:08  bigatti
// -- changed RingZ.H, RingQ.H --> RingZZ.H, RingQQ.H
//
// Revision 1.28  2012/02/08 16:14:07  bigatti
// -- changed Z,Q --> ZZ,QQ
//
// Revision 1.27  2012/01/26 16:51:13  bigatti
// -- changed back_inserter into insert
//
// Revision 1.26  2011/11/09 14:09:53  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.25  2011/08/24 10:28:49  bigatti
// -- renamed QQ --> BigRat
// -- sorted #include
//
// Revision 1.24  2011/08/14 15:52:17  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.23  2011/08/12 15:22:15  abbott
// Added two missing const keywords.
//
// Revision 1.22  2011/06/23 16:04:47  abbott
// Added IamExact mem fn for rings.
// Added myRecvTwinFloat mem fn for rings.
// Added first imple of RingHom from RingTwinFloat to other rings.
//
// Revision 1.21  2011/05/24 14:54:29  abbott
// Consquential changes from removing several ctors for principal ideals.
//
// Revision 1.20  2011/05/19 14:42:09  abbott
// Replaced calls to DefaultResiduesAreSymm by calls to DefaultResidueSetting.
// Added myIsDouble.
//
// Revision 1.19  2011/03/22 20:20:19  abbott
// Now uses the new fns for testing whether an args is suitable for
// pseudo-ctor for RingFp, RingFpLog, RingFpDouble.
//
// Revision 1.18  2011/03/11 21:52:00  abbott
// Added defn of myPowerBigExp (o/w test-RingFp1 fails).
//
// Revision 1.17  2011/03/10 16:39:34  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.16  2010/10/06 14:10:24  abbott
// Added increments to the ref count in ring and PPMonoid ctors to make
// them exception safe.
//
// Revision 1.15  2010/10/01 15:45:17  bigatti
// -- added mySymbolValue
//
// Revision 1.14  2010/06/10 08:00:02  bigatti
// -- fixed naming conventions
//
// Revision 1.13  2010/02/04 09:57:11  bigatti
// -- added "mul" for ideals.  Implemented only for SparsePolyRing
//
// Revision 1.12  2009/12/03 17:40:36  abbott
// Added include directives for ZZ.H (as a consequence of removing
// the directive from ring.H).
//
// Revision 1.11  2009/09/24 15:36:00  abbott
// Removed some unnecessary "std::" prefixes.
// Removed some pointless include directives.
//
// Revision 1.10  2009/07/24 14:19:23  abbott
// Cleaned up include directives, and added fwd decls to the header file.
//
// Revision 1.9  2009/07/02 16:32:11  abbott
// Consequential changes stemming from new class BigRat, and modified interface to the member
// function RingBase::myIsRational.  Also some new conversion functions.
//
// Revision 1.8  2009/05/22 10:31:41  bigatti
// -- added myIsPrintAtom, myIsPrintedWithMinus
// -- improved myOutput
//
// Revision 1.7  2009/05/14 09:39:29  abbott
// Added possibility to specify "symmetric" or "non-negative" residues
// in quotients of ZZ.  Affects printing of elements in quotients of ZZ
// (also changed printing of elements in general quotient rings).
// Consequent changes in several tests.
//
// Revision 1.6  2008/12/17 12:11:52  abbott
// Changed type from long to MachineInt in operations which use a machine integer
// in place of a RingElem.  The change is "superficial" but affects many files.
//
// Revision 1.5  2008/10/07 15:45:22  abbott
// Changed ErrorInfo objects so they include the name of their own error ID.
// Changed catch statements to catch const objects.
// Removed calls to the member fn which accessed the error ID member of an
// ErrorInfo; now you simply compare directly with the error ID (makes the
// code easier to read).
//
// Revision 1.4  2007/10/30 17:14:07  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.3  2007/09/25 16:32:30  abbott
// Several minor changes to silence gcc-4.3:
//    more #includes,
//    and fixed a template problemm in RegisterServerOps.C
//
// Revision 1.2  2007/05/22 22:45:14  abbott
// Changed fn name IsUnit to IsInvertible.
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
// Revision 1.11  2007/03/02 10:47:53  cocoa
// First stage of RingZ modifications -- tests do not compile currently, Anna will fix this.
//
// Revision 1.10  2007/02/26 15:00:54  bigatti
// -- added placeholders for new syntax based on unique Z implementation
//
// Revision 1.9  2007/01/15 15:47:57  cocoa
// -- added prefix "raw" to RawPtr arguments names
// -- changed rhs into rawx, n, or N
//
// Revision 1.8  2007/01/13 14:14:34  cocoa
// Overhaul of RingHom code: it nows uses SmartPtrIRC, and printing is more logical.
// Have not yet updated the documentation.
//
// Revision 1.7  2006/12/21 13:48:32  cocoa
// Made all increment/decrement calls prefix (except where the must be postfix).
//
// Revision 1.6  2006/12/18 15:16:05  cocoa
// Added definition of CanonRepr.
//
// Revision 1.5  2006/12/07 11:57:53  cocoa
// -- style: RawPtr args are now called "raw.."
//
// Revision 1.4  2006/11/08 16:21:59  cocoa
// Structural cleaning of RingHom; many consequential changes.
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
// Revision 1.9  2006/05/29 16:22:37  cocoa
// Third time lucky???
// Added myIsInteger member function to all rings (NYI for RingFloat).
//
// Revision 1.8  2006/05/12 16:10:58  cocoa
// Added OpenMathFwd.H, and tidied OpenMath.H.
// Many consequential but trivial changes.
//
// Revision 1.7  2006/04/21 15:01:36  cocoa
// Changed default implementation of RingBase::myGcd -- it now gives a SERIOUS
// error.  All fields must now handle a call to gcd explicitly: they can use
// the new myGcdInField function.  It's now cleaner than it was.
//
// Revision 1.6  2006/03/27 12:21:25  cocoa
// Minor silly changes to reduce number of complaints from some compiler or other.
//
// Revision 1.5  2006/03/21 09:43:14  cocoa
// Changed names of some member fns of ideals (dealing with setting and testing
// the flags for primeness and maximality).  Hope icc will complain less now.
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
// Revision 1.5  2005/10/14 15:25:07  cocoa
// Major tidying and cleaning to small prime finite fields.
// Several consequential changes.  Improved their documentation.
//
// Added Makefile and script to include/CoCoA/ directory to
// keep library.H up to date.
//
// Revision 1.4  2005/10/11 16:37:30  cocoa
// Added new small prime finite field class (see RingFpDouble).
//
// Cleaned makefiles and configuration script.
//
// Tidied PPMonoid code (to eliminate compiler warnings).
//
// Fixed bug in RingFloat::myIsInteger.
//
// Revision 1.3  2005/09/30 15:03:39  cocoa
// Minor cleaning and tidying.
// DistrMPolyInlPP: use of summands now rather cleaner.
//
// Revision 1.2  2005/07/08 15:09:28  cocoa
// Added new symbol class (to represent names of indets).
// Integrated the new class into concrete polynomial rings
// and PPMonoid -- many consequential changes.
// Change ctors for the "inline" sparse poly rings: they no
// longer expect a PPMonoid, but build their own instead
// (has to be a PPMonoidOv).
//
// Revision 1.1.1.1  2005/05/03 15:47:31  cocoa
// Imported files
//
// Revision 1.5  2005/04/20 15:40:48  cocoa
// Major change: modified the standard way errors are to be signalled
// (now via a macro which records filename and line number).  Updated
// documentation in error.txt accordingly.
//
// Improved the documentation in matrix.txt (still more work to be done).
//
// Revision 1.4  2005/04/19 14:06:04  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.3  2005/04/01 16:18:20  cocoa
// Friday check-in.  Fixed a bug in the ctor for GeneralQuotientRingImpl.
//
// Revision 1.2  2005/02/11 14:15:20  cocoa
// New style ring elements and references to ring elements;
// I hope I have finally got it right!
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.19  2004/11/19 15:44:27  cocoa
// Changed names of "casting" functions which convert a ring into
// one with a more special structure (e.g. FractionField).  These
// functions now have names starting with "As".  There were several
// consequential changes.
//
// Revision 1.18  2004/11/18 18:33:41  cocoa
// Now every ring know its own "one" element (as well as "zero").
// Several consequential changes.
//
// Revision 1.17  2004/11/12 15:49:29  cocoa
// Tidying prior to 0.90 release.
// (a) documentation improved (or marked as poor)
// (b) sundry minor improvements to the code
//
// Revision 1.16  2004/11/09 15:50:54  cocoa
// -- changed GeneralQuotientRing --> GeneralQuotientRingImpl
// -- added: NewZmod(Z, N)  with N in ZZ
//
// Revision 1.15  2004/11/04 18:47:43  cocoa
// (1) Ring member functions which previously expected mpz_t args
//     now expect ZZ args.  Numerous minor consequential changes.
// (2) Renamed function which gives access to the mpz_t value inside
//     a ZZ object: previously was raw(...), now is mpzref(...).
//     Plenty of calls had to be altered.
//
// Revision 1.14  2004/07/27 16:03:39  cocoa
// Added IsCommutative test and IamCommutative member function
// to all rings.  Tidied geobuckets a little.
//
// Revision 1.13  2004/07/16 15:45:12  cocoa
// First stage of new RingElem implementation completed.
//
// Revision 1.12  2004/07/13 16:32:27  cocoa
// First stage of major revamp of ring elements.
// Implementation of RingFp has been split into "ring interface"
// and "algorithms plus data structures".
//
// Revision 1.11  2004/06/30 16:46:06  cocoa
// End of day check-in: mostly tidying up, and correcting some
// overly lax access permissions.
//
// Revision 1.10  2004/06/29 17:10:22  cocoa
// Partially tidied use of "protected" and "private" in various
// base classes.  Checking in at the end of the day -- it works,
// and I wouldn't want it to be lost next time point's disk
// misbehaves.
//
// Revision 1.9  2004/06/29 15:34:29  cocoa
// Tidied the functions dealing with bool3 values.
//
// Revision 1.8  2004/05/24 15:52:14  cocoa
// Major update:
//   new error mechanism
//   many fixes
//   RingHoms almost work now
//   RingFloat much improved
//
// Revision 1.7  2004/04/08 15:33:34  cocoa
// Added function IsInteger, and the related RingBase::myIsInteger
// virtual function, plus all necessary implementations.
//
// Revision 1.6  2004/03/20 17:46:10  cocoa
// Check in prior to departure to RWCA
//
// Revision 1.5  2004/02/03 16:16:21  cocoa
// Removed pointless IamGCDDomain functions from several concrete rings.
// Added IamOrderedDomain functions where appropriate.
// Tidied ctors for the small finite fields.
//
// Revision 1.4  2004/01/28 15:32:48  cocoa
// Tidied use of C++ std library components, and reordered some #includes.
//
// Revision 1.3  2003/10/17 10:51:06  cocoa
// Major cleaning, and new naming convention.
//
// Revision 1.2  2003/10/09 12:16:38  cocoa
// New coding convention for rings.
//
// Revision 1.2  2003/06/23 16:48:54  abbott
// Minor cleaning prior to public release.
// Mostly consequential changes.
// Fixed NewQuotientRing to work properly in limit cases (e.g. R/ideal(0)).
//
// Revision 1.1  2003/05/14 17:12:01  abbott
// Initial revision
//
//
