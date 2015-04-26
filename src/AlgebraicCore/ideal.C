//   Copyright (c)  2005  John Abbott

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


// Source code for classes ideal and IdealBase

#include "CoCoA/ideal.H"

#include "CoCoA/OpenMath.H"
#include "CoCoA/convert.H"
#include "CoCoA/error.H"
#include "CoCoA/VectorOperations.H"  // for HasUniqueOwner
#include "CoCoA/utils.H"  // for len

#include <algorithm>
using std::copy;
#include <iostream>
using std::ostream;
//#include <vector>
using std::vector;


namespace CoCoA
{

  // C++ needs this function to be defined
  IdealBase::~IdealBase()
  {}


  //---------------------------------------------------------------------------

  ideal::ideal(IdealBase* IPtr):
      myPtr(IPtr)
  {
    IPtr->myRefCountInc();
  }


  ideal::ideal(ConstRefRingElem r)
  {
    vector<RingElem> gens;
    gens.push_back(RingElem(r));
    ideal tmp = owner(r)->myIdealCtor(gens);
    myPtr = tmp.myPtr;
    myPtr->myRefCountInc();
  }


  ideal::ideal(ConstRefRingElem r1, ConstRefRingElem r2)
  {
    vector<RingElem> gens;
    gens.push_back(RingElem(r1));
    gens.push_back(RingElem(r2));
    if (!HasUniqueOwner(gens))
      CoCoA_ERROR(ERR::MixedRings, "ideal(r1, r2)");
    ideal tmp = owner(r1)->myIdealCtor(gens);
    myPtr = tmp.myPtr;
    myPtr->myRefCountInc();
  }


  ideal::ideal(ConstRefRingElem r1, ConstRefRingElem r2, ConstRefRingElem r3)
  {
    vector<RingElem> gens;
    gens.push_back(RingElem(r1));
    gens.push_back(RingElem(r2));
    gens.push_back(RingElem(r3));
    if (!HasUniqueOwner(gens))
      CoCoA_ERROR(ERR::MixedRings, "ideal(r1, r2, r3)");
    ideal tmp = owner(r1)->myIdealCtor(gens);
    myPtr = tmp.myPtr;
    myPtr->myRefCountInc();
  }


  ideal::ideal(ConstRefRingElem r1, ConstRefRingElem r2, ConstRefRingElem r3, ConstRefRingElem r4)
  {
    vector<RingElem> gens;
    gens.push_back(RingElem(r1));
    gens.push_back(RingElem(r2));
    gens.push_back(RingElem(r3));
    gens.push_back(RingElem(r4));
    if (!HasUniqueOwner(gens))
      CoCoA_ERROR(ERR::MixedRings, "ideal(r1, r2, r3, r4)");
    ideal tmp = owner(r1)->myIdealCtor(gens);
    myPtr = tmp.myPtr;
    myPtr->myRefCountInc();
  }


  ideal::ideal(const std::vector<RingElem>& gens)
  {
    if (gens.empty()) CoCoA_ERROR(ERR::Empty, "ideal(gens)");
    if (!HasUniqueOwner(gens)) CoCoA_ERROR(ERR::MixedRings, "ideal(gens)");
    ideal tmp = owner(gens[0])->myIdealCtor(gens);
    myPtr = tmp.myPtr;
    myPtr->myRefCountInc();
  }


  ideal::ideal(const ring& R, const std::vector<RingElem>& gens)
  {
    if (!gens.empty())
      if (owner(gens[0]) != R || !HasUniqueOwner(gens))
        CoCoA_ERROR(ERR::MixedRings, "ideal(R, gens)");
    ideal tmp = R->myIdealCtor(gens);
    myPtr = tmp.myPtr;
    myPtr->myRefCountInc();
  }


  ideal::ideal(const ideal& copy)
  {
    copy->myRefCountInc();
    myPtr = copy.myPtr;
  }


  ideal& ideal::operator=(const ideal& rhs)
  {
    if (RingOf(*this) != RingOf(rhs))
      CoCoA_ERROR(ERR::MixedRings, "ideal=ideal (assignment)");
    rhs->myRefCountInc();
    myPtr->myRefCountDec();
    myPtr = rhs.myPtr;
    return *this;
  }


  ideal::~ideal()
  {
    myPtr->myRefCountDec();
  }


  IdealBase* MakeUnique(ideal& I)
  {
    if (I->myRefCountIsOne()) return I.myPtr;
    IdealBase* NewPtr(I->myClone());
    I->myRefCountDec();  // after myClone for exc. safety
    I.myPtr = NewPtr;
    I.myPtr->myRefCount = 1;
    return I.myPtr;
  }


  //---------------------------------------------------------------------------
  // Syntactic sugar functions

  RingElem operator%(const MachineInt& n, const ideal& I)
  {
    RingElem ans(RingOf(I), n);
    I->myReduceMod(raw(ans));
    return ans;
  }


  RingElem operator%(const BigInt& N, const ideal& I)
  {
    RingElem ans(RingOf(I), N);
    I->myReduceMod(raw(ans));
    return ans;
  }


  RingElem operator%(const BigRat& Q, const ideal& I)
  {
    RingElem ans(RingOf(I), Q);
    I->myReduceMod(raw(ans));
    return ans;
  }


  RingElem operator%(ConstRefRingElem r, const ideal& I)
  {
    if (owner(r) != RingOf(I))
      CoCoA_ERROR(ERR::MixedRings, "r%I  -- reduction of RingElem modulo an ideal");
    RingElem ans(r);
    I->myReduceMod(raw(ans));
    return ans;
  }


  // Two separate impls in case ring is not commutative
  const ideal operator*(ConstRefRingElem r, const ideal& I)
  {
    if (owner(r) != RingOf(I))
      CoCoA_ERROR(ERR::MixedRings, "r*I, product of RingElem and ideal");
    if (IsZero(r)) return ideal(r);
    if (IsInvertible(r)) return I;
    vector<RingElem> g = gens(I);
    const long n = len(g);
    for (long i=0; i < n; ++i)
      g[i] = r*g[i];
    return ideal(g);
  }

  const ideal operator*(const ideal& I, ConstRefRingElem r)
  {
    if (owner(r) != RingOf(I))
      CoCoA_ERROR(ERR::MixedRings, "I*r, product of ideal and RingElem");
    if (IsZero(r)) return ideal(r);
    if (IsInvertible(r)) return I;
    vector<RingElem> g = gens(I);
    const long n = len(g);
    for (long i=0; i < n; ++i)
      g[i] = g[i]*r;
    return ideal(g);
  }


  const ideal operator+(const ideal& I, const ideal& J)
  {
    if (RingOf(I) != RingOf(J))
      CoCoA_ERROR(ERR::MixedRings, "ideal+ideal");
    ideal ans(I);
    MakeUnique(ans)->myAdd(J);
    return ans;
  }


  ideal& operator+=(ideal& I, const ideal& J)
  {
    if (RingOf(I) != RingOf(J))
      CoCoA_ERROR(ERR::MixedRings, "ideal+=ideal");
    MakeUnique(I)->myAdd(J);
    return I;
  }


  const ideal operator*(const ideal& I, const ideal& J)
  {
    if (RingOf(I) != RingOf(J))
      CoCoA_ERROR(ERR::MixedRings, "ideal*ideal");
    ideal ans(I);
    MakeUnique(ans)->myMul(J);
    return ans;
  }


  ideal& operator*=(ideal& I, const ideal& J)
  {
    if (RingOf(I) != RingOf(J))
      CoCoA_ERROR(ERR::MixedRings, "ideal*=ideal");
    MakeUnique(I)->myMul(J);
    return I;
  }


  ideal power(const ideal& I, const MachineInt& n)
  {
    const ring R = RingOf(I);
    if (IsNegative(n)) CoCoA_ERROR(ERR::NotNonNegative, "power(I,n)");
    unsigned long N = AsUnsignedLong(n);
    if (N == 0) return ideal(one(R));

    // An iterative implementation of binary powering.
    unsigned long bit = 1; while (bit <= N/2) bit <<= 1;
    ideal ans = I;
    while (bit != 1)
    {
      ans *= ans;
      bit >>= 1;
      if (N&bit) ans *= I;
    }
    return ans;
  }


  ideal power(const ideal& I, const BigInt& N)
  {
    long n;
    if (!IsConvertible(n, N)) CoCoA_ERROR(ERR::ExpTooBig, "power(I,N)");
    return power(I, n);
  }


  const ideal intersect(const ideal& I, const ideal& J)
  {
    if (RingOf(I) != RingOf(J))
      CoCoA_ERROR(ERR::MixedRings, "intersect(ideal,ideal)");
    if (IsZero(I) || IsZero(J))
      return ideal(RingOf(I), vector<RingElem>(0));
    // case IsOne in general can be very expensive
    ideal ans(I);
    MakeUnique(ans)->myIntersect(J);
    return ans;
  }


  const ideal colon(const ideal& I, const ideal& J)
  {
    if (RingOf(I) != RingOf(J))
      CoCoA_ERROR(ERR::MixedRings, "colon(ideal,ideal)");
    ideal ans(I);
    MakeUnique(ans)->myColon(J);
    return ans;
  }


  bool IsContained(const ideal& I, const ideal& J)
  {
    const vector<RingElem>& gensI = gens(I);
    const long NumGensI = len(gensI);
    for (long i=0; i < NumGensI; ++i)
      if (!IsElem(gensI[i], J)) return false;
    return true;
  }


  bool operator==(const ideal& I, const ideal& J)
  {
    //??? check first whether the myIdealPtrs are equal???
    return IsContained(I, J) && IsContained(J, I);
  }


  //---------------------------------------------------------------------------
  // Functions to do with IdealBase
  // The functions for user assertions must check consistency of the assertions:


  bool IdealBase::IamOne() const
  {
    return IhaveElem(raw(one(myRing())));
  }

  //   maximality implies primality

  bool IdealBase::IamMaximal() const
  {
    if (IsUncertain3(IamMaximal3Flag))
      myMaximalTest();
    return IsTrue3(IamMaximal3Flag);
  }


  void IdealBase::mySetMaximalFlag(bool b) const
  {
    if (b)
    { // user asserts maximal
      if (IsFalse3(IamMaximal3Flag))
        CoCoA_ERROR("Contradictory user assertion", "mySetMaximalFlag");
      IamMaximal3Flag = true3;
      IamPrime3Flag = true3;
      return;
    }

    // Otherwise user asserts not maximal...
    if (IsTrue3(IamMaximal3Flag))
      CoCoA_ERROR("Contradictory user assertion", "mySetMaximalFlag");
    IamMaximal3Flag = false3;
  }


  bool IdealBase::IamPrime() const
  {
    if (IsUncertain3(IamPrime3Flag))
      myPrimeTest();
    return IsTrue3(IamPrime3Flag);
  }


  void IdealBase::mySetPrimeFlag(bool b) const
  {
    if (b)
    { // user asserts prime...
      if (IsFalse3(IamPrime3Flag))
        CoCoA_ERROR("Contradictory user assertion", "IamPrime");
      IamPrime3Flag = true3;
      return;
    }

    // Otherwise user asserts not prime...
    if (IsTrue3(IamPrime3Flag))
      CoCoA_ERROR("Contradictory user assertion", "IamPrime");
    IamPrime3Flag = false3;
    IamMaximal3Flag = false3;
  }


  // Simplistic default definition
  void IdealBase::myOutputSelf(std::ostream& out) const
  {
    out << "ideal(";
    const vector<RingElem>& g = myGens();
    const long NumGens = len(g);
    for (long i=0; i < NumGens; ++i)
    {
      out << g[i];
      if (i != NumGens-1) out << ", ";
    }
    out << ")";
  }


  std::ostream& operator<<(std::ostream& out, const ideal& I)
  {
    I->myOutputSelf(out);
    return out;
  }


  void IdealBase::myOutputSelf(OpenMathOutput& OMOut) const
  {
    OMOut->mySendApplyStart();
    OMOut << OpenMathSymbol("cocoa", "ideal");
    const vector<RingElem>& G = myGens();
    const long NumGens = len(G);
    OMOut << NumGens;                // Number of gens, should be an attribute???
    for (long i=0; i < NumGens; ++i) // To be reconsidered ???
      OMOut << G[i];                 // ???
    OMOut->mySendApplyEnd();
  }


  OpenMathOutput& operator<<(OpenMathOutput& OMOut, const ideal& I)
  {
    OMOut->mySendApplyStart();
    OMOut << OpenMathSymbol("cocoa", "ideal");
    OMOut << RingOf(I);
    I->myOutputSelf(OMOut);
    OMOut->mySendApplyEnd();
    return OMOut;
  }


  bool IsElem(ConstRefRingElem r, const ideal& I)
  {
    if (owner(r) != RingOf(I)) CoCoA_ERROR(ERR::MixedRings, "IsElem(r, I)");
    return I->IhaveElem(raw(r));
  }


  bool IsElem(const BigInt& r, const ideal& I) {return IsElem(RingElem(RingOf(I),r), I);}
  bool IsElem(const BigRat& r, const ideal& I) {return IsElem(RingElem(RingOf(I),r), I);}
  bool IsElem(const MachineInt& r, const ideal& I)  {return IsElem(RingElem(RingOf(I),r), I);}


  void IdealBase::myElim(const std::vector<RingElem>& /*v*/)
  { // default implementation
    CoCoA_ERROR(ERR::NYI, "myElim (only for SparsePolyRing)");
  }

  void IdealBase::mySaturate(const ideal& /*v*/)
  { // default implementation
    CoCoA_ERROR(ERR::NYI, "myElim (only for SparsePolyRing)");
  }

  void IdealBase::myMinimalize()
  { // default implementation
    CoCoA_ERROR(ERR::NYI, "myMinimalize (only for SparsePolyRing)");
  }

} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/ideal.C,v 1.29 2014/07/31 14:45:19 abbott Exp $
// $Log: ideal.C,v $
// Revision 1.29  2014/07/31 14:45:19  abbott
// Summary: Merged io.H and UtilsTemplate.H into new header VectorOperations.H
// Author: JAA
//
// Revision 1.28  2014/07/30 14:12:29  abbott
// Summary: Changed name AmbientRing --> RingOf
// Author: JAA
//
// Revision 1.27  2014/07/14 15:09:50  abbott
// Summary: Changed include of tmp.H into UtilsTemplate.H
// Author: JAA
//
// Revision 1.26  2014/04/30 16:27:17  abbott
// Summary: Replaced size_t by long; Replaced X.size() by len(X)
// Author: JAA
//
// Revision 1.25  2014/03/27 17:15:28  abbott
// Summary: Added products ideal*RingElem and RingElem*ideal
// Author: JAA
//
// Revision 1.24  2014/03/27 14:57:23  bigatti
// -- added myMinimalize
//
// Revision 1.23  2014/03/21 13:08:46  bigatti
// -- added check in intersect
//
// Revision 1.22  2013/06/03 10:53:05  bigatti
// -- just sorted includes
//
// Revision 1.21  2013/02/21 17:18:05  bigatti
// -- now using ERR::Empty in ideal ctor
//
// Revision 1.20  2012/05/30 16:04:55  bigatti
// -- applied "3" convention on bool3 functions and member fields
//
// Revision 1.19  2012/05/29 07:45:23  abbott
// Implemented simplification change to bool3:
//  changed names of the constants,
//  changed names of the testing fns.
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
// Revision 1.15  2011/05/24 14:52:39  abbott
// Removed several "old style" pseudo ctors for principal ideals.
// Cleaned defn of power.
//
// Revision 1.14  2011/05/23 12:35:49  bigatti
// -- added power(ideal, long/ZZ)
//
// Revision 1.13  2011/03/11 10:54:56  bigatti
// -- added mySaturate
// -- added IsElem for integer values
//
// Revision 1.12  2011/01/28 17:58:07  bigatti
// -- added myElim
//
// Revision 1.11  2010/07/16 13:10:34  bigatti
// -- MakeUnique is now exception safe
//
// Revision 1.10  2010/06/10 08:00:02  bigatti
// -- fixed naming conventions
//
// Revision 1.9  2010/03/18 13:55:25  abbott
// Added pseudo-ctors for principal ideals from BigRats.
//
// Revision 1.8  2010/02/04 09:57:11  bigatti
// -- added "mul" for ideals.  Implemented only for SparsePolyRing
//
// Revision 1.7  2009/12/03 17:35:42  abbott
// Moved IsElem fn from .H to .C file, so that .H does not need
// to include CoCoA/error.H.
//
// Revision 1.6  2009/07/30 15:36:28  bigatti
// -- added convenience constructor for ideals with 2, 3, 4 generators.
//
// Revision 1.5  2009/07/24 15:12:33  bigatti
// -- new constructors: ideal(r), ideal(gens)
//
// Revision 1.4  2008/12/17 12:11:52  abbott
// Changed type from long to MachineInt in operations which use a machine integer
// in place of a RingElem.  The change is "superficial" but affects many files.
//
// Revision 1.3  2007/12/07 15:27:01  bigatti
// -- default implementation of "IamOne" in ideal.C
//
// Revision 1.2  2007/10/30 17:14:06  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.3  2007/03/08 18:22:28  cocoa
// Just whitespace cleaning.
//
// Revision 1.2  2006/11/02 13:25:43  cocoa
// Simplification of header files: the OpenMath classes have been renamed.
// Many minor consequential changes.
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.5  2006/05/12 16:10:58  cocoa
// Added OpenMathFwd.H, and tidied OpenMath.H.
// Many consequential but trivial changes.
//
// Revision 1.4  2006/03/21 09:43:13  cocoa
// Changed names of some member fns of ideals (dealing with setting and testing
// the flags for primeness and maximality).  Hope icc will complain less now.
//
// Revision 1.3  2006/03/15 18:09:31  cocoa
// Changed names of member functions which print out their object
// into myOutputSelf -- hope this will appease the Intel C++ compiler.
//
// Revision 1.2  2006/03/12 21:28:33  cocoa
// Major check in after many changes
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.1.1.1  2005/05/03 15:47:31  cocoa
// Imported files
//
// Revision 1.5  2005/04/29 15:42:02  cocoa
// Improved documentation for GMPAllocator.
// Added example program for GMPAllocator.
// Added example program for simple ops on polynomials.
// Added two new ctors for (principal) ideals (from long, and from ZZ).
// Added (crude) printing for PPMonoids.
// Updated library.H (#included GMPAllocator.H).
//
// Revision 1.4  2005/04/20 15:40:47  cocoa
// Major change: modified the standard way errors are to be signalled
// (now via a macro which records filename and line number).  Updated
// documentation in error.txt accordingly.
//
// Improved the documentation in matrix.txt (still more work to be done).
//
// Revision 1.3  2005/04/19 14:06:03  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.2  2005/02/11 14:15:20  cocoa
// New style ring elements and references to ring elements;
// I hope I have finally got it right!
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.9  2004/11/12 15:49:29  cocoa
// Tidying prior to 0.90 release.
// (a) documentation improved (or marked as poor)
// (b) sundry minor improvements to the code
//
// Revision 1.8  2004/11/05 15:30:57  cocoa
// Separated FieldIdealImpl from the "abstract" classes ideal and
// IdealBase.
//
// Revision 1.7  2004/07/20 15:04:06  cocoa
// The next step in the new "ring element" conversion process:
// handling the case of creating a "const RefRingElem" object
// (since C++ refuses to do this properly itself).
//
// Revision 1.6  2004/06/29 17:10:22  cocoa
// Partially tidied use of "protected" and "private" in various
// base classes.  Checking in at the end of the day -- it works,
// and I wouldn't want it to be lost next time point's disk
// misbehaves.
//
// Revision 1.5  2004/05/27 16:14:02  cocoa
// Minor revision for new coding conventions.
//
// Revision 1.4  2004/03/20 17:46:10  cocoa
// Check in prior to departure to RWCA
//
// Revision 1.3  2003/10/17 10:51:06  cocoa
// Major cleaning, and new naming convention.
//
// Revision 1.2  2003/10/09 12:16:38  cocoa
// New coding convention for rings.
//
// Revision 1.3  2003/06/23 17:00:34  abbott
// Minor cleaning prior to public release.
//
// Revision 1.2  2003/05/30 12:00:22  abbott
// Correct ctor for principal ideal (was missing const).
//
// Revision 1.1  2003/04/24 14:43:29  abbott
// Initial revision
//
//
