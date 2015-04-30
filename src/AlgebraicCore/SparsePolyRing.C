//   Copyright (c)  2005,2009  John Abbott

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


// Source code for abstract class SparsePolyRing and friends

#include "CoCoA/SparsePolyRing.H"

#include "CoCoA/BigRat.H"
#include "CoCoA/FGModule.H"  // for myGcd
#include "CoCoA/IntOperations.H"
#include "CoCoA/SmallFpImpl.H"
#include "CoCoA/DUPFp.H"
#include "CoCoA/OpenMath.H"
#include "CoCoA/TeX.H"
#include "CoCoA/PPMonoidHom.H"
#include "CoCoA/CanonicalHom.H"
#include "CoCoA/QuotientRing.H" // for IsQuotientRing
#include "CoCoA/ReductionCog.H"
#include "CoCoA/RingDistrMPolyClean.H" // for NewPolyRing_DMP
#include "CoCoA/RingFp.H" // for IsRingFp
#include "CoCoA/RingQQ.H" // for IsQQ
#include "CoCoA/RingTwinFloat.H" // for IsRingTwinFloat
#include "CoCoA/RingZZ.H" // for IsZZ
#include "CoCoA/factor.H"  // for myGcd
#include "CoCoA/TmpGOperations.H"  // for myIntersect, my Elim..
#include "CoCoA/TmpMonomialIdeal.H"
#include "CoCoA/assert.H"
#include "CoCoA/convert.H"
#include "CoCoA/degree.H"
#include "CoCoA/error.H"
#include "CoCoA/geobucket.H" // for myMul
#include "CoCoA/ideal.H"     // for myGcd
#include "CoCoA/matrix.H" // for OrdMat
#include "CoCoA/module.H"    // for myGcd
#include "CoCoA/submodule.H"  // for myGcd
#include "CoCoA/symbol.H"


#include <algorithm>
using std::max;     // for MaxExponent, StdDeg
using std::remove;  // for myColon
using std::sort;    // for AreGoodIndetNames
#include <functional>
using std::not1;    // for AreLPPSquareFree
using std::ptr_fun; // for AreLPPSquareFree
#include <iostream>
// using std::ostream in SparsePolyRingBase::myOutput
#include <iterator>
using std::back_inserter;
#include <list>
#include <map>
using std::map;
// using std::list;
#include <utility>
using std::make_pair;
using std::pair;
//#include <vector>
using std::vector;

namespace CoCoA
{

  namespace // anonymous
  {
    // This fn is needed in a call to std::transform
    CoeffPP CoeffPPCtor(const pair<PPMonoidElem, RingElem>& arg)
    {
      return CoeffPP(arg.second, arg.first);
    }
  } // end of namespace anonymous


  SparsePolyRing NewPolyRing(const ring& CoeffRing, const PPMonoid& PPM)
  {
    return NewPolyRing_DMP(CoeffRing, PPM);
  }

  SparsePolyRing NewPolyRing(const ring& CoeffRing, const std::vector<symbol>& IndetSyms, const PPOrdering& ord)
  {
    return NewPolyRing_DMP(CoeffRing, IndetSyms, ord);
  }

  SparsePolyRing NewPolyRing(const ring& CoeffRing, const std::vector<symbol>& IndetSyms, const PPOrderingCtor& ord)
  {
    return NewPolyRing_DMP(CoeffRing, IndetSyms, ord);
  }

  SparsePolyRing NewPolyRing(const ring& CoeffRing, const std::vector<symbol>& IndetSyms)
  {
    return NewPolyRing_DMP(CoeffRing, IndetSyms);
  }

  SparsePolyRing NewPolyRing(const ring& CoeffRing, long NumIndets)
  {
    return NewPolyRing_DMP(CoeffRing, NumIndets);
  }

  SparsePolyRing NewPolyRing(const ring& CoeffRing, long NumIndets, const PPOrderingCtor& ord)
  {
    return NewPolyRing_DMP(CoeffRing, NumIndets, ord);
  }

  bool AreGoodIndetNames(const ring& CoeffRing, const std::vector<symbol>& IndetNames)
  {
    // inefficient: we might know that SymbolsCR and IndetNames are good
    // does it matter?
    vector<symbol> syms = symbols(CoeffRing);
    syms.insert(syms.end(), IndetNames.begin(), IndetNames.end());
    sort(syms.begin(), syms.end());
    const long NumSyms = len(syms);
    for (long i=0; i < NumSyms-1; ++i)
    {
      if (syms[i] == syms[i+1]) return false;
      if (head(syms[i]) == head(syms[i+1]) &&
          NumSubscripts(syms[i]) != NumSubscripts(syms[i+1]))
        return false;
    }
    return true;
  }


  matrix OrdMat(const SparsePolyRing& Rx)
  { return OrdMat(PPM(Rx)); }

  matrix GradingMat(const SparsePolyRing& Rx)
  { return GradingMat(PPM(Rx)); }


  //---- Functions for creating/building polynomials

  RingElem monomial(const SparsePolyRing& P, ConstRefRingElem c, ConstRefPPMonoidElem pp)
  {
    if (owner(c) != CoeffRing(P)) CoCoA_ERROR(ERR::MixedCoeffRings, "monomial(P,c,pp)");
    if (owner(pp) != PPM(P)) CoCoA_ERROR(ERR::MixedPPMs, "monomial(P,c,pp)");
    if (IsZero(c)) return zero(P);
    return P->myMonomial(raw(c), raw(pp));
  }

  RingElem monomial(const SparsePolyRing& P, const MachineInt& n, ConstRefPPMonoidElem pp)
  {
    return monomial(P, RingElem(CoeffRing(P), n), pp);
  }

  RingElem monomial(const SparsePolyRing& P, const BigInt& N, ConstRefPPMonoidElem pp)
  {
    return monomial(P, RingElem(CoeffRing(P), N), pp);
  }

  RingElem monomial(const SparsePolyRing& P, const BigRat& Q, ConstRefPPMonoidElem pp)
  {
    return monomial(P, RingElem(CoeffRing(P), Q), pp);
  }

  RingElem monomial(const SparsePolyRing& P, ConstRefRingElem c, const std::vector<long>& expv)
  {
    return monomial(P, c, PPMonoidElem(PPM(P), expv));
  }

  RingElem monomial(const SparsePolyRing& P, const MachineInt& n, const std::vector<long>& expv)
  {
    return monomial(P, RingElem(CoeffRing(P), n), PPMonoidElem(PPM(P), expv));
  }

  RingElem monomial(const SparsePolyRing& P, const BigInt& N, const std::vector<long>& expv)
  {
    return monomial(P, RingElem(CoeffRing(P), N), PPMonoidElem(PPM(P), expv));
  }

  RingElem monomial(const SparsePolyRing& P, const BigRat& Q, const std::vector<long>& expv)
  {
    return monomial(P, RingElem(CoeffRing(P), Q), PPMonoidElem(PPM(P), expv));
  }


  RingElem SparsePolyRingBase::mySymbolValue(const symbol& s) const
  {
    std::vector<symbol> syms = symbols(myPPM());
    syms.push_back(s);
    if (!AreDistinct(syms))
      return myMonomial(raw(one(myCoeffRing())), raw(myPPM()->mySymbolValue(s)));
    if (AreArityConsistent(syms))
      return myCoeffEmbeddingHomCtor()(myCoeffRing()->mySymbolValue(s));
    CoCoA_ERROR(ERR::BadIndetNames, "SparsePolyRingBase::mySymbolValue");
    return myZero(); // just to keep the compiler quiet
  }


  namespace // for functions local to this file/compilation unit.
  {
    inline void CheckCompatible(ConstRefRingElem x, ConstRefRingElem y, const char* const FnName)
    {
      if (owner(x) != owner(y))  CoCoA_ERROR(ERR::MixedRings, FnName);
    }

    inline void CheckElemSparsePolyRing(ConstRefRingElem f, const char* const FnName)
    {
      if (!IsSparsePolyRing(owner(f))) CoCoA_ERROR(ERR::NotElemSparsePolyRing, FnName);
    }

    void CheckCoeffExpv(const SparsePolyRing& P,
                        ConstRefRingElem c, const std::vector<long>& expv,
                        const char* const FnName)
    {
      if (CoeffRing(P) != owner(c))    CoCoA_ERROR(ERR::MixedCoeffRings, FnName);
      if (NumIndets(P) != len(expv))   CoCoA_ERROR(ERR::BadArraySize, FnName);
    }
    
    void CheckCoeffPP(const SparsePolyRing& P,
                      ConstRefRingElem c, ConstRefPPMonoidElem pp,
                      const char* const FnName)
    {
      if (CoeffRing(P) != owner(c))    CoCoA_ERROR(ERR::MixedCoeffRings, FnName);
      if (PPM(P) != owner(pp))         CoCoA_ERROR(ERR::MixedPPMs, FnName);
    }
  }
  
  
  RingElem& PushFront(RingElem& f, ConstRefRingElem c, const std::vector<long>& expv) /// SHOULD BE vector<BigInt> ????
  {
    CheckElemSparsePolyRing(f, "PushFront(f, c, expv)");
    const SparsePolyRing Rx = owner(f);
    CheckCoeffExpv(Rx, c, expv, "PushFront(f, c, expv)");
    PPMonoidElem pp(PPM(Rx), expv);
    if (!IsZero(f) && pp <= LPP(f))
      CoCoA_ERROR(ERR::PPOrder, "PushFront(f, c, expv)");
    Rx->myPushFront(raw(f), raw(c), raw(pp)); // OK 'cos makes a copy of raw(pp)
    return f;
  }


  RingElem& PushFront(RingElem& f, ConstRefRingElem c, ConstRefPPMonoidElem pp)
  {
    CheckElemSparsePolyRing(f, "PushFront(f, c, pp)");
    const SparsePolyRing Rx = owner(f);
    CheckCoeffPP(Rx, c, pp, "PushFront(f, c, pp)");
    if (!IsZero(f) && pp <= LPP(f)) CoCoA_ERROR(ERR::PPOrder, "PushFront(f, c, pp)");
    Rx->myPushFront(raw(f), raw(c), raw(pp));
    return f;
  }


  RingElem& PushBack(RingElem& f, ConstRefRingElem c, const std::vector<long>& expv) /// SHOULD BE vector<BigInt> ????
  {
    CheckElemSparsePolyRing(f, "PushBack(f, c, expv)");
    const SparsePolyRing Rx = owner(f);
    CheckCoeffExpv(Rx, c, expv, "PushBack(f, c, expv)");
    PPMonoidElem pp(PPM(Rx), expv);
    Rx->myPushBack(raw(f), raw(c), raw(pp)); // OK 'cos makes a copy of raw(pp)
    return f;
  }


  RingElem& PushBack(RingElem& f, ConstRefRingElem c, ConstRefPPMonoidElem pp)
  {
    CheckElemSparsePolyRing(f, "PushBack(f, c, pp)");
    const SparsePolyRing Rx = owner(f);
    CheckCoeffPP(Rx, c, pp, "PushBack(f, c, pp)");
    Rx->myPushBack(raw(f), raw(c), raw(pp));
    return f;
  }


  RingElem ClearDenom(const SparsePolyRing& ZZx, const RingElem& f)
  {
    const ring& P = owner(f);
    if (!IsSparsePolyRing(P) ||
        !IsFractionField(CoeffRing(P)) ||
        BaseRing(CoeffRing(P)) != CoeffRing(ZZx) ||
        PPM(P) != PPM(ZZx))
      CoCoA_ERROR(ERR::BadArg, "ClearDenom(NewRing, f)");
    const RingElem D = CommonDenom(f);
    RingElem ans(ZZx);
    for (SparsePolyIter it=BeginIter(f); !IsEnded(it); ++it)
      PushBack(ans, num(coeff(it))*(D/den(coeff(it))), PP(it));
    return ans;
  }



  /*----------------------------------------------------------------------
    Member functions every concrete SparsePolyRing implementation
    must have in addition to those of PolyRingBase.
    ----------------------------------------------------------------------*/

  bool SparsePolyRingBase::myIsValid(ConstRawPtr rawf) const
  {
    if (myIsZero(rawf)) return true;
    SparsePolyIter itf = myBeginIter(rawf);
    if (IsZero(coeff(itf))) return false;
    PPMonoidElem PrevPP = PP(itf);
    for (++itf; !IsEnded(itf); ++itf)
    {
      if (IsZero(coeff(itf))) return false;
      if (PrevPP <= PP(itf)) return false;
      PrevPP = PP(itf);
    }
    return true;
  }


  // ANNA: add check if ordering is StdDeg compatible and return StdDeg(LPP)
  long SparsePolyRingBase::myStdDeg(ConstRawPtr rawf) const
  {
    if (myIsZero(rawf)) CoCoA_ERROR(ERR::ZeroRingElem, "myStdDeg(rawf)");
    long PolyDegree = 0;
    for (SparsePolyIter i=myBeginIter(rawf); !IsEnded(i); ++i)
      PolyDegree = max(PolyDegree, StdDeg(PP(i)));
    return PolyDegree;
  }


  long SparsePolyRingBase::myDeg(ConstRawPtr rawf, long index) const
  {
    if (myIsZero(rawf)) CoCoA_ERROR(ERR::ZeroRingElem, "myDeg(rawf, index)");
    long res = 0;
    for (SparsePolyIter i=myBeginIter(rawf); !IsEnded(i); ++i)
      res = max(res, exponent(PP(i), index));
    return res;
  }


  void SparsePolyRingBase::myContent(RawPtr rawcontent, ConstRawPtr rawf) const
  {
    const ring& R = myCoeffRing();
    CoCoA_ASSERT(IsTrueGCDDomain(R));
    // Compute answer in local var to avoid aliasing problems; also exception clean.
    RingElem ans(R);
    for (SparsePolyIter i=myBeginIter(rawf); !IsEnded(i); ++i)
    {
      R->myGcd(raw(ans), raw(ans), raw(coeff(i))); // ans = GCD(ans, coeff(i));
      if (IsOne(ans)) break;
    }
    // Finally, swap answer into rawcontent.
    R->mySwap(rawcontent, raw(ans));
    return;
  }

  void SparsePolyRingBase::myContentFrF(RawPtr rawcontent, ConstRawPtr rawf) const
  {
    CoCoA_ASSERT(IsFractionFieldOfGCDDomain(myCoeffRing()));
    const FractionField R = myCoeffRing();
    const ring& S = BaseRing(R);
    RingElem N(S);
    RingElem D(S,1);
    for (SparsePolyIter i=myBeginIter(rawf); !IsEnded(i); ++i)
    {
      N = gcd(N, num(coeff(i)));
      D = lcm(D, den(coeff(i)));
//      S->myGcd(raw(ans), raw(ans), raw(num(coeff(i)))); // ans = GCD(ans, num(coeff(i)));
//      if (IsOne(ans)) break;
    }
    RingHom phi = EmbeddingHom(R);
    RingElem ans = phi(N)/phi(D);
    // Finally, swap answer into rawcontent.
    R->mySwap(rawcontent, raw(ans));
  }


  void SparsePolyRingBase::myCommonDenom(RawPtr rawcontent, ConstRawPtr rawf) const
  {
    CoCoA_ASSERT(IsFractionFieldOfGCDDomain(myCoeffRing()));
    const ring& R = BaseRing(myCoeffRing());
    // Compute result in local "ans" to avoid aliasing problems; also exception clean.
    RingElem ans = one(R);
    for (SparsePolyIter it=myBeginIter(rawf); !IsEnded(it); ++it)
    {
      const RingElem D = den(coeff(it));
      ans *= D/gcd(ans, D);
    }
    // Finally, swap answer into rawcontent -- cheaper than assignment.
    R->mySwap(rawcontent, raw(ans));
  }


  void SparsePolyRingBase::myClearDenom(RawPtr rawg, ConstRawPtr rawf) const
  {
    CoCoA_ASSERT(IsFractionFieldOfGCDDomain(myCoeffRing()));
    const ring& R = BaseRing(myCoeffRing());
    RingElem c(R);
    myCommonDenom(raw(c), rawf);
    const RingElem coeff = EmbeddingHom(myCoeffRing())(c);
    RingElem ans = RingElemAlias(ring(this),rawf);
    myMulByCoeff(raw(ans), raw(coeff));
    // Finally, swap answer into rawg -- cheaper than assignment.
    mySwap(rawg, raw(ans));
  }


  void SparsePolyRingBase::myRemoveBigContent(RawPtr rawf) const
  {
    CoCoA_ASSERT(IsTrueGCDDomain(myCoeffRing()));
    CoCoA_ASSERT(!myIsZero(rawf));
    RingElem cont(myCoeffRing());
    myContent(raw(cont), rawf);
    myDivByCoeff(rawf, raw(cont));
  }


  /*----------------------------------------------------------------------
    Member functions inherited from ring with a single implementation
    for all SparsePolyRing implementations
    ----------------------------------------------------------------------*/

  void SparsePolyRingBase::myOutputSelf(OpenMathOutput& OMOut) const
  {
    OMOut->mySendApplyStart();
    OMOut << OpenMathSymbol("polyd", "poly_ring_d");
    OMOut << myCoeffRing();
    OMOut << myNumIndets(); //???? losing the ordering and grading info here!!!
    OMOut->mySendApplyEnd();
  }


  void SparsePolyRingBase::myMul(RawPtr rawlhs, ConstRawPtr rawf, ConstRawPtr rawg) const
  {
    if (myIsZero(rawf) || myIsZero(rawg)) { myAssignZero(rawlhs); return; }
    if (myIsConstant(rawf))
    {
      RingElem ans(RingElemAlias(ring(this), rawg));
      myMulByCoeff(raw(ans), raw(myLC(rawf)));  // weak exc guarantee, but not a problem here
      mySwap(rawlhs, raw(ans));
      return;
    }
    if (myIsConstant(rawg) && myCoeffRing()->IamCommutative())  // CoeffRing should be comm
    {
      myMul(rawlhs, rawg, rawf);
      return;
    }
    const long gLen = myNumTerms(rawg);
    if (IamCommutative() && myNumTerms(rawf) > gLen) { myMul(rawlhs, rawg, rawf); return; }
    const SparsePolyRing P(this);
    RingElemAlias g(P, rawg);

    if (myIsMonomial(rawf))
    {
      RingElem ans(P);
      myAddMul(raw(ans), rawf, rawg);
      mySwap(raw(ans), rawlhs);
      return;
    }
    geobucket gbk(P);
    for (SparsePolyIter itf=myBeginIter(rawf); !IsEnded(itf); ++itf)
      gbk.myAddMul(monomial(P, coeff(itf), PP(itf)), g, gLen);
    RingElem ans(P);
    AddClear(ans, gbk); // this is for exception safety
    mySwap(raw(ans), rawlhs); // really an assignment
  }


  bool SparsePolyRingBase::myIsDivisible(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    if (myIsZero(rawy)) { return false; }
    if (myIsZero(rawx)) { myAssignZero(rawlhs); return true; }
    if (myIsConstant(rawy))
    {
      RingElem ans(RingElemAlias(ring(this), rawx));
      if (!myDivByCoeff(raw(ans), raw(myLC(rawy)))) // if myDivByCoeff throws then ans may be meaningless but it's a local variable
        return false;
      mySwap(rawlhs, raw(ans));
      return true;
    }
    if (!IamCommutative()) { CoCoA_ERROR(ERR::NYI, "SparsePolyRingBase::myDiv non commutative"); }
    const SparsePolyRing P(this);
    RingElem xCopy(RingElemAlias(P, rawx));
    geobucket gbk(P);
    gbk.myAddClear(xCopy, NumTerms(xCopy));
    RingElemAlias y(P, rawy);
    const long yLen = NumTerms(y);
    RingElem ans(P);
    const ring& R = myCoeffRing();
    RingElem coeff(R);
    // if not divisible will throw ERR::BadQuot
    while ( !IsZero(gbk) )
    {
      if (!R->myIsDivisible(raw(coeff), raw(LC(gbk)), raw(myLC(rawy)))) return false;
      if (!IsDivisible(LPP(gbk), LPP(y))) return false;      
      RingElem m(monomial(P, coeff, LPP(gbk)/LPP(y)));
      // xCopy -= m*y;
      gbk.myAddMul(-m, y, yLen);
      // ans += m;
      myAppendClear(raw(ans), raw(m));
    }
    mySwap(raw(ans), rawlhs); // really an assignment
    return true;
  }


  void SparsePolyRingBase::myDeriv(RawPtr rawlhs, ConstRawPtr rawf, ConstRawPtr rawx) const
  {
    if (myIsOne(rawx)) { myAssign(rawlhs, rawf); return; }
    const long n = myNumIndets();
    const SparsePolyRing P(this);
    vector<long> expv(n);
    exponents(expv, myLPP(rawx));

    /* Differentiating x^n -> n*x^(n-1) requires multiplying monomials
     * by exponents.
     *
     * Two cases to consider: the exponents are integers, or not
     * (generic RingElems).
     *
     * The cases are handled separately due to the difficulty of
     * multiplying RingElems from different rings, though we can
     * always multiply by integers without explicit homomorphisms.
     * For non-integers, CanonicalHom is used to inject exponents
     * into the polynomial ring.
     */

    RingElem ans(P);
    for (SparsePolyIter itf=myBeginIter(rawf); !IsEnded(itf); ++itf)
    {
      RingElem scale(P, 1);

      for (long indet=0; indet < n; ++indet)
        if (expv[indet] != 0)
        {
          const RingElem d = RingElemExponent(PP(itf), indet);
	  BigInt D;
	  long dd;
	  if (IsInteger(D, d) && IsConvertible(dd,D)) {
	    if (dd < expv[indet]) { scale = 0; break; }
	    scale *= RangeFactorial(dd-expv[indet]+1, dd);
	  } else {
	    for (long i=0; i<expv[indet]; i++) {
	      scale *= CanonicalHom(owner(d), P)(d-i);
	    }
	  }
        }
      if (IsZero(scale)) continue;

      RingElem m(scale * monomial(P, coeff(itf), PP(itf)/myLPP(rawx)));
      if (!IsZero(m)) myAppendClear(raw(ans), raw(m));
    }
    mySwap(raw(ans), rawlhs); // really an assignment
  }


  void SparsePolyRingBase::mySymbols(std::vector<symbol>& SymList) const
  {
    myCoeffRing()->mySymbols(SymList);
//    myPPM()->mySymbols(SymList);
    const vector<symbol>& PPs = symbols(myPPM());
    SymList.insert(SymList.end(), PPs.begin(), PPs.end());
  }


//   void SparsePolyRingBase::mySymbolsAndValues(std::vector<symbol>& SymList, std::vector<RingElem>& v) const
//   {
//     std::vector<RingElem>& TmpV;
//     myCoeffRing()->mySymbolsAndValues(SymList, TmpV);
//     myPPM()->mySymbols(SymList);
//     v.clear();
//     for (vector<RingElem>::const_iterator it=TmpV.begin(); it!=TmpV.end() ; ++it)
//       v.push_back(myCoeffEmbeddingHomCtor()(*it));
//     for (vector<PPMonoidElem>::const_iterator it=indets(myPPM()).begin(); it!=indets(myPPM()).end() ; ++it)
//       v.push_back(myMonomial(1, raw(*it)));
//   }


  void SparsePolyRingBase::myOutput(std::ostream& out, ConstRawPtr rawf) const
  {
    if (myIsZero(rawf)) { out << '0'; return; }


    const ring& R = myCoeffRing();
    const PPMonoid PPM = myPPM();
    const bool IsQuotientOfZZ = IsQuotientRing(R) && IsZZ(BaseRing(R));
    //const bool IsQuotientOfZ = IsRingFp(R);
    const bool IsNumberRing = IsZZ(R) || IsQuotientOfZZ || IsRingTwinFloat(R); // || IsQQ(R) 

    bool IsFirstCoeff = true;
    for (SparsePolyIter itf=myBeginIter(rawf); !IsEnded(itf) ; ++itf, IsFirstCoeff = false)
    {
      bool PrintStar = true;
      RingElemAlias c = coeff(itf);
      ConstRefPPMonoidElem pp = PP(itf);
      const bool IsWithMinus = R->myIsPrintedWithMinus(raw(c));

      // ---- coefficient ----
      if (!IsFirstCoeff)
      {
        out << ' ';
        if (!IsWithMinus)  out << '+';
      }
      if (IsOne(pp)) { out << c; continue; }
      // ---------- PP != 1 ----------
      if (IsOne(c))
      { // Do not print "1 * "...
        PrintStar = false;
        goto PrintPP;
      }
      if ( IsWithMinus && IsMinusOne(c) ) // in some Z/(n) "-1" prints as n-1
      { // Do not print "-1 * "...
        out << '-';
        PrintStar = false;
        goto PrintPP;
      }
      // General case: coeff is neither +1 nor -1
      if (IsNumberRing || R->myIsPrintAtom(raw(c)) ||
          (IsWithMinus && R->myIsPrintAtom(raw(-c))) )
      {
        out << c;
        goto PrintPP;
      }
      if (!IsFirstCoeff && IsWithMinus) out << '+'; // not printed before
      out << '(' << c << ')';

    PrintPP:      // ---- PP ----
      if (PrintStar && !TeX_mode(out))  out << '*';
      out << pp;
    }
  }


  bool SparsePolyRingBase::myIsPrintAtom(ConstRawPtr rawx) const
  {
    long NoUse;
    if (myIsIndet(NoUse, rawx)) return true;
    if (myIsMonomial(rawx)) return true;
    if (!myIsConstant(rawx)) return false;
    return myCoeffRing()->myIsPrintAtom(raw(myLC(rawx)));
  }


  bool SparsePolyRingBase::myIsPrintedWithMinus(ConstRawPtr rawx) const
  {
    //    if (IsMinusOne(myLC(rawx)) && IsIndet(myLPP(rawx))) return true;
    //    if (!myIsConstant(rawx)) return false;
    return myCoeffRing()->myIsPrintedWithMinus(raw(myLC(rawx)));
  }


  void SparsePolyRingBase::myOutput(OpenMathOutput& OMOut, ConstRawPtr rawx) const
  {
    OMOut->mySendApplyStart();
    OMOut << OpenMathSymbol("cocoa", "DMPSummands");
    OMOut << myNumTerms(rawx);

    for (SparsePolyIter itf=myBeginIter(rawx); !IsEnded(itf) ; ++itf)
    {
      OMOut << coeff(itf);
      OMOut << PP(itf);
//      R->myOutput(OMOut, it->myCoeff);
//      ordering(PPM)->myOutput(OMOut, it->myOrdv);
    }
    OMOut->mySendApplyEnd();
  }


  bool SparsePolyRingBase::myIsOne(ConstRawPtr rawf) const
  {
    if (!myIsMonomial(rawf)) return false;
    if (!IsOne(myLPP(rawf))) return false;
    return IsOne(myLC(rawf));
  }


  bool SparsePolyRingBase::myIsMinusOne(ConstRawPtr rawf) const
  {
    if (!myIsMonomial(rawf)) return false;
    if (!IsOne(myLPP(rawf))) return false;
    return IsMinusOne(myLC(rawf));
  }


  bool SparsePolyRingBase::myIsConstant(ConstRawPtr rawf) const
  {
    if (myIsZero(rawf)) return true;
    if (!myIsMonomial(rawf)) return false;
    return (IsOne(myLPP(rawf)));
  }


  bool SparsePolyRingBase::myIsIndet(long& IndetIndex, ConstRawPtr rawf) const
  {
    if (!myIsMonomial(rawf)) return false;
    if (!IsOne(myLC(rawf))) return false;
    return IsIndet(IndetIndex, myLPP(rawf));
  }


  bool SparsePolyRingBase::myIsIndetPosPower(ConstRawPtr rawf) const
  {
    if (!myIsMonomial(rawf)) return false;
    if (!IsOne(myLC(rawf))) return false;
    return IsIndetPosPower(myLPP(rawf));
  }


  bool SparsePolyRingBase::myIsHomogPartial(ConstRawPtr rawf, long n) const  // assumes 0 <= n <= GrDim
  {
    CoCoA_ASSERT(0 <= n && n <= myGradingDim());
    if (myIsZero(rawf)) { return true; }
    SparsePolyIter itf=myBeginIter(rawf);
    const PPMonoidElem FirstPP=PP(itf);
    for (++itf; !IsEnded(itf); ++itf)
    {
      CoCoA_ASSERT( cmp(FirstPP, PP(itf)) > 0 ); // assert f is correctly sorted
      if ( CmpWDegPartial(FirstPP, PP(itf), n) != 0 )  return false;
    }
    return true;
  }

  bool SparsePolyRingBase::myIsHomog(ConstRawPtr rawf) const
  {
    if (myGradingDim()==0)
      CoCoA_ERROR(ERR::ZeroGradingDim, "SparsePolyRingBase::myIsHomog");
    if (myIsZero(rawf)) { return true; }
    SparsePolyIter itf=myBeginIter(rawf);
    const PPMonoidElem FirstPP=PP(itf);
    for (++itf; !IsEnded(itf); ++itf)
    {
      CoCoA_ASSERT( cmp(FirstPP, PP(itf))>0 ); // assert f is correctly sorted
      if ( CmpWDeg(FirstPP, PP(itf))!=0 )  return false;
    }
    return true;
  }


  void SparsePolyRingBase::myHomog(RawPtr rawfHom, ConstRawPtr rawf, ConstRawPtr rawh) const
  {
    const SparsePolyRing P(this);
    RingElemAlias h(P,rawh);
    CoCoA_ASSERT( myGradingDim()==1 );
    CoCoA_ASSERT( IsIndet(h) );
    CoCoA_ASSERT( IsOne(wdeg(h)[0]) );
    geobucket gbk(P);  // accumulate result into a geobucket for speed
    RingElemAlias f(P,rawf);
    ConstRefPPMonoidElem PPh = LPP(h);
    const BigInt d = wdeg(f)[0];
    for (SparsePolyIter it=BeginIter(f) ; !IsEnded(it) ; ++it )
    {
      RingElem term = monomial(P, coeff(it),PP(it) * power(PPh,d-wdeg(PP(it))[0]));
      gbk.myAddClear(term,1);
    }
    RingElem tmp(P);  // for exception safety
    AddClear(tmp, gbk);
    mySwap(rawfHom, raw(tmp));
  }


  bool SparsePolyRingBase::myIsInvertible(ConstRawPtr rawx) const
  {
    return !myIsZero(rawx) && myIsConstant(rawx) && IsInvertible(myLC(rawx));
  }


  // code for R[x,y]: compute gcd in FractionField(R)[x,y]
  void SparsePolyRingBase::myGcd(RawPtr rawlhs, ConstRawPtr rawf, ConstRawPtr rawg) const
  {
    CoCoA_ASSERT(myNumIndets() != 0);
    if ( myIsInvertible(rawf) || myIsInvertible(rawg) ) { myAssign(rawlhs, raw(myOne())); return; }
    if ( myIsZero(rawf) ) { myAssign(rawlhs, rawg); return; }
    if ( myIsZero(rawg) ) { myAssign(rawlhs, rawf); return; }
    const SparsePolyRing P(this);
    RingElemAlias f(P, rawf);
    RingElemAlias g(P, rawg);
    if (IsZZ(myCoeffRing()) || IsQQ(myCoeffRing()))
    {
      RingElem ans = GCD_DMPZ(f, g);
      mySwap(rawlhs, raw(ans));
      return;
    }
    if (!IsField(myCoeffRing()))
    {
      if (!IsTrueGCDDomain(myCoeffRing()))
        CoCoA_ERROR("NYI gcd of poly with coeffs not in field or TrueGCDDomain", "myGcd");
      else
      {
        FractionField K = NewFractionField(myCoeffRing());
        SparsePolyRing Kx = NewPolyRing_DMP(K, PPM(P));
        RingHom phi = PolyRingHom(P, Kx, CoeffEmbeddingHom(Kx)(EmbeddingHom(K)), indets(Kx));
        RingElem h = gcd(phi(f), phi(g));
        Kx->myClearDenom(raw(h), raw(h));
        RingElem ans(P);
        for (SparsePolyIter it=BeginIter(h) ; !IsEnded(it) ; ++it )
          ans += monomial(P, num(coeff(it)), PP(it));
        myDivByCoeff(raw(ans), raw(content(ans)));
        myMulByCoeff(raw(ans), raw(gcd(content(f), content(g))));
        P->mySwap(rawlhs, raw(ans));
        return;
      }
    }
    // From here onwards myCoeffRing is a *FIELD*
    if (myNumIndets() == 1)
    {
      if (characteristic(myCoeffRing()) > 0 && characteristic(myCoeffRing()) <= SmallFpImpl::ourMaxModulus())
      {
        SmallFpImpl ModP(ConvertTo<long>(characteristic(myCoeffRing())));
        RingElem ans = ConvertFromDUPFp(myIndets()[0], gcd(ConvertToDUPFp(ModP, f), ConvertToDUPFp(ModP, g)));
//        RingElem ans = GCD_DUPFF(f,g);
        P->mySwap(rawlhs, raw(ans));
        return;
      }
      const ideal I = ideal(f, g);
      const vector<RingElem>& GB = GBasis(I);
      if (len(GB) != 1)
        CoCoA_ERROR("Unable to compute GCD", "SparsePolyRingBase::gcd");
      myAssign(rawlhs, raw(GB[0]));
      return;
    }
    const vector<ModuleElem> v = gens(SyzOfGens(ideal(f,g)));
    if ( len(v) != 1 ) 
      CoCoA_ERROR("Unable to compute GCD", "SparsePolyRingBase::gcd");
    RingElem ans = f/((v[0])[1]);
//    if (IsInvertible(ans))  myAssign(rawlhs, raw(myOne()));
//    else P->mySwap(rawlhs, raw(ans));  // exception safe
    ans = monic(ans);
    P->mySwap(rawlhs, raw(ans));  // exception safe
  }


  void SparsePolyRingBase::myNormalizeFracNoGcd(RawPtr rawnum, RawPtr rawden) const
  {
    CoCoA_ASSERT(!myIsZero(rawden));
    // Handle case of 0 specially; later code fails otherwise.
    if (myIsZero(rawnum)) { myAssign(rawden, 1); return; }

    const ring& k = myCoeffRing();
    if (IsTrueGCDDomain(k))
    {
      // Coeff ring is a (true) GCD domain
      if (!IsOrderedDomain(k)) return;
      if (myLC(rawden) > 0) return;
      myNegate(rawnum, rawnum);
      myNegate(rawden, rawden);
      return;
    }
    if (IsFractionFieldOfGCDDomain(k))
    {
      // Coeff ring is a FractionField
      const ring R = BaseRing(k);
      const RingHom embed = EmbeddingHom(k);
      RingElemAlias N(ring(this), rawnum);
      RingElemAlias D(ring(this), rawden);

      const RingElem ContN = content(N);
      const RingElem ContD = content(D);
      const RingElem ContQuot = ContN/ContD;
      myMulByCoeff(rawnum, raw(embed(num(ContQuot))/ContN));
      myMulByCoeff(rawden, raw(embed(den(ContQuot))/ContD));

//       RingElem scale(k);
//       scale = embed(lcm(CommonDenom(N), CommonDenom(D)));
//       if (!IsOne(scale))
//       {
//         myMulByCoeff(rawnum, raw(scale));
//         myMulByCoeff(rawden, raw(scale));
//       }
//       scale = embed(gcd(content(N), content(D)));
//       if (!IsOne(scale))
//       {
//         myDivByCoeff(rawnum, raw(scale));
//         myDivByCoeff(rawden, raw(scale));
//       }
      if (!IsOrderedDomain(k)) return;
      if (myLC(rawden) > 0) return;
      myNegate(rawnum, rawnum);
      myNegate(rawden, rawden);
      return;
    }
    // This must come after the case handling FractionFields!
    if (IsField(k))
    {
      if (IsOne(myLC(rawden))) return;
      myDivByCoeff(rawnum, raw(myLC(rawden))); // not exc safe :-(
      myDivByCoeff(rawden, raw(myLC(rawden)));
      return;
    }
    if (IsOrderedDomain(k))
    {
      if (myLC(rawden) > 0) return;
      myNegate(rawnum, rawnum);
      myNegate(rawden, rawden);
      return;
    }
    CoCoA_ERROR(ERR::NYI, "SparsePolyRing::myNormalizeFracNoGcd");
  }


  bool SparsePolyRingBase::myIsInteger(BigInt& N, ConstRawPtr rawx) const
  {
    if (myIsZero(rawx)) { N = 0; return true; }
    if (myNumTerms(rawx) > 1) return false;
    if (!IsOne(myLPP(rawx))) return false;
    return IsInteger(N, myLC(rawx));
  }


  bool SparsePolyRingBase::myIsRational(BigRat& Q, ConstRawPtr rawx) const
  {
    if (myIsZero(rawx)) { Q = 0; return true; }
    if (myNumTerms(rawx) > 1) return false;
    if (!IsOne(myLPP(rawx))) return false;
    return IsRational(Q, myLC(rawx));
  }


  void SparsePolyRingBase::mySquare(RawPtr rawlhs, ConstRawPtr rawx) const
  {
    if (IamCommutative())
    {
      const long NT = myNumTerms(rawx);
      if (NT > 4)
      {
        RingElem f(ring(this));
        RingElem g(ring(this));
        long count=0;
        for (SparsePolyIter it=myBeginIter(rawx) ; !IsEnded(it) ; ++it )
        {
          if (++count < NT/2) myPushBack(raw(f), raw(coeff(it)), raw(PP(it)));
          else myPushBack(raw(g), raw(coeff(it)), raw(PP(it)));
        }
        RingElem h = power(f, 2) + power(g, 2) + 2*f*g;
        mySwap(rawlhs, raw(h));
        return;
      }
    }
    mySequentialPower(rawlhs, rawx, 2); //??? BUG/SLUG myBinaryPower better if univariate or coeffs are finite field
  }


  void SparsePolyRingBase::myPowerSmallExp(RawPtr rawlhs, ConstRawPtr rawx, long n) const
  {
    // Assert that we have a genuinely non-trivial case.
    CoCoA_ASSERT(n > 1);
    CoCoA_ASSERT(!myIsZero(rawx) && !myIsOne(rawx) && !myIsMinusOne(rawx));

    // anna 3 apr 2008: under testing
    if (myIsMonomial(rawx) && IamCommutative())
    {
      RingElem m = monomial(SparsePolyRing(this), power(myLC(rawx), n), power(myLPP(rawx), n));
      mySwap(rawlhs, raw(m));
      return;
    }
    if (myIsMonomial(rawx)) { myBinaryPower(rawlhs, rawx, n); return; }
    if (n==2)
    {
      mySquare(rawlhs, rawx);
      return;
    }
    mySequentialPower(rawlhs, rawx, n); //??? BUG/SLUG myBinaryPower better if univariate or coeffs are finite field
  }

  void SparsePolyRingBase::myPowerRingElemExp(RawPtr rawlhs, ConstRawPtr rawx, ConstRefRingElem pow) const
  {
    BigInt N;

    if (IsInteger(N, pow)) {
      long n;
      if (IsConvertible(n, N)) {
	myPowerSmallExp(rawlhs, rawx, n);
      } else {
	CoCoA_ERROR(ERR::ExpTooBig, "SparsePolyRing::power(RingElem, RingElem)");
      }
    } else {
      if (myIsMonomial(rawx) && IamCommutative()) {
	RingElem m = monomial(SparsePolyRing(this), power(myLC(rawx), pow), power(myLPP(rawx), pow));
	mySwap(rawlhs, raw(m));
      } else {
	CoCoA_ERROR(ERR::NYI, "SparsePolyRing::power(RingElem, RingElem) non-monomial case");
      }
    }
  }


  //---- Special functions on RingElem owned by SparsePolyRing

  long UnivariateIndetIndex(ConstRefRingElem f)  // returns < 0 if not univariate
  {
    const SparsePolyRing P = owner(f);
    const long nvars = NumIndets(P);
    if (IsZero(f) || StdDeg(f) == 0 || nvars == 0) return 0; // f is constant, or there is only 1 indet.
    vector<long> expv(nvars);
    exponents(expv, LPP(f));
    long ans = 0;
    while (expv[ans] == 0) ++ans;
    for (SparsePolyIter it=BeginIter(f); !IsEnded(it); ++it)
    {
      exponents(expv, PP(it));
      for (int i=0; i < nvars; ++i)
        if (i != ans && expv[i] != 0) return -1;
    }
    return ans;
  }


  degree wdeg(ConstRefRingElem f)
  {
    CheckElemSparsePolyRing(f, "wdeg(f)");
    if (IsZero(f)) CoCoA_ERROR(ERR::ZeroRingElem, "wdeg(f)");
    return wdeg(LPP(f)); // yes! see [KR] introduction Sec.4.3
  }


  int CmpWDeg(ConstRefRingElem f1, ConstRefRingElem f2)
  {
    CheckCompatible(f1, f2, "CmpWDeg(f1,f2)");
    CheckElemSparsePolyRing(f1, "CmpWDeg(f1,f2)");
    if (IsZero(f1) || IsZero(f2)) CoCoA_ERROR(ERR::ZeroRingElem, "CmpWDeg(f1,f2)");
    // Now we know that both f1 and f2 are non-zero
    //    return SparsePolyRingPtr(owner(f1))->myCmpWDeg(raw(f1), raw(f2));
    return CmpWDeg(LPP(f1), LPP(f2));
  }


  int CmpWDegPartial(ConstRefRingElem f1, ConstRefRingElem f2, long i)  // assumes 0 <= i
  {
    CheckCompatible(f1, f2, "CmpWDegPartial(f1,f2,i)");
    CheckElemSparsePolyRing(f1, "CmpWDegPartial(f1,f2,i)");
    if (IsZero(f1) || IsZero(f2)) CoCoA_ERROR(ERR::ZeroRingElem, "CmpWDegPartial(f1,f2,i)");
    // Now we know that both f1 and f2 are non-zero
    //    return SparsePolyRingPtr(owner(f1))->myCmpWDeg(raw(f1), raw(f2));
    return CmpWDegPartial(LPP(f1), LPP(f2), i);
  }


  bool IsHomog(ConstRefRingElem f)
  {
    CheckElemSparsePolyRing(f, "IsHomog(f)");
    if (GradingDim(owner(f))==0)
      CoCoA_ERROR(ERR::ZeroGradingDim, "IsHomog(RingElem)");
    return SparsePolyRingPtr(owner(f))->myIsHomog(raw(f));
  }


  bool IsHomogPartial(ConstRefRingElem f, long n)  // assumes n >= 0
  {
    CheckElemSparsePolyRing(f, "IsHomogPartial(f,n)");
    return SparsePolyRingPtr(owner(f))->myIsHomogPartial(raw(f), n);
  }


  RingElem homog(ConstRefRingElem f, ConstRefRingElem h)
  {
    const char* const FnName = "homog(RingElem, RingElem)";
    CheckCompatible(f, h, FnName);
    const SparsePolyRing P = owner(f);
    if ( GradingDim(P)!=1 )
      CoCoA_ERROR("GrDim must be 1", FnName);
    if ( !IsIndet(h) )
      CoCoA_ERROR("second arg must be an indeterminate", FnName);
    if ( wdeg(h)[0]!=1 )
      CoCoA_ERROR("degree of hom.indet must be 1", FnName);
    RingElem fHom(P);
    P->myHomog(raw(fHom), raw(f), raw(h));
    return fHom;
  }


  RingElem LF(ConstRefRingElem f)
  {
    if (!IsPolyRing(owner(f)))
      CoCoA_ERROR(ERR::NotElemPolyRing, "LC(f)");
    const SparsePolyRing P = owner(f);
    if ( GradingDim(P)==0 )
      CoCoA_ERROR("GradingDim must be non-0", "LF(f)");
    if (IsZero(f)) CoCoA_ERROR(ERR::ZeroRingElem, "LF(f)");
    RingElem F(P);
    ConstRefPPMonoidElem LPPf(LPP(f));
    for (SparsePolyIter it=BeginIter(f) ; !IsEnded(it) ; ++it )
      if (wdeg(PP(it)) == wdeg(LPPf)) F += monomial(P, coeff(it), PP(it));
    return F;
  }


  //----------------------------------------------------------------------
  //??? the following functions to compute NR will be replaced by GBMill

  int FindReducerIndex(ConstRefPPMonoidElem pp, const vector<RingElem>& v)
  {
    const long nelems = len(v);
    for (long i=0; i < nelems; ++i)
      if (IsDivisible(pp, LPP(v[i])))
        return i;
    return -1;
  }


  inline int FindReducerIndex(const ReductionCog& F, const vector<RingElem>& v)
  {
    if ( IsActiveZero(F) ) return -1;
    return FindReducerIndex(ActiveLPP(F), v);
  }


  void ReduceActiveLM(ReductionCog& F, const vector<RingElem>& v)
  {
    int i;
    while ( (i = FindReducerIndex(F, v) ) != -1)
      F->myReduce(v[i]);
  }


  void reduce(ReductionCog& F, const vector<RingElem>& v)
  {
    ReduceActiveLM(F, v);
    while ( !IsActiveZero(F) )
    {
      F->myMoveToNextLM();
      ReduceActiveLM(F, v);
    }
  }
  //--------------------------------------------

  RingElem NR(ConstRefRingElem f, const vector<RingElem>& v)
  {
    if ( IsZero(f) ) return f;
    RingElem ans(f);
    ReductionCog F = NewRedCogGeobucketField(owner(ans));
    F->myAssignReset(ans);
    reduce(F, v);
    F->myRelease(ans);
    return ans;
  }


  SparsePolyIter BeginIter(ConstRefRingElem f)
  {
    CheckElemSparsePolyRing(f, "BeginIter(f)");
    return SparsePolyRingPtr(owner(f))->myBeginIter(raw(f));
  }


  SparsePolyIter EndIter(ConstRefRingElem f)
  {
    CheckElemSparsePolyRing(f, "EndIter(f)");
    return SparsePolyRingPtr(owner(f))->myEndIter(raw(f));
  }


  std::vector<long> indets(ConstRefRingElem f)
  {
    const SparsePolyRing P = owner(f);
    const long nvars = NumIndets(P);
    long NumIndetsSeen = 0;
    vector<bool> bitmap(nvars);
    for (SparsePolyIter it=BeginIter(f); !IsEnded(it); ++it)
    {
      // Next few lines just do ans |= indets(PP(it));
      const vector<long> v = indets(PP(it));
      const long n = len(v);
      for (long i=0; i < n; ++i)
      {
        if (bitmap[v[i]]) continue;
        bitmap[v[i]] = true;
        ++NumIndetsSeen;
      }
      if (NumIndetsSeen == nvars) break;
    }
    const long NumIndets = count(bitmap.begin(), bitmap.end(), true);
    vector<long> ans; ans.reserve(NumIndets);
    for (long i=0; i < nvars; ++i)
      if (bitmap[i])
        ans.push_back(i);
    return ans;
  }


  namespace // anonymous for file local defns
  {
    class ByDecreasingPP
    {
    public:
      bool operator()(const PPMonoidElem& A, const PPMonoidElem& B) const
        {
          return A > B;
        }
    };
  }

  std::ostream& operator<<(std::ostream& out, const CoeffPP& term)
  {
    out << "[coeff:=" << term.myCoeff << ", PP:=" << term.myPP << "]" ;
    return out;
  }


  std::vector<CoeffPP> CoefficientsWRT(ConstRefRingElem f, const std::vector<long>& indets)
  {
    if (!IsSparsePolyRing(owner(f)))
      CoCoA_ERROR(ERR::NotSparsePolyRing, "CoefficientsWRT(f,indets)");
    const SparsePolyRing P = owner(f);
    for (long i=0; i < len(indets); ++i)
      if (indets[i] < 0 || indets[i] >= NumIndets(P))
        CoCoA_ERROR(ERR::BadIndetIndex, "CoefficientsWRT(f,indets)");

    // Force the sorting criterion in the map
    typedef map<PPMonoidElem, RingElem, ByDecreasingPP> CoeffTable_t;
    CoeffTable_t CoeffTable;
    PPMonoidHom projection = RestrictionHom(PPM(P), indets);
    for (SparsePolyIter it = BeginIter(f); !IsEnded(it); ++it)
    {
      const PPMonoidElem t = projection(PP(it));
      CoeffTable_t::iterator pos = CoeffTable.find(t);
      if (pos == CoeffTable.end())
      {
        CoeffTable.insert(make_pair(t, monomial(P, coeff(it), PP(it)/t)));
        continue;
      }
      pos->second += monomial(P, coeff(it), PP(it)/t);
    }
    vector<CoeffPP> ans; ans.reserve(len(CoeffTable));
    transform(CoeffTable.begin(), CoeffTable.end(), back_inserter(ans), CoeffPPCtor);
    // NOTE: CoeffTable will automatically be in decreasing order by PP!
    // (see 23.2.4/10 in C++11)
    return ans;
  }


  std::vector<CoeffPP> CoefficientsWRT(ConstRefRingElem f, ConstRefRingElem x)
  {
    if (owner(f) != owner(x)) CoCoA_ERROR(ERR::MixedRings, "CoefficientsWRT(f,x)");
    if (!IsSparsePolyRing(owner(f))) CoCoA_ERROR(ERR::NotSparsePolyRing, "CoefficientsWRT(f,x)");
    vector<long> indices(1);
    if (!IsIndet(indices[0], x))
      CoCoA_ERROR(ERR::NotIndet, "CoefficientsWRT(f, x)");
    return CoefficientsWRT(f, indices);
  }


  std::vector<RingElem> CoeffVecWRT(ConstRefRingElem f, ConstRefRingElem x)
  {
    if (owner(f) != owner(x)) CoCoA_ERROR(ERR::MixedRings, "CoeffVecWRT(f,x)");
    if (!IsSparsePolyRing(owner(f))) CoCoA_ERROR(ERR::NotSparsePolyRing, "CoeffVecWRT(f,x)");
    if (IsZero(f)) return vector<RingElem>();
    vector<long> indices(1);
    if (!IsIndet(indices[0], x))
      CoCoA_ERROR(ERR::NotIndet, "CoeffVecWRT(f, x)");
    vector<CoeffPP> CoeffList = CoefficientsWRT(f, x);
    long Degf = 0;
    for (int i=0; i < len(CoeffList); ++i)
      Degf = max(Degf, StdDeg(CoeffList[i].myPP));
    vector<RingElem> ans(1+Degf, zero(owner(f)));
    for (vector<CoeffPP>::iterator it=CoeffList.begin(); it != CoeffList.end(); ++it)
      swap(ans[StdDeg(it->myPP)], it->myCoeff); // swap avoids a wasteful copy!
    return ans;
  }


  RingElem ContentWRT(ConstRefRingElem f, const std::vector<long>& indets)
  {
    if (!IsSparsePolyRing(owner(f)))
      CoCoA_ERROR(ERR::NotSparsePolyRing, "ContentWRT(f,indets)");
    const SparsePolyRing P = owner(f);
    for (long i=0; i < len(indets); ++i)
      if (indets[i] < 0 || indets[i] >= NumIndets(P))
        CoCoA_ERROR(ERR::BadIndetIndex, "ContentWRT(f,indets)");

    const vector<CoeffPP> CoeffList = CoefficientsWRT(f, indets);
    RingElem ans(P);
    const long n = len(CoeffList);
    for (long i=0; i < n; ++i)
      ans = gcd(ans, CoeffList[i].myCoeff);
    return ans;
  }

  RingElem ContentWRT(ConstRefRingElem f, ConstRefRingElem x)
  {
    if (owner(f) != owner(x))
      CoCoA_ERROR(ERR::MixedRings, "ContentWRT(f,x)");
    vector<long> indices(1);
    if (!IsIndet(indices[0], x))
      CoCoA_ERROR(ERR::NotIndet, "ContentWRT(f,x)");
    return ContentWRT(f, indices);
  }



  //-- ideals ----------------------------------------

  bool IsZeroDim(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))
      CoCoA_ERROR(ERR::NotSparsePolyRing, "IsZeroDim(I)");
    if (IsZero(I) || IsOne(I)) return false;
    // Now we know I is non-trivial.
    const SparsePolyRing P = RingOf(I);
    const vector<RingElem>& GB = TidyGens(I);
    const long GBsize = len(GB); // MUST BE A REDUCED GBASIS !!!
    long NumIndetPowers = 0;
    for (long i=0; i < GBsize; ++i)
      if (IsIndetPosPower(LPP(GB[i])))
        ++NumIndetPowers;
    return (NumIndetPowers == NumIndets(P));
  }


  bool IsHomog(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))
      CoCoA_ERROR(ERR::NotSparsePolyRing, "IsHomog(ideal)");
    if (GradingDim(RingOf(I))==0)
      CoCoA_ERROR(ERR::ZeroGradingDim, "IsHomog(ideal)");
    if (IsZero(I) || IsOne(I)) return true;
    // Now we know I is non-trivial.
    const SparsePolyRing P = RingOf(I);
    const vector<RingElem>& GB = TidyGens(I);
    const long GBsize = len(GB); // MUST BE A REDUCED GBASIS !!!
    for (long i=0; i < GBsize; ++i)
      if (!IsHomog(GB[i]))  return false;
    return true;
  }


  const std::vector<RingElem>& GBasis(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))
      CoCoA_ERROR(ERR::NotSparsePolyRing, "GBasis(I)");
    return TidyGens(I);
  }


  const std::vector<RingElem>& MinGens(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))
      CoCoA_ERROR(ERR::NotSparsePolyRing, "MinGens(I)");
    const SparsePolyRingBase::IdealImpl* const ptrI = 
      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    return ptrI->myMinGens();
  }


  std::vector<ideal> PrimaryDecomposition(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))
      CoCoA_ERROR(ERR::NotSparsePolyRing, "PrimaryDecomposition(I)");
    if (!AreGensMonomial(I))
      CoCoA_ERROR(ERR::NYI, "PrimaryDecomposition(I)");
    return PrimaryDecompositionMonId(I);
  }


  ideal LT(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))
      CoCoA_ERROR(ERR::NotSparsePolyRing, "LT(I)");
    std::vector<RingElem> GB = TidyGens(I);
    std::vector<RingElem> v;
    const SparsePolyRing P = RingOf(I);
    for (long i=0; i<len(GB); ++i)
      v.push_back(monomial(P, 1, LPP(GB[i])));
    return ideal(P, v);
  }


  ideal LF(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))
      CoCoA_ERROR(ERR::NotSparsePolyRing, "LF(I)");
    const SparsePolyRing P = RingOf(I);
    if ( GradingDim(P)==0 ) CoCoA_ERROR("GradingDim must be non-0", "LF(I)");
    std::vector<RingElem> GB = TidyGens(I);
    std::vector<RingElem> v;
    for (long i=0; i<len(GB); ++i)  v.push_back(LF(GB[i]));
    return ideal(P, v);
  }


  ideal homog(const ideal& I, ConstRefRingElem x)
  {
    if (!IsSparsePolyRing(RingOf(I)))
      CoCoA_ERROR(ERR::NotSparsePolyRing, "homog(I, x)");
    if (AreGensMonomial(I)) return I;
    if (IsZero(I)) return I;
    std::vector<RingElem> HomogIdealGens;
    std::vector<RingElem> v(1,x);
    ComputeHomogenization(HomogIdealGens, gens(I), v);
    return ideal(HomogIdealGens);
  }


  // Anna: should be friend
  bool AreGensMonomial(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))
      CoCoA_ERROR(ERR::NotSparsePolyRing, "AreGensMonomial(I)");
    const SparsePolyRingBase::IdealImpl* const ptrI = 
      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    return ptrI->IhaveMonomialGens();
  }
  

  // Anna: should be friend
  bool AreGensSquareFreeMonomial(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))
      CoCoA_ERROR(ERR::NotSparsePolyRing, "AreGensSquareFreeMonomial(I)");
    const SparsePolyRingBase::IdealImpl* const ptrI = 
      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    return ptrI->IhaveSquareFreeMonomialGens();
  }
  

  // Anna: should be friend
  void SetGBasisAsGens(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))
      CoCoA_ERROR(ERR::NotSparsePolyRing, "GBasisAsGens(I)");
    const SparsePolyRingBase::IdealImpl* const ptrI = 
      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    ptrI->mySetGBasisAsGens();
  }
  

  namespace  //anonymous  for QuotientBasis
  {

    bool IsDivisible(ConstRefPPMonoidElem pp, const std::list<PPMonoidElem>& ByL)
    {
      for (std::list<PPMonoidElem>::const_iterator i=ByL.begin(); i != ByL.end(); ++i)
        if (IsDivisible(pp, *i)) return true;
      return false;
    }
    
    // Apparently STL bind2nd cannot have (const?) reference arguments.
    // Apparently BOOST bind would work: what about C++-11?
    // This is a hack for STL
    class IsDivHack
    {
    public:
      IsDivHack (const std::list<PPMonoidElem>& L): myL(L) {}
      bool operator () (ConstRefPPMonoidElem pp) {return IsDivisible(pp, myL);}
      
    private:
      const std::list<PPMonoidElem> & myL;
    };


    void QuotientBasisRec(std::vector<PPMonoidElem>& ans, 
                          const std::list<PPMonoidElem>& L, 
                          ConstRefPPMonoidElem prefix, 
                          long idx)
    {
      PPMonoid PPM = owner((L.front()));
      const PPMonoidElem& X = indets(PPM)[idx];
      PPMonoidElem prefixXd(prefix);  // prefix * x[idx]^d
      int MaxDeg = 0;
      if (idx == NumIndets(PPM)-1)
      {
        MaxDeg = exponent((L.front()), idx);
        for (int d=0; d < MaxDeg;  ++d, prefixXd *= X)  ans.push_back(prefixXd);
        return;
      }
      for (std::list<PPMonoidElem>::const_iterator it=L.begin(); it != L.end() ; ++it)
        if (exponent((*it),idx) > MaxDeg)  MaxDeg = exponent((*it),idx);
      std::list<PPMonoidElem> CutOff, tmp;
      PPMonoidElem Xd(PPM);  // x[idx]^d
      for (int d=0; d < MaxDeg; ++d, prefixXd *= X, Xd *= X)
      {
        for (std::list<PPMonoidElem>::const_iterator it=L.begin(); it != L.end(); ++it)
          if (exponent(*it,idx) == d)  tmp.push_back(((*it)/Xd));
        CutOff.remove_if(IsDivHack(tmp));
        CutOff.splice(CutOff.end(), tmp);
        QuotientBasisRec(ans, CutOff, prefixXd, idx+1);
      }
    }
  } // anonymous namespace


  std::vector<PPMonoidElem> QuotientBasis(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))
      CoCoA_ERROR(ERR::NotSparsePolyRing, "QuotientBasis(I)");
    if (!IsZeroDim(I))
      CoCoA_ERROR("ideal must be 0-dimensional", "QuotientBasis");
    vector<RingElem> GB = GBasis(I);
    std::list<PPMonoidElem> LeadingPPs;
    vector<PPMonoidElem> ans;
    for (long i=0; i < len(GB); ++i)  LeadingPPs.push_back(LPP(GB[i]));
    QuotientBasisRec(ans, LeadingPPs, PPMonoidElem(PPM(RingOf(I))), 0);
    return ans;
  }


  //-- IdealImpl ----------------------------------------

  ideal SparsePolyRingBase::myIdealCtor(const std::vector<RingElem>& gens) const
  {
    return ideal(new IdealImpl(SparsePolyRing(this), gens)); //??? ugly ???
  }


  SparsePolyRingBase::IdealImpl::IdealImpl(const SparsePolyRing& P, const std::vector<RingElem>& gens):
      myP(P),
      myGensValue(gens),
      myGBasisIsValid(false)
  {
    // IhaveMonomialGens3Value = uncertain3; // default for bool3
    // IhaveSquareFreeMonomial3Gens = uncertain3; // default for bool3
    if (!IsField(CoeffRing(P)))
      CoCoA_ERROR("NYI ideal of polynomials with coeffs not in a field", "ideal(SparsePolyRing, gens)");//???
  }


  IdealBase* SparsePolyRingBase::IdealImpl::myClone() const
  {
    return new IdealImpl(*this);
  }


  const SparsePolyRing& SparsePolyRingBase::IdealImpl::myRing() const
  {
    return myP;
  }


  bool SparsePolyRingBase::IdealImpl::IamZero() const
  {
    for (long i=0; i<len(myGens()); ++i)
      if (!IsZero(myGens()[i])) return false;
    return true;
  }


  
  namespace // anonymous
  {
    bool IsRadLPP(ConstRefRingElem f)
    { return IsRadical(LPP(f)); }
    
    bool AreLPPSquareFree(const std::vector<RingElem>& v)
    {
      const long n = len(v);
      for (long i=0; i < n; ++i)
	if (!IsRadLPP(v[i])) return false;
      return true;
//   We *DO NOT USE* STL algorithm because std::ptr_fun does not work if the fn has formal params which are of reference type
//       return find_if(v.begin(), v.end(),
//                      not1(ptr_fun(CoCoA::IsRadLPP)))
// 	//                     not1(ptr_fun(static_cast<bool(*)(ConstRefRingElem)>(CoCoA::IsRadLPP))))
//         == v.end(); 
    }

    } // anonymous namespace
  

  bool SparsePolyRingBase::IdealImpl::IhaveMonomialGens() const
  {
    if (IsUncertain3(IhaveMonomialGens3Value))
      IhaveMonomialGens3Value = AreMonomials(myGensValue);
    return IsTrue3(IhaveMonomialGens3Value);
  }


  bool SparsePolyRingBase::IdealImpl::IhaveSquareFreeMonomialGens() const
  {
    if (IsUncertain3(IhaveSquareFreeMonomialGens3Value))
    {
      if (!IhaveMonomialGens()) IhaveSquareFreeMonomialGens3Value = false3;
      else IhaveSquareFreeMonomialGens3Value = AreLPPSquareFree(myGensValue);
    }
    return IsTrue3(IhaveSquareFreeMonomialGens3Value);
  }


  void SparsePolyRingBase::IdealImpl::mySetGBasisAsGens() const
  {
    myGBasisValue = myGensValue;
    myGBasisIsValid = true;
  }


//   //------- functions for monomial ideals: where should they go? ------
//   void MonomialInterreduce(std::vector<RingElem>& res, const std::vector<RingElem>& l)
//   {
//     if (l.empty()) { res.clear(); return; }
//     const SparsePolyRing P = owner(l[0]);
//     PPVector g(PPM(P), NewDivMaskEvenPowers());
//     for (vector<RingElem>::const_iterator it=l.begin(); it!=l.end() ; ++it)
//       PushBack(g, LPP(*it));
//     interreduce(g);
//     convert(res, P, g);
//   }


//   void MonomialIntersection(std::vector<RingElem>& res, const std::vector<RingElem>& l1,  const std::vector<RingElem>& l2)
//   {
//     if (l1.empty()) { res.clear(); return; }
//     if (l2.empty()) { res.clear(); return; }
//     const SparsePolyRing P = owner(l1[0]);
//     PPVector g1(PPM(P), NewDivMaskEvenPowers());
//     for (vector<RingElem>::const_iterator it=l1.begin(); it!=l1.end() ; ++it)
//       PushBack(g1, LPP(*it));
//     PPVector g2(PPM(P), NewDivMaskEvenPowers());
//     for (vector<RingElem>::const_iterator it=l2.begin(); it!=l2.end() ; ++it)
//       PushBack(g2, LPP(*it));
//     PPVector g(PPM(P), NewDivMaskEvenPowers());
//     g.myLcms(g1, g2);
//     convert(res, P, g);
//   }
//   //------- functions for monomial ideals: end ------


  void SparsePolyRingBase::IdealImpl::myMaximalTest() const
  {
    if (NumIndets(myP) == 1 && IsField(CoeffRing(myP)))
    {
      // Simple case: we are in univariate poly ring
      const vector<RingElem>& GB = myGBasis();
      if (len(GB) == 0) { mySetMaximalFlag(false); return; }// zero ideal
      mySetMaximalFlag(IsIrred(GB[0]));
      mySetPrimeFlag(IsTrue3(IamMaximal3Flag));
      return;
    }
    // Harder case: poly rings has at least 2 indets.
    CoCoA_ERROR(ERR::NYI, "SparsePolyRingBase::IdealImpl::myMaximalTest() - multivariate case");//???
  }


  void SparsePolyRingBase::IdealImpl::myPrimeTest() const
  {
    if (NumIndets(myP) == 1 && IsField(CoeffRing(myP)))
    {
      myMaximalTest();
      return;
    }
    CoCoA_ERROR(ERR::NYI, "SparsePolyRingBase::IdealImpl::myPrimeTest()");//???
  }


  void SparsePolyRingBase::IdealImpl::myReduceMod(RingElemRawPtr rawf) const
  {
    //??? very basic default implementation
    RingElem tmp = NR(RingElemAlias(myP, rawf), myGBasis());
    myP->mySwap(rawf, raw(tmp));
  }


  bool SparsePolyRingBase::IdealImpl::IhaveElem(RingElemConstRawPtr rawf) const
  {
    RingElem g = RingElemAlias(myP, rawf);
    myReduceMod(raw(g));
    return IsZero(g);
  }


  const SparsePolyRingBase::IdealImpl* SparsePolyRingBase::IdealImpl::ourGetPtr(const ideal& I)
  {
    return dynamic_cast<const SparsePolyRingBase::IdealImpl*>(I.myIdealPtr());
  }


  void SparsePolyRingBase::IdealImpl::myClearGBasis()
  {
    myGBasisIsValid = false;
    myGBasisValue.clear();
  }


  void SparsePolyRingBase::IdealImpl::myAdd(const ideal& Jin)
  {
    const IdealImpl* const J = ourGetPtr(Jin);
    myGensValue.insert(myGensValue.end(), gens(Jin).begin(), gens(Jin).end());
    if (IsTrue3(IhaveMonomialGens3Value))
      IhaveMonomialGens3Value = J->IhaveMonomialGens3Value;
    if (IsTrue3(IhaveSquareFreeMonomialGens3Value))
      IhaveSquareFreeMonomialGens3Value = J->IhaveSquareFreeMonomialGens3Value;
    myClearGBasis();
  }


  // 20140506  What was this fn for???  Apparently unused, so commenting out.
  // namespace{
  //   void MonomialMulTrick(std::vector<RingElem>& Vout,
  //                         const std::vector<RingElem>& v1,
  //                         const std::vector<RingElem>& v2)
  //   {
  //     RingElem out = sum(v1)*sum(v2);
  //     Vout.clear();
  //     const SparsePolyRing P = owner(out);
  //     RingElem OneR = one(CoeffRing(P));
  //     for (SparsePolyIter it=BeginIter(out) ; !IsEnded(it) ; ++it )
  //       Vout.push_back(monomial(P, OneR, PP(it)));
  //   }
  // }
  
  void SparsePolyRingBase::IdealImpl::myMul(const ideal& Jin)
  {
    if (IhaveMonomialGens() && AreGensMonomial(Jin))
    {
      myMulMonId(Jin);
      return;
    }
    vector<RingElem> tmpV;
    const SparsePolyRingBase::IdealImpl* const J = ourGetPtr(Jin);
    for (vector<RingElem>::const_iterator itI=myGensValue.begin(); itI!=myGensValue.end(); ++itI)
      for (vector<RingElem>::const_iterator itJ=J->myGensValue.begin(); itJ!=J->myGensValue.end(); ++itJ)
        tmpV.push_back((*itI)*(*itJ));
    swap(tmpV, myGensValue);  // ANNA does this make copies???  2010-02-03
    if (characteristic(myRing()) != 0)
    {
      IhaveMonomialGens3Value = uncertain3;
      IhaveSquareFreeMonomialGens3Value = uncertain3;
    }
    myClearGBasis();
    //    IamPrime3Flag = IamMaximal3Flag = uncertain3;
  }


  void SparsePolyRingBase::IdealImpl::myIntersect(const ideal& J)
  {
    if (IamZero()) return;
    CoCoA_ASSERT(!IsZero(J));
    if (IhaveMonomialGens() && AreGensMonomial(J))
    {
      myIntersectMonId(J);
      return;
    }
    IhaveMonomialGens3Value = uncertain3;
    IhaveSquareFreeMonomialGens3Value = uncertain3;
    ComputeIntersection(myGensValue, myGensValue, gens(J));
    myClearGBasis();
  }


  void SparsePolyRingBase::IdealImpl::myColon(const ideal& J)
  {
    if (IsZero(J))
      myGensValue = vector<RingElem>(1, one(myRing()));
    else
    {
      if (IhaveMonomialGens() && AreGensMonomial(J))
      {
        myColonMonId(J);
        return;
      }
      const RingElem Z(zero(myRing()));
      myGensValue.erase(remove(myGensValue.begin(), myGensValue.end(),Z),
                        myGensValue.end());
      ComputeCColon(myGensValue, myGensValue, gens(J));
    }
    IhaveMonomialGens3Value = uncertain3;
    IhaveSquareFreeMonomialGens3Value = uncertain3;
    myClearGBasis();
  }


  void SparsePolyRingBase::IdealImpl::mySaturate(const ideal& J)
  {
    ComputeSSaturation(myGensValue,myGensValue, gens(J));
    IhaveMonomialGens3Value = uncertain3;
    IhaveSquareFreeMonomialGens3Value = uncertain3;
    myClearGBasis();
  }


  void SparsePolyRingBase::IdealImpl::myMinimalize()
  {
    if (GradingDim(myRing())==0 || !IsHomog(myGensValue))
      CoCoA_ERROR("Input is not homogeneous", "myMinimalize");
    myGBasis(); // this sets GBasis and MinGens
    myGensValue = myMinGens(); // if monomial ideal min gens are in GBasis
    if (IsFalse3(IhaveMonomialGens3Value))
      IhaveMonomialGens3Value = uncertain3;
    if (IsFalse3(IhaveSquareFreeMonomialGens3Value))
      IhaveSquareFreeMonomialGens3Value = uncertain3;
  }


  void SparsePolyRingBase::IdealImpl::myElim(const std::vector<RingElem>& ElimIndets)
  {
    //    if (myGBasisIsValid) and elim ordering ...
    if (IhaveMonomialGens())
    {
      myElimMonId(ElimIndets);
      return;
    }
    const RingElem Z(zero(myRing()));
    myGensValue.erase(remove(myGensValue.begin(), myGensValue.end(),Z),
                      myGensValue.end());
    PPMonoidElem ElimIndetsProd(myRing()->myPPM());
    const long n = len(ElimIndets);
    for (long i=0 ; i<n ; ++i)
    {
      if (!IsIndet(ElimIndets[i]))
        CoCoA_ERROR(ERR::NotIndet, "myElim");
      ElimIndetsProd *= LPP(ElimIndets[i]);
    }
    ComputeElim(myGensValue,myGensValue, ElimIndetsProd);
    IhaveMonomialGens3Value = uncertain3;
    IhaveSquareFreeMonomialGens3Value = uncertain3;
    myClearGBasis();
  }


  bool SparsePolyRingBase::IdealImpl::myDivMod(RingElemRawPtr rawlhs, RingElemConstRawPtr rawnum, RingElemConstRawPtr rawden) const
  {
    const SparsePolyRing P = myRing();
    if (IsField(CoeffRing(P)) && P->myIsConstant(rawden))
    {
      P->myDiv(rawlhs, rawnum, rawden);
      return true;
    }
// The following should be a general solution...
///    auto tmp = GenRepr(rawnum, ideal(rawden)+this);
///    return tmp[0];

// !!! DOES NOT COMPILE BECAUSE const!!! CoCoA_ASSERT(IsZeroDim(ideal(this)));
//     vector<PPMonoidElem> QB = QuotientBasis(???);
//     const int LenQB = len(QB);
//     for (int i=0; i < LenQB; ++i)
//     {
//       RingElem tmp = monomial(P, one(CoeffRing(P)), QB[i])*RingElemAlias(P,rawden);

//     }
    CoCoA_ERROR(ERR::NYI, "SparsePolyRingBase::IdealImpl::myDivMod");
    return false; // just to keep compiler quiet!!
  }


  const std::vector<RingElem>& SparsePolyRingBase::IdealImpl::myGens() const
  { return myGensValue; }


  const std::vector<RingElem>& SparsePolyRingBase::IdealImpl::myTidyGens() const
  { return myGBasis(); }


  // ??? copied from ex-QuotientBasis.C
  bool IsDivisible(ConstRefPPMonoidElem pp, const vector<PPMonoidElem>& ByL)
  {
    const long n = len(ByL);
    for (long i=0; i < n; ++i)
      if ( IsDivisible(pp, ByL[i]) ) return true;
    return false;
  }


  const std::vector<RingElem>& SparsePolyRingBase::IdealImpl::myGBasis() const
  {
    if (IhaveMonomialGens()) return myGBasisMonId();
    if (myGBasisIsValid) return myGBasisValue;
    CoCoA_ASSERT(myGBasisValue.empty());
    if (IamZero()) return myGBasisValue;
    vector<RingElem> MinGens;
    ComputeGBasis(myGBasisValue, MinGens, myGensValue);
    if (!MinGens.empty())  myMinGensValue = MinGens;
    myGBasisIsValid = true;
    return myGBasisValue;
  }


  const std::vector<RingElem>& SparsePolyRingBase::IdealImpl::myMinGens() const
  {
    if (IhaveMonomialGens()) return myGBasisMonId(); // interreduced
    if (!myMinGensValue.empty()) return myMinGensValue;
    if (GradingDim(myRing())==0 || !IsHomog(myGensValue))
      CoCoA_ERROR("Input is not homogeneous", "myMinGens");
    myGBasis();
    return myMinGensValue;
  }


  //-- HomImpl ----------------------------------------

  SparsePolyRingBase::HomImpl::HomImpl(const SparsePolyRing& domain, const ring& codomain, const RingHom& CoeffHom, const vector<RingElem>& IndetImages):
      RingHomBase(domain, codomain),
      myCoeffHom(CoeffHom),
      myIndetImages(IndetImages)
  {
    // No need to check anything: checks already made when CoeffHom was built.
  }

namespace
{
  // ??? appropriate use of inheritance here?  this is getting pretty hugly

  // assume image==0
  void ApplySPRCodomain(RingElem& image, ConstRefRingElem arg, const RingHom CoeffHom, const vector<RingElem>& IndetImages)
  {
    const SparsePolyRing S = owner(image);
    geobucket gbk(S);

    const long NumInd = len(IndetImages);
    for (SparsePolyIter i=BeginIter(arg); !IsEnded(i); ++i)
    {
      RingElem SummandImage = CoeffHom(coeff(i));
      CoCoA_ASSERT(owner(SummandImage) == S);
      if (IsZero(SummandImage)) continue; // efficiency hack????
      ConstRefPPMonoidElem t(PP(i));
      for (long ind=0; ind < NumInd; ++ind)
      {
        const long d = exponent(t, ind); // ??? should we compute exponents?
        if (d == 0) continue;
        SummandImage *= power(IndetImages[ind], d);
      }
      //        SparsePolyRingPtr(S)->myAddClear(raw(ans), raw(SummandImage));
      gbk.myAddClear(SummandImage, NumTerms(SummandImage));
    }
    AddClear(image, gbk);
  }


  void ApplyGeneral(RingElem& image, ConstRefRingElem arg, const RingHom CoeffHom, const vector<RingElem>& IndetImages)
  {
    ring S = owner(image);
    const long NumInd = len(IndetImages);
    for (SparsePolyIter i=BeginIter(arg); !IsEnded(i); ++i)
    {
      RingElem SummandImage = CoeffHom(coeff(i));
      CoCoA_ASSERT(owner(SummandImage) == S);
      if (IsZero(SummandImage)) continue; // efficiency hack????
      ConstRefPPMonoidElem t(PP(i));
      for (long ind=0; ind < NumInd; ++ind)
      {
        const long d = exponent(t, ind); // ??? should we compute exponents?
        if (d == 0) continue;
        SummandImage *= power(IndetImages[ind], d);
      }
      S->myAdd(raw(image), raw(image), raw(SummandImage));
    }
  }
}  // end of anonymous namespace

  void SparsePolyRingBase::HomImpl::myApply(RingElemRawPtr rawimage, RingElemConstRawPtr rawarg) const
  {
    RingElem ans(myCodomain);  // Putting result into ans is exception safe and avoids aliasing problems.
    if ( IsSparsePolyRing(myCodomain) )
      ApplySPRCodomain(ans, RingElemAlias(myDomain, rawarg), myCoeffHom, myIndetImages);
    else
      ApplyGeneral(ans, RingElemAlias(myDomain, rawarg), myCoeffHom, myIndetImages);
    myCodomain->mySwap(rawimage, raw(ans));
  }


  void SparsePolyRingBase::HomImpl::myOutputSelfDetails(std::ostream& out) const
  {
    const SparsePolyRing P = myDomain;
    if (NumIndets(P) == 0) return;
    out << " sending "
        << "(" << indet(P, 0) << " |--> " << myIndetImages[0] << ")";
    const long n = len(myIndetImages);
    for (long i=1; i < n; ++i)
    {
      out << " & (" << indet(P, i) << " |--> " << myIndetImages[i] << ")";
    }
  }


  RingHom SparsePolyRingBase::myCoeffEmbeddingHomCtor() const
  {
    return RingHom(new CoeffEmbeddingHomImpl(SparsePolyRing(this)));
  }


  RingHom SparsePolyRingBase::myHomCtor(const ring& codomain, const RingHom& CoeffHom, const std::vector<RingElem>& IndetImages) const
  {
    // Args already sanity checked by PolyRingHom (see PolyRing.C)
// DON'T KNOW IF I REALLY WANT TO MAKE THIS CHECK...
//       // Check to see if we're building an identity homomorphism
//       if (ring(this) == codomain && IsIdentity(CoeffHom))
//       {
//         bool IndetsFixed = true;
//         for (long i=0; i < myNumIndetsValue; ++i)
//           IndetsFixed &= (myIndetVector[i] == IndetImages[i]);
//         if (IndetsFixed) return IdentityHom(ring(this));
//       }
      // General case
    return RingHom(new HomImpl(SparsePolyRing(this), codomain, CoeffHom, IndetImages));
  }


  RingHom SparsePolyRingBase::myCompose(const RingHom& phi, const RingHom& theta) const
  {
    vector<RingElem> IndetImages;
    for (long var=0; var < myNumIndets(); ++var)
      IndetImages.push_back(phi(theta(myIndets()[var])));

    return myHomCtor(codomain(phi), phi(theta(myCoeffEmbeddingHomCtor())), IndetImages);
  }


  bool SparsePolyRingBase::myImageLiesInSubfield(const RingHom& phi) const
  {
    CoCoA_ERROR(ERR::NYI, "SparsePolyRingBase::myImageLiesInSubfield");
    return false;
  }


  //-- CoeffEmbeddingHomImpl ----------------------------------------

  //---------------------------------------------------------------------------
  // Functions for the class SparsePolyRingBase::CoeffEmbeddingHomImpl


  SparsePolyRingBase::CoeffEmbeddingHomImpl::CoeffEmbeddingHomImpl(const SparsePolyRing& P):
    RingHomEmbeddingBase(CoeffRing(P), P)
  {}


  void SparsePolyRingBase::CoeffEmbeddingHomImpl::myApply(RingElemRawPtr rawimage, RingElemConstRawPtr rawarg) const
  {
    const SparsePolyRing P = myCodomain;
    RingElem ans(P);  // don't use image here for aliasing
    // ??? ANNA profile this:  (probably better to have myMonomial)
    if (!myDomain->myIsZero(rawarg))
      ans = monomial(P, RingElemAlias(myDomain, rawarg), one(PPM(P)));
    P->mySwap(rawimage, raw(ans));
  }


} // end of namespace CoCoA

// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/SparsePolyRing.C,v 1.150 2014/07/31 13:10:45 bigatti Exp $
// $Log: SparsePolyRing.C,v $
// Revision 1.150  2014/07/31 13:10:45  bigatti
// -- GetMatrix(PPO) --> OrdMat(PPO)
// -- added OrdMat and GradingMat to PPOrdering, PPMonoid, SparsePolyRing
//
// Revision 1.149  2014/07/30 14:09:59  abbott
// Summary: Changed name AmbientRing --> RingOf
// Author: JAA
//
// Revision 1.148  2014/07/28 15:52:10  abbott
// Summary: Redesign: ringhoms no longer cached in rings (caused ref count trouble); added general defn of myCoeffEmbeddingHomCtor
// Author: JAA
//
// Revision 1.147  2014/07/11 15:37:16  bigatti
// -- default implementation of myOutputSelf(OpenMathOutput& OMOut)
//
// Revision 1.146  2014/07/09 11:42:41  abbott
// Summary: Removed some cruft
// Author: JAA
//
// Revision 1.145  2014/07/08 15:23:34  abbott
// Summary: Updated comment
// Author: JAA
//
// Revision 1.144  2014/07/08 13:14:41  abbott
// Summary: Removed AsQuotientRing; added new defn of BaseRing
// Author: JAA
//
// Revision 1.143  2014/07/08 09:02:09  abbott
// Summary: Corrected silly typo
// Author: JAA
//
// Revision 1.142  2014/07/08 08:38:39  abbott
// Summary: Removed AsFractionField
// Author: JAA
//
// Revision 1.141  2014/07/07 17:11:57  abbott
// Summary: [MAJOR CHANGE] Removed AsSparsePolyRing; added SparsePolyRingPtr
// Author: JAA
//
// Revision 1.140  2014/07/04 12:55:30  abbott
// Summary: Revised following redefinition of symbols for a PPM
// Author: JAA
//
// Revision 1.139  2014/05/16 13:24:42  bigatti
// -- fixed problem with elim(.., ideal(0)) (?)
//
// Revision 1.138  2014/05/06 16:01:26  abbott
// Summary: Commented out unused "suspect" function MonomialMulTrick
// Author: JAA
//
// Revision 1.137  2014/05/06 13:20:41  abbott
// Summary: Changed names (my)MaxExponents into (my)Deg
// Author: JAA
//
// Revision 1.136  2014/04/30 16:13:56  abbott
// Summary: Replaced X.size() by len(X)
// Author: JAA
//
// Revision 1.135  2014/04/22 12:45:16  abbott
// Summary: Added commente out "general soln" for myDivMod
// Author: JAA
//
// Revision 1.134  2014/04/15 13:28:55  abbott
// Summary: Added new ClearDenom fn (with 2 args)
// Author: JAA
//
// Revision 1.133  2014/04/11 15:07:39  abbott
// Summary: Renamed TmpFactor to factor in include
// Author: JAA
//
// Revision 1.132  2014/03/27 14:57:23  bigatti
// -- added myMinimalize
//
// Revision 1.131  2014/03/21 13:07:51  bigatti
// -- added check in myIntersect
// -- cosmetics
//
// Revision 1.130  2014/02/26 15:35:57  abbott
// Summary: Changed some 1 char strings into chars
// Author: JAA
//
// Revision 1.129  2014/01/28 13:29:25  bigatti
// -- moved back check for IsInvertible in gcd (so that gcd(inv,0) --> 1)
//
// Revision 1.128  2014/01/28 13:09:06  bigatti
// -- added LF for ideal
//
// Revision 1.127  2014/01/28 09:47:56  abbott
// Improved impls of PushFront & PushBack (avoids double creation of PP).
//
// Revision 1.126  2013/10/28 13:21:58  bigatti
// -- fixed error message
//
// Revision 1.125  2013/06/28 17:03:51  abbott
// Modified semantics of IdealBase::myDivMod;
// it now returns a boolean.
// Several consequential changes.
//
// Revision 1.124  2013/06/27 16:56:09  abbott
// Changed arg checking fn to expect a "const char*" const instead of a "std::string"
// (as the ctor cost for a string is paid whenever the checking fn is called, even if
// there is no error to report).
//
// Revision 1.123  2013/06/12 08:49:17  bigatti
// -- added MinGens
//
// Revision 1.122  2013/05/30 13:13:58  bigatti
// -- added SetGBasisAsGens
//
// Revision 1.121  2013/05/30 07:24:10  bigatti
// -- improved myPrimeTest (waiting for IsMaximal3)
//
// Revision 1.120  2013/05/28 13:29:54  abbott
// Added printing for struct CoeffPP.
//
// Revision 1.119  2013/05/27 17:16:37  abbott
// Added missing const.
//
// Revision 1.118  2013/05/27 16:25:56  abbott
// Fixed "typo" (forgot to rename a couple of variables).
// Renamed args and local variables in myGCD.
//
// Revision 1.117  2013/05/27 15:00:11  abbott
// Added some arg checks to CoefficientsWRT and ContentWRT.
//
// Revision 1.116  2013/05/27 07:07:49  bigatti
// -- fixed all cases with 0's
//
// Revision 1.115  2013/05/24 17:55:19  bigatti
// -- fixed myColon for input with 0
//
// Revision 1.114  2013/03/25 17:04:19  abbott
// Major clean-up of interface to SmallFpImpl/SmallFpLogImpl/SmallFpDoubleImpl
// (underlying impl remains much the same).  Removed lots of cruft.
// Consequential changes to RingFp* classes; small change to SparsePolyRing.
//
// Revision 1.113  2013/02/21 17:35:00  bigatti
// -- simplified gcd (call to SyzOfGens instead of ComputeSyz)
//
// Revision 1.112  2013/02/21 12:51:42  abbott
// Added new fn UnivariateIndetIndex.
//
// Revision 1.111  2013/02/12 16:29:15  bigatti
// -- added ERR::ZeroGradingDim in IsHomog
//
// Revision 1.110  2013/01/18 18:16:02  abbott
// Modified CoefficientsWRT so that it now gives result guaranteeing that the entries are in decsreasing order of PP.
//
// Revision 1.109  2013/01/17 15:11:19  abbott
// Added new fn CoeffVecWRT.
// Added new virt mem fns myImageLiesInSubfield & IamPartial.
//
// Revision 1.108  2012/10/24 13:38:52  abbott
// Changed ConstRefRingElem into RingElemAlias in ctor calls and some local variable types.
//
// Revision 1.107  2012/10/17 09:40:15  abbott
// Replaced  RefRingElem  by  RingElem&
// (plus a few consequential changes)
//
// Revision 1.106  2012/10/05 10:21:39  bigatti
// -- added LF (leading form)
//
// Revision 1.105  2012/10/03 12:28:46  bigatti
// -- fixed error message (homog)
//
// Revision 1.104  2012/10/02 16:44:04  bigatti
// -- added homog for ideal
//
// Revision 1.103  2012/05/30 13:44:11  bigatti
// -- renamed IhaveMonomialGensB3Value --> IhaveMonomialGens3Value
//
// Revision 1.102  2012/05/29 07:45:23  abbott
// Implemented simplification change to bool3:
//  changed names of the constants,
//  changed names of the testing fns.
//
// Revision 1.101  2012/05/28 09:18:20  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.100  2012/05/24 14:49:23  bigatti
// -- changed symbol "index" into "subscripts"
//
// Revision 1.99  2012/05/22 10:02:37  abbott
// Removed IsGCDDomain; substituted by IsTrueGCDDomain.
// Added IsFractionFieldOfGCDDomain.
//
// Revision 1.98  2012/04/27 15:06:48  abbott
// Corrected myMaximalTest
//
// Revision 1.97  2012/04/15 20:10:13  abbott
// Added quick hacks to use old "faster" code for polynomial GCD in some cases.
//
// Revision 1.96  2012/04/13 16:27:51  abbott
// Added a working simple case in myDivMod; enough for Laura to proceed.
//
// Revision 1.95  2012/04/12 12:22:40  abbott
// Corrected handling of 0 in myNormalizeFracNoGcd;
// previously  0*x/(-1)  produced the result 0/(-1).
//
// Revision 1.94  2012/04/11 10:47:20  abbott
// Added special handling for the case 0/1 in myNormalizeFracNoGcd
// (without it the code wrongly divides by zero).
//
// Revision 1.93  2012/04/04 08:47:53  bigatti
// -- improved error messages in myGcd
//
// Revision 1.92  2012/04/03 16:12:46  abbott
// Changed slightly the fn signatures of CoefficientsWRT.
// Added some things to *.C file -- will complete later.
//
// Revision 1.91  2012/03/02 14:21:46  bigatti
// -- added ContentWRT, CoefficientsWRT(f, x)
//
// Revision 1.90  2012/02/24 13:10:16  abbott
// Removed cruft (commented out useless code).
//
// Revision 1.89  2012/02/14 15:16:14  bigatti
// -- commented out unused code
//
// Revision 1.88  2012/02/10 17:09:38  abbott
// Added new fns  indets, CoefficientsWRT, ContentWRT.
//
// Revision 1.87  2012/02/10 10:29:07  bigatti
// -- changed RingZ.H, RingQ.H --> RingZZ.H, RingQQ.H
//
// Revision 1.86  2012/02/08 15:14:13  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.85  2012/01/26 16:47:16  bigatti
// -- changed back_inserter into insert
//
// Revision 1.84  2011/12/05 16:32:11  bigatti
// -- fixed bug about saturation (by non-principal ideal)
//
// Revision 1.83  2011/11/09 14:29:37  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.82  2011/11/07 10:55:23  bigatti
// -- AreMonomials is now public
//
// Revision 1.81  2011/08/24 10:29:55  bigatti
// -- renamed QQ --> BigRat
// -- sorted #include
//
// Revision 1.80  2011/08/14 15:52:16  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.79  2011/08/12 15:22:49  abbott
// Improved a comment.
//
// Revision 1.78  2011/07/07 10:03:54  abbott
// Removed "inline" before ourGetPtr (as compilation fails with -O2).
//
// Revision 1.77  2011/07/05 15:02:17  bigatti
// -- added AlexanderDual
// -- added ad-hoc functions for colon, elim on monomial ideals
//
// Revision 1.76  2011/06/27 13:30:37  bigatti
// -- new file for monomial ideals:
// -- some functions moved there, but some have to be declared in SparsePolyRing.H
//
// Revision 1.75  2011/06/23 16:04:46  abbott
// Added IamExact mem fn for rings.
// Added myRecvTwinFloat mem fn for rings.
// Added first imple of RingHom from RingTwinFloat to other rings.
//
// Revision 1.74  2011/05/26 16:33:47  bigatti
// -- improved implementation of IamZero()
// -- check for IamZero before computing GBasis
//
// Revision 1.73  2011/05/20 16:12:59  bigatti
// -- added QuotientBasis
//
// Revision 1.72  2011/04/27 08:21:48  bigatti
// -- added gcd with coefficients in GCDDomain
//
// Revision 1.71  2011/04/19 13:59:49  bigatti
// -- added AreGensMonomial
//
// Revision 1.70  2011/04/12 09:52:54  bigatti
// -- added IsHomog(ideal), LT(ideal)
//
// Revision 1.69  2011/03/22 15:28:20  bigatti
// -- fixed IsIndetPosPower
//
// Revision 1.68  2011/03/16 15:40:22  bigatti
// -- added myIsIndetPosPower(f), IsIndetPosPower(f)
//
// Revision 1.67  2011/03/16 13:21:23  abbott
// Added comments for myIsHomogPartial & myCmpWDegPartial.
// Cleaned up impls of myIsHomog & myIsHomogPartial.
// Corrected typo in impl of myElim.
//
// Revision 1.66  2011/03/11 10:54:03  bigatti
// -- added mySaturate
//
// Revision 1.65  2011/03/10 16:39:33  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.64  2011/03/01 14:10:47  bigatti
// -- added ClearDenom/myClearDenom
//
// Revision 1.63  2011/02/28 14:16:08  bigatti
// -- added GBasis(ideal)  -- only for SparsePolyRing
// -- error for myGcd when CoeffRing is not a field
//
// Revision 1.62  2011/02/23 16:04:50  bigatti
// -- more fixes to get gcd==1 when coprime
//
// Revision 1.61  2011/02/23 15:02:22  bigatti
// -- gcd for polys now returns 1 when coprime (i.e. not any invertible constant)
//
// Revision 1.60  2011/02/04 11:05:06  bigatti
// -- changed printing for coefficients in QQ: added parentheses
//
// Revision 1.59  2011/01/28 17:58:07  bigatti
// -- added myElim
//
// Revision 1.58  2011/01/28 11:41:42  bigatti
// -- added IsPrintedWithMinus
// -- improved myOutput
// -- fixed bug in IsMinusOne
//
// Revision 1.57  2011/01/18 14:38:43  bigatti
// -- moved **_forC5 functions into CoCoA-5/CoCoALibSupplement:
//    myMonomials_forC5, mySupport_forC5, monomials_forC5, support_forC5,
//    LPP_forC5, LT_forC5, LM_forC5
//
// Revision 1.56  2010/11/30 11:28:27  bigatti
// -- moved IndetsCalled into unique implementation in PolyRing
//
// Revision 1.55  2010/11/25 12:31:22  bigatti
// -- added myIndetsCalled
//
// Revision 1.54  2010/11/05 15:59:31  bigatti
// -- added myMonomials_forC5, mySupport_forC5
//
// Revision 1.53  2010/11/02 15:33:23  bigatti
// -- fixed myIsPrintAtom
//
// Revision 1.52  2010/10/08 11:39:52  abbott
// Renamed DistrMPoly to DistrMPolyClean.
//
// Revision 1.51  2010/10/01 15:52:23  bigatti
// -- added mySymbolValue
//
// Revision 1.50  2010/06/10 08:00:02  bigatti
// -- fixed naming conventions
//
// Revision 1.49  2010/03/30 16:06:10  bigatti
// -- fixed call to characteristic
//
// Revision 1.48  2010/03/30 15:17:51  bigatti
// -- added nasty trick for simplifying multiplication of monomial ideals
//
// Revision 1.47  2010/03/18 13:55:56  abbott
// Added pseudo-ctors for monomials with QQ coeffs.
//
// Revision 1.46  2010/03/05 18:43:48  abbott
// Added pseudo-ctors allowing polynomial rings to be created specifying
// the ordering using a PPOrderingCtor object.
//
// Revision 1.45  2010/02/04 10:14:38  bigatti
// -- fixed "mul" for ideal
//
// Revision 1.44  2010/02/04 09:57:11  bigatti
// -- added "mul" for ideals.  Implemented only for SparsePolyRing
//
// Revision 1.43  2010/02/03 14:15:58  bigatti
// -- moved IsMonomial(vector<RingElem>) into anomymous namespace
//
// Revision 1.42  2010/01/21 13:17:08  bigatti
// -- commented out unnecessary default definition for bool3
//
// Revision 1.41  2010/01/20 16:49:52  bigatti
// -- fixed IhaveMonomialGens3Value in some functions modifying the ideal
//
// Revision 1.40  2009/11/26 17:21:33  bigatti
// -- added PushFront/PushBack(f, c, pp)
// -- added in .C inline functions:
// --   CheckCompatible, CheckElemSparsePolyRing, CheckCoeffExpv, CheckCoeffPP
//
// Revision 1.39  2009/10/29 18:43:23  abbott
// Added necessary include directive for iostream
// (previously was wrongly commented out).
//
// Revision 1.38  2009/10/02 13:27:26  bigatti
// -- unique implementation of myDiv in PolyRing.C
//
// Revision 1.37  2009/09/28 16:19:43  bigatti
// -- unique implementation for myDeriv
//
// Revision 1.36  2009/09/25 13:09:27  bigatti
// -- fixed: reinserted #include QuotientRing
//
// Revision 1.35  2009/09/25 13:02:09  bigatti
// -- myDiv with one implementation in SparsePolyRing
//
// Revision 1.34  2009/09/24 16:24:32  abbott
// Added include directive (after removing it from RingFp.H).
//
// Revision 1.33  2009/09/22 14:07:33  bigatti
// -- added CmpWDegPartial and IsHomogPartial
//
// Revision 1.32  2009/07/30 15:43:56  bigatti
// -- new implementation of IsMonomial(vector<RingElem>) using STL functions
//
// Revision 1.31  2009/07/24 12:26:42  abbott
// Added CommonDenom function for polynomials.
//
// Revision 1.30  2009/07/02 16:32:11  abbott
// Consequential changes stemming from new class QQ, and modified interface to the member
// function RingBase::myIsRational.  Also some new conversion functions.
//
// Revision 1.29  2009/05/22 10:27:08  bigatti
// -- fixed myOutput and myIsPrintAtom
//
// Revision 1.28  2009/05/21 14:49:33  abbott
// Cleaned up myOutput: now handles -1 in small finite fields correctly
// (if residues are symmetric).
//
// Revision 1.27  2008/12/17 12:11:52  abbott
// Changed type from long to MachineInt in operations which use a machine integer
// in place of a RingElem.  The change is "superficial" but affects many files.
//
// Revision 1.26  2008/11/18 15:20:09  bigatti
// -- added const to myGBasis return value
// -- added myIdealCtor to RingWeyl for proper inheritance
//
// Revision 1.25  2008/10/09 15:08:42  bigatti
// -- added MonomialIntersection
//
// Revision 1.24  2008/07/04 09:11:04  bigatti
// -- new PPVector class
//
// Revision 1.23  2008/05/29 16:02:32  bigatti
// -- improved myMul for the case [coeff * poly]
//
// Revision 1.22  2008/05/27 16:43:35  bigatti
// -- minor change
//
// Revision 1.21  2008/05/27 16:25:59  bigatti
// -- moved code for monomial ideals into TmpPPVector
// -- small improvement for myMul
//
// Revision 1.20  2008/04/15 15:53:07  bigatti
// -- new for SparsePolyRing: mySquare, first draft
//
// Revision 1.19  2008/04/10 16:36:03  bigatti
// -- prototype for squaring polynomials
//
// Revision 1.18  2008/04/03 13:07:16  bigatti
// -- improved "myPowerSmallExp" for monomial input
//
// Revision 1.17  2008/03/12 16:35:18  bigatti
// -- changed: IsHomogeneous --> IsHomog
// -- changed: ERR:ZeroPoly --> ERR::ZeroRingElem
//
// Revision 1.16  2007/12/07 15:27:01  bigatti
// -- default implementation of "IamOne" in ideal.C
//
// Revision 1.15  2007/12/05 11:06:24  bigatti
// -- changed "size_t StdDeg/myStdDeg(f)" into "long"  (and related functions)
// -- changed "log/myLog(f, i)" into "MaxExponent/myMaxExponent(f, i)"
// -- fixed bug in "IsOne(ideal)" in SparsePolyRing.C
//
// Revision 1.14  2007/12/04 14:27:06  bigatti
// -- changed "log(pp, i)" into "exponent(pp, i)"
//
// Revision 1.13  2007/10/30 17:14:07  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.12  2007/10/05 14:37:42  bigatti
// -- just fixed a comment
//
// Revision 1.11  2007/09/24 14:25:14  abbott
// Changed IsIndetPower to IsIndetPosPower.
// Added parens to shut up gcc-4.3.
//
// Revision 1.10  2007/06/25 11:38:15  bigatti
// -- reverted IsIndetPosPower --> IsIndetPower
//
// Revision 1.9  2007/06/21 21:29:47  abbott
// Changed name of RingFloat into RingTwinFloat.
//
// Revision 1.8  2007/05/31 16:34:37  bigatti
// -- Changed IsValid (now returns true of false and does not throw an error)
// -- using IsValid for sanity check in PushBack
//
// Revision 1.7  2007/05/31 15:43:56  bigatti
// -- added mySymbols and AreGoodIndetNames
//
// Revision 1.5  2007/05/22 22:45:14  abbott
// Changed fn name IsUnit to IsInvertible.
//
// Revision 1.4  2007/05/21 14:46:34  bigatti
// -- PushFront now accepts zero coefficient
//
// Revision 1.3  2007/03/12 16:00:29  bigatti
// -- moved myLog(F, index) into unique implementation in SparsePolyRing
//
// Revision 1.2  2007/03/09 18:56:56  bigatti
// -- added Tmp prefix to Groebner related files
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.32  2007/03/08 18:49:20  bigatti
// -- removed TmpGcd.C: moved code into SparsePolyRing.C
//
// Revision 1.31  2007/03/08 18:22:29  cocoa
// Just whitespace cleaning.
//
// Revision 1.30  2007/03/08 11:07:12  cocoa
// Made pseudo ctors for polynomial rings more uniform.  This allowed me to
// remove an include of CoCoA/symbol.H  from the RingDistrM*.H files, but then
// I had to put the include in several .C files.
//
// Revision 1.29  2007/03/07 17:04:31  cocoa
// -- several changes by M.Caboara: more operations on ideals,
//    exception cleaner, coding conventions, WSugar, dynamic
//
// Revision 1.28  2007/03/07 14:08:52  bigatti
// -- minor: commented argument names for -Wextra
//
// Revision 1.27  2007/03/05 21:06:07  cocoa
// New names for homomorphism pseudo-ctors: removed the "New" prefix.
//
// Revision 1.26  2007/03/02 10:47:53  cocoa
// First stage of RingZ modifications -- tests do not compile currently, Anna will fix this.
//
// Revision 1.25  2007/02/28 13:51:59  bigatti
// -- added function IsMonomial
//
// Revision 1.24  2007/02/26 15:00:01  bigatti
// -- just a comment fix
//
// Revision 1.23  2007/02/10 18:44:03  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.22  2007/01/20 14:07:25  bigatti
// -- moved code for homomorphism into common implementation in SparsePolyRing
//
// Revision 1.21  2006/12/07 17:36:19  cocoa
// -- migrated  myRemoveBigContent myContent myPowerSmallExp  into
//    single implementation in SparsePolyRing
// -- removed  content  from DistrMPoly(..)
//
// Revision 1.20  2006/11/24 17:06:10  cocoa
// -- reorganized includes of header files
//
// Revision 1.19  2006/11/23 17:48:43  cocoa
// -- minor change
//
// Revision 1.18  2006/11/22 17:51:31  cocoa
// -- moved printing functions into unified implementation in SparsePolyRing
//
// Revision 1.17  2006/11/22 15:12:48  cocoa
// -- minor cleaning (indicated by Intel compiler)
//
// Revision 1.16  2006/11/21 18:09:23  cocoa
// -- added myIsMonomial
// -- implemented myIsOne, myIsMinusOne, myIsConstant, myIsIndet in SparsePolyRing
// -- removed the 4 functions from DistrMPoly(..) and RingDistrMPoly(..)
// -- changed all names of RawPtr arguments into "raw(..)"
//
// Revision 1.15  2006/11/17 18:14:01  cocoa
// -- fixed: now Hilbert computes TidyGens instead of cheating
// -- myGBasis is much more efficient on monomial input
//
// Revision 1.14  2006/11/09 17:46:58  cocoa
// -- version 0.9712:
// --   IdealImpl moved to SparsePolyRing from concrete rings
// -- PolyList in GTypes is now vector<RingElem> (was list)
// -- "my" coding convention applied to DistrMPoly
//
// Revision 1.13  2006/10/16 23:18:59  cocoa
// Corrected use of std::swap and various special swap functions.
// Improved myApply memfn for homs of RingDistrMPolyInlPP.
//
// Revision 1.12  2006/10/06 14:04:14  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.11  2006/09/27 14:30:45  cocoa
// -- changed f into rawf for ConstRawPtr and RawPtr variables
//
// Revision 1.10  2006/08/17 09:45:07  cocoa
// -- added: homogenization
//
// Revision 1.9  2006/08/07 21:23:25  cocoa
// Removed almost all publicly visible references to SmallExponent_t;
// changed to long in all PPMonoid functions and SparsePolyRing functions.
// DivMask remains to sorted out.
//
// Revision 1.8  2006/07/20 17:06:08  cocoa
// -- moved myStdDeg into SparsePolyRing
//
// Revision 1.7  2006/07/20 16:51:38  cocoa
// -- added common implementation of myStdDeg
//
// Revision 1.6  2006/07/20 14:24:12  cocoa
// -- added special cases for myMul (0) and myGcd (unit)
//
// Revision 1.5  2006/07/17 19:32:55  cocoa
// Added arg check to PushFront.
//
// Revision 1.4  2006/07/17 11:05:53  cocoa
// -- added: myIsValid, myIsHomogeneous, IsHomogeneous
//
// Revision 1.3  2006/06/20 17:25:27  cocoa
// -- added function geobucket::myAddMul(monom, g, gLen);   [without SkipLMFlag]
//
// Revision 1.2  2006/06/08 16:45:27  cocoa
// -- RingDistrMPoly*.H  have been "moved" into RingDistrMPoly*.C
// -- some coding conventions fixed in DistrMPoly*
// -- functions wdeg and CmpWDeg have a common implementation in SparsePolyRing
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.11  2006/05/29 16:22:37  cocoa
// Third time lucky???
// Added myIsInteger member function to all rings (NYI for RingFloat).
//
// Revision 1.10  2006/05/12 17:00:34  cocoa
// -- added some function whose implementation in RingDistr*** were identical
// -- added multivariate gcd
//
// Revision 1.9  2006/05/12 13:16:30  cocoa
// Added functions for preprocessing approximate points.
//
// Revision 1.8  2006/04/27 14:08:24  cocoa
// -- sorry: missed parenthesis
//
// Revision 1.7  2006/04/27 14:05:29  cocoa
// -- changed CmpWDeg: I think it should throw an error if polys are 0
//
// Revision 1.6  2006/04/26 16:44:53  cocoa
// -- myMul has now a single implementation in SparsePolyRing
// -- myMul and mul in RingDistrMPoly* and DistrMPoly* have been disabled
//
// Revision 1.5  2006/04/21 14:58:04  cocoa
// Removed myWDeg member function: it is no longer needed.
//
// Revision 1.4  2006/04/11 14:16:29  cocoa
// -- reorganization of fns between reduce,SparsePolyRing,GPoly
// -- added optional "len" argument to myAssignReset in ReductionCog
//
// Revision 1.3  2006/02/13 13:19:01  cocoa
// -- fixed: "const PPMonoidElem&" --> "ConstRefPPMonoidElem"
//
// Revision 1.2  2006/01/19 16:34:42  cocoa
// -- added NF, myReduceMod functions (not yet tested)
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.3  2005/08/08 16:36:32  cocoa
// Just checking in before going on holiday.
// Don't really recall what changes have been made.
// Added IsIndet function for RingElem, PPMonoidElem,
// and a member function of OrdvArith.
// Improved the way failed assertions are handled.
//
// Revision 1.2  2005/07/19 15:30:20  cocoa
// A first attempt at iterators over sparse polynomials.
// Main additions are to SparsePolyRing, DistrMPoly*.
// Some consequential changes to PPMonoid*.
//
// Revision 1.1  2005/07/01 16:08:15  cocoa
// Friday check-in.  Major change to structure under PolyRing:
// now SparsePolyRing and DUPolyRing are separated (in preparation
// for implementing iterators).
//
// A number of other relatively minor changes had to be chased through
// (e.g. IndetPower).
//
