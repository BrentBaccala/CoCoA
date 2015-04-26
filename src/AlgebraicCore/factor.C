//   Copyright (c)  2006,2009  John Abbott

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


// These are from directory TmpFactorDir/
#include "mpz_alias.h"
#include "primes.h"
#include "FF.h"
#include "DUPFF.h"
#include "DUPFFfactor.h"
//#include "DUPFFlist.h"
#include "WGD.h"

// These are from directory TmpFactorDir/multivariate/
#include "DMPZgcd.h"
#include "DMPZfactor.h"
#include "DMPZfactor_modp.h"

/***************************************************************************/
/* CoCoALib includes */

#include "CoCoA/factor.H"

#include "CoCoA/BigInt.H"
#include "CoCoA/IntOperations.H"
#include "CoCoA/CanonicalHom.H"
#include "CoCoA/DUPFp.H"
#include "CoCoA/FractionField.H"
#include "CoCoA/NumTheory.H"
#include "CoCoA/PolyRing.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/convert.H"
#include "CoCoA/utils.H"

//#include <iostream>
#include <cstdlib>
using std::malloc;// someone actually uses malloc!!!!
#include <vector>
using std::vector;

/***************************************************************************/

namespace CoCoA
{


  DUPFF ConvertToDUPFF(ConstRefRingElem f)
  {
    CoCoA_ASSERT(IsSparsePolyRing(owner(f)));
    const SparsePolyRing Rx = owner(f);
    const ring R = CoeffRing(Rx);
    CoCoA_ASSERT(IsField(R) && IsQuotientRing(R) && IsZZ(BaseRing(R)));
    const BigInt P = characteristic(R);
    const long p = ConvertTo<long>(P);
    CoCoA_ASSERT(p < 32768); /// *** BUG BUG BUG ***

    FF Fp = FFctor(p);
    FFselect(Fp);
    const long d = StdDeg(f);
    DUPFF ans = DUPFFnew(d);
    { FFelem* CoeffVec = ans->coeffs; for (long i=0; i<=d; ++i) CoeffVec[i] = 0; }
    for (SparsePolyIter it = BeginIter(f); !IsEnded(it); ++it)
    {
      BigInt c;
      if (!IsInteger(c, coeff(it))) { FFdtor(Fp); CoCoA_ERROR("coeff not integer","ConvertToDUPFF"); }
      long coeff = ConvertTo<long>(c);
      if (coeff < 0) coeff += p;
      const long exp = StdDeg(PP(it));
      ans->coeffs[exp] = coeff;
    }
    ans->deg = d;
///??? FFdtor(Fp); // AMB 2013-01-29: leaks, but o/w lots of "Invalid read of size 4"
    return ans;
  }


  RingElem ConvertDUPFFToRingElem(DUPFF F, ConstRefRingElem x)
  {
    const SparsePolyRing P = owner(x);
    PPMonoidElem t = LPP(x);
    RingElem ans(P);
    const int D = DUPFFdeg(F);
    for (int i=D; i >= 0; --i)
    {
      if (F->coeffs[i] == 0) continue;
      ans += monomial(P, F->coeffs[i], power(t,i));
    }
    return ans;
  }

/////?????  DUPZ ConvertToDUPZ(ConstRefRingElem f) { ... }

// Exception safe???  You're kidding, right?
  DMPZ ConvertToDMPZ(ConstRefRingElem f)
  {
    CoCoA_ASSERT(IsSparsePolyRing(owner(f)));
    const SparsePolyRing P = owner(f);
    const ring R = CoeffRing(P);
    CoCoA_ASSERT(characteristic(R) == 0);
    CoCoA_ASSERT(IsZZ(R) || IsFractionFieldOfGCDDomain(R));

    const long nvars = NumIndets(P);
    DMPZ g = NULL;
    std::vector<long> expv(nvars);
    for (SparsePolyIter it = BeginIter(f); !IsEnded(it); ++it)
    {
      BigInt c;
      if (!IsInteger(c, coeff(it))) CoCoA_ERROR("coeff not in Z","factorize");
      exponents(expv, PP(it));
      int *exps = (int*)malloc(nvars*sizeof(int)); // block absorbed by g
      std::copy(expv.begin(), expv.end(), exps);
      g = DMPZprepend(mpzref(c), exps, g);
    }

    return DMPZreverse(g); /* restore original order of the terms */
  }


  void ConvertDMPZToRingElem(RingElem& f, DMPZ poly)
  {
    CoCoA_ASSERT(IsSparsePolyRing(owner(f)));
    const SparsePolyRing P = owner(f);
    const PPMonoid M = PPM(P);
    PPMonoidElem t(M);
    f = 0;
    const long nvars = NumIndets(P);
    std::vector<long> exps(nvars);
    for (DMPZ iter = poly; iter; iter = iter->next)
    {
      if (mpz_sgn(iter->coeff) == 0) continue;
      const BigInt c(iter->coeff);
      std::copy(iter->exps, iter->exps+nvars, exps.begin());
      f += monomial(P, c, PPMonoidElem(M, exps));
    }
  }

  factorization<RingElem> ConvertToFactorList(DMPZfactors facpows, ring Rx)
  {
    RingElem content = RingElem(Rx, BigInt(facpows->content));
    RingElem factor(Rx);
    vector<RingElem> facs;
    vector<long> exps;

    DMPZlist iter;
    int i;
    for (i=1, iter=facpows->list; iter; ++i, iter = iter->next)
    {
      ConvertDMPZToRingElem(factor, iter->poly);
      facs.push_back(factor);
      exps.push_back(iter->deg);
    }
    DMPZfactors_dtor(facpows);  // NASTY - we destroy one of our args!!
    return factorization<RingElem>(facs, exps, content);
  }


  RingElem GCD_DUPFF(ConstRefRingElem f, ConstRefRingElem g)
  {
    CoCoA_ASSERT(owner(f) == owner(g));
    CoCoA_ASSERT(IsSparsePolyRing(owner(f)));
    const SparsePolyRing P = owner(f);
    const vector<long> index = indets(f);
    const RingElem x = indet(P, index[0]);
    DUPFF F = ConvertToDUPFF(f);
    DUPFF G = ConvertToDUPFF(g);
    DUPFFgcd2(F,G);
    RingElem ans = ConvertDUPFFToRingElem(F, x);
    DUPFFfree(G);
    DUPFFfree(F);
    return ans;
  }


  RingElem GCD_DMPZ(ConstRefRingElem f, ConstRefRingElem g)
  {
    CoCoA_ASSERT(owner(f) == owner(g));
    CoCoA_ASSERT(IsSparsePolyRing(owner(f)));
    const SparsePolyRing P = owner(f);
    CoCoA_ASSERT(IsZZ(CoeffRing(P)) || IsQQ(CoeffRing(P)));

    if (IsQQ(CoeffRing(P)))
    {
      // case P = QQ[x,y,z,...]
      RingHom coeff = CoeffEmbeddingHom(P);
      RingElem ContF = content(f);
      RingElem ContG = content(g);
      NVARS = NumIndets(P); // BUG BUG BUG:  NVARS is a global!!!
      DMPZ F = ConvertToDMPZ(f/coeff(ContF));
      DMPZ G = ConvertToDMPZ(g/coeff(ContG));
      DMPZ H = DMPZgcd(F,G);
      DMPZdtor(G);
      DMPZdtor(F);

      RingElem ans(P);
      ConvertDMPZToRingElem(ans, H);
      DMPZdtor(H);
      return ans;
    }
    if (IsZZ(CoeffRing(P)))
    {
      // case P = ZZ[x,y,z,...]
      NVARS = NumIndets(P); // BUG BUG BUG:  NVARS is a global!!!
      DMPZ F = ConvertToDMPZ(f);
      DMPZ G = ConvertToDMPZ(g);
      DMPZ H = DMPZgcd(F,G);
      DMPZdtor(G);
      DMPZdtor(F);

      RingElem ans(P);
      ConvertDMPZToRingElem(ans, H);
      DMPZdtor(H);
      return ans;
    }

    CoCoA_ERROR(ERR::SERIOUS, "GCD_DMPZ");
    return f; // just to keep compiler quiet
  }


  factorization<RingElem> factor(ConstRefRingElem f)
  {
    ring R = owner(f);
    if (IsField(R) || !IsTrueGCDDomain(R))
      CoCoA_ERROR("factorization not supported", "factor(RingElem)");
    if (IsZero(f))
      CoCoA_ERROR(ERR::NotNonZero, "factor(x)");
    if (IsZZ(R))
    {
      // Factorization is actually over ZZ, so use factor(BigInt) to do the work!
      const factorization<BigInt> IntFacs = factor(ConvertTo<BigInt>(f));
      const vector<BigInt>& prime = IntFacs.myFactors();
      const long nprimes = len(prime);
      vector<RingElem> facs;
      for (long i=0; i < nprimes; ++i)
        facs.push_back(RingElem(R, prime[i]));
      return factorization<RingElem>(facs, IntFacs.myMultiplicities(), RingElem(R, IntFacs.myRemainingFactor()));
    }
    if (!IsPolyRing(R))
      CoCoA_ERROR(ERR::NotPolyRing, "factor(RingElem)");
/////    if (!IsSparsePolyRing(R))
/////      CoCoA_ERROR(ERR::NYI, "factor(RingElem)  NYI over DUP");
    const SparsePolyRing P = owner(f);
    if (IsQQ(CoeffRing(P)))
    {
      // case P = QQ[x,y,z,...]
      const FractionField Q = CoeffRing(P);
      const RingHom coeff = CoeffEmbeddingHom(P)(EmbeddingHom(Q));
      const RingElem D = CommonDenom(f);
      NVARS = NumIndets(P); // NVARS is a global!!!
      const DMPZ g = ConvertToDMPZ(f*coeff(D));
      factorization<RingElem> ans = ConvertToFactorList(DMPZfactor(g), P); // BUG BUG BUG leaks the factorlist
      DMPZdtor(g);
      if (!IsOne(D)) ans.myNewRemainingFactor(ans.myRemainingFactor()/CanonicalHom(RingZZ(), P)(D));
      return ans;
    }
    if (IsZZ(CoeffRing(P)))
    {
      // We do not attempt to factorize the content in ZZ[x,y,z] as it could easily be very costly.
      PolyRing QQx = NewPolyRing(RingQQ(), NumIndets(P));
      RingHom phi = PolyRingHom(P, QQx, ZZEmbeddingHom(RingQQ()), indets(QQx)); // maps ZZ[x,y,z] --> QQ[x,y,z]
      RingElem FinQQx = phi(f);
      factorization<RingElem> FactorizationInQQx = factor(FinQQx);
      const vector<RingElem>& FacsInQQx = FactorizationInQQx.myFactors();
      vector<RingElem> FacsInP;
      const RingElem ContentInP = CoeffEmbeddingHom(P)(num(FactorizationInQQx.myRemainingFactor()));
      RingHom psi = PolyRingHom(QQx, P, QQEmbeddingHom(RingZZ()), indets(P)); // maps QQ[x,y,z] --> ZZ[x,y,z]
      const long NumFacs = len(FacsInQQx);
      for (long i=0; i < NumFacs; ++i)
      {
        FacsInP.push_back(psi(FacsInQQx[i]));
      }
      return factorization<RingElem>(FacsInP, FactorizationInQQx.myMultiplicities(), ContentInP);
    }
    if (IsFiniteField(CoeffRing(P)))
    {
      // Case FF(q)[x,y,z,...]
      ring FFq = CoeffRing(P);
      if (LogCardinality(FFq) > 1)
        CoCoA_ERROR(ERR::NYI, "Factorization over alg extn of finite field");
      long var;
      if ((var = UnivariateIndetIndex(f)) >= 0)
      {
        const RingElem x = indet(P, var);
        SmallFpImpl ModP(ConvertTo<long>(characteristic(P)));
        DUPFp F = ConvertToDUPFp(ModP, f);
//std::clog<<"F="<<F<<std::endl;
        const factorization<DUPFp> facs = factor(F);
        const vector<DUPFp>& DUPfac = facs.myFactors(); // convenient alias
//std::clog<<"facs="<<facs<<std::endl; // doesn't compile for some reason...?
        const long n = len(facs.myFactors());
//std::clog<<"numfacs="<<n<<std::endl;
        vector<RingElem> irreds; irreds.reserve(n);
        vector<long> multiplicity; multiplicity.reserve(n);
        for (long i=0; i < n; ++i)
        {
//std::clog<<"LOOP fac="<<ConvertFromDUPFp(x, facs.myFactor(i))<<std::endl;
          irreds.push_back(ConvertFromDUPFp(x, DUPfac[i]));
        }
        return factorization<RingElem>(irreds, facs.myMultiplicities(), CoeffEmbeddingHom(P)(LC(f)));

//         DUPFF F = ConvertToDUPFF(f);
//         DUPFFlist facs = DUPFFfactor(F);
//         DUPFFfree(F);
//         const long n = DUPFFlist_length(facs);
//         vector<RingElem> irreds; irreds.reserve(n);
//         vector<long> multiplicity; multiplicity.reserve(n);
//         for (DUPFFlist curr=facs; curr != 0/*nullptr*/; curr = curr->next)
//         {
//           irreds.push_back(monic(ConvertDUPFFToRingElem(curr->poly, x)));
//           multiplicity.push_back(curr->deg);
//         }
//         DUPFFlist_dtor(facs);
//         return factorization<RingElem>(irreds, multiplicity, CoeffEmbeddingHom(P)(LC(f)));
      }
      CoCoA_ERROR(ERR::NYI, "Multivariate factorization over finite field");
    }
    CoCoA_ERROR(ERR::NYI, "Factorization not in QQ[x,y,...]  or in ZZ/(p)[x,y,...]");
    return factorization<RingElem>(vector<RingElem>(),vector<long>(),RingElem(P)); // NEVER EXECUTED, just to keep compiler quiet
  }


  // temporary: by AMB 2008-01-25
  bool IsIrredPoly(ConstRefRingElem f)
  { 
    if (!IsPolyRing(owner(f)))
      CoCoA_ERROR(ERR::NotPolyRing, "IsIrredPoly");
 // ??? warning: factorize syntax might change: now only for multivariate
    if (IsConstant(f)) return IsIrred(LC(f));
    // ANNA: should do it for Z as well (taking care of constant factors)
    const PolyRing P = owner(f);
    if (!IsField(CoeffRing(P)))
      CoCoA_ERROR(ERR::NYI, "IsIrredPoly, !IsField(CoeffRing)");
    if ( IsSparsePolyRing(P) )
    {
      const factorization<RingElem> facs = factor(f);
      const long NumFacs = len(facs.myFactors());
      const long mult = facs.myMultiplicities()[0];
      return (IsInvertible(facs.myRemainingFactor()) && NumFacs == 1 && mult == 1);
    }
    else
    {
      PolyRing SPR = NewPolyRing(CoeffRing(P), NumIndets(P));
      RingHom phi(PolyAlgebraHom(P, SPR, indets(SPR)));
      const factorization<RingElem> facs = factor(phi(f));
      const long NumFacs = len(facs.myFactors());
      const long mult = facs.myMultiplicities()[0];
      return (IsInvertible(facs.myRemainingFactor()) && NumFacs == 1 && mult == 1);
    }

  }


  factorization<RingElem> SqFreeFactor_generic(ConstRefRingElem x)
  {
    return factor(x);
  }


  // Simple rather than super-efficient.
  RingElem IteratedPthRoot(RingElem f)
  {
    CoCoA_ASSERT(IsSparsePolyRing(owner(f)) && IsFiniteField(CoeffRing(owner(f))));
    while (IsPthPower(f))
      f = PthRoot(f);
    return f;
  }

  factorization<RingElem> SqFreeFactorPosDerChar0(ConstRefRingElem f, long x)
  {
    const RingElem derivf = deriv(f,x);
    RingElem RemainingFacs = monic(gcd(f,derivf));
    RingElem rad = f/RemainingFacs;
    RingElem DerivCofactor = derivf/RemainingFacs;

    RingElem U = DerivCofactor - deriv(rad, x);
    long i = 1;

    vector<RingElem> facs;
    vector<long> exps;
    while (!IsZero(U))
    {
      const RingElem G = monic(gcd(rad, U));
      if (!IsConstant(G))
      {
        facs.push_back(G);
        exps.push_back(i);
      }
      ++i;
      rad /= G;
      RemainingFacs /= rad;
      DerivCofactor = U/G;
      U = DerivCofactor - deriv(rad,x);
    }
    if (!IsConstant(rad))
    { facs.push_back(monic(rad)); exps.push_back(i); }

    return factorization<RingElem>(facs, exps, RemainingFacs);
  }


  factorization<RingElem> SqFreeFactorChar0(RingElem f)
  {
    const PolyRing P = owner(f);
    const long nvars = NumIndets(P);
    const RingElem LCF = CoeffEmbeddingHom(P)(LC(f));
    f /= LCF; // good idea or bad???
    vector<RingElem> facs;
    vector<long> exps;
    for (long x=0; x < nvars; ++x)
    {
///???      if (deg(f,x) == 0) continue;
      if (deg(f,x) == 0) continue;
      const factorization<RingElem> PartialFacs = SqFreeFactorPosDerChar0(f, x);
      f = PartialFacs.myRemainingFactor();
      facs.insert(facs.end(), PartialFacs.myFactors().begin(), PartialFacs.myFactors().end());
      exps.insert(exps.end(), PartialFacs.myMultiplicities().begin(), PartialFacs.myMultiplicities().end());
    }
    return factorization<RingElem>(facs, exps, LCF);
  }

  factorization<RingElem> SqFreeFactorPosDerCharP(ConstRefRingElem f, long x)
  {
    const BigInt P = characteristic(owner(f));
    const RingElem derivf = deriv(f,x);
    RingElem RemainingFacs = gcd(f,derivf);
    RingElem rad = f/RemainingFacs;
    RingElem DerivCofactor = derivf/RemainingFacs;

    RingElem U = DerivCofactor - deriv(rad, x);
    long i = 1;

    vector<RingElem> facs;
    vector<long> exps;
    while (i < P && !IsZero(U))
    {
      const RingElem G = monic(gcd(rad, U));
      if (!IsConstant(G))
      {
        facs.push_back(G);
        exps.push_back(i);
      }
      ++i;
      rad /= G;
      RemainingFacs /= rad;
      DerivCofactor = U/G;
      U = DerivCofactor - deriv(rad,x);
    }
    if (!IsConstant(rad))
    { facs.push_back(monic(rad)); exps.push_back(i); }

    return factorization<RingElem>(facs, exps, RemainingFacs);
  }


  factorization<RingElem> SqFreeFactorCharP(RingElem f)
  {
    const BigInt P = characteristic(owner(f));
    const long NumVars = NumIndets(owner(f));
    const RingElem LCF = CoeffEmbeddingHom(owner(f))(LC(f));
    f /= LCF; // good idea or bad???

    vector<RingElem> facs;
    vector<long> mults;
    for (long x=0; x < NumVars; ++x)
    {
////?????      if (deg(f,x) == 0) continue;
      if (deg(f,x) == 0) continue;
      const factorization<RingElem> tmp = SqFreeFactorPosDerCharP(f, x);
      facs.insert(facs.end(), tmp.myFactors().begin(), tmp.myFactors().end());
      mults.insert(mults.end(), tmp.myMultiplicities().begin(), tmp.myMultiplicities().end());
      f = tmp.myRemainingFactor();
    }

    if (deg(f) == 0)
    {
      return factorization<RingElem>(facs, mults, LCF);
    }
    // else deg(f) > 0, so f is a p-th power.
    long q = deg(f);
    f = IteratedPthRoot(f);
    q /= deg(f); // a power of p
    const factorization<RingElem> RootFacs = SqFreeFactorCharP(f);
    vector<RingElem> PthRootFacs = RootFacs.myFactors(); // QUICK HACK now that myFactors is read-only!
    // Now create a sort of GCDfreeBasis; we know that the entries of facs
    // are coprime, and also RootsFacs.myFactors are coprime.
    const long n = len(facs);
    for (long i=0; i < n; ++i)
    {
      for (long j=0; j < len(PthRootFacs); ++j)
      {
        if (IsConstant(PthRootFacs[j])) continue;
        const RingElem Q = monic(gcd(facs[i], PthRootFacs[j]));
        if (IsConstant(Q)) continue;
        facs[i] /= Q;
        PthRootFacs[j] /= Q;
        facs.push_back(Q);
        mults.push_back(mults[i]+q*RootFacs.myMultiplicities()[j]);
        if (IsConstant(facs[i])) break;
      }
    }

    // Remove all constant factors
    long last = len(facs)-1;
    while (last >= 0 && IsConstant(facs[last]))
      --last;
    for (long i=0; i < last; ++i)
    {
      if (!IsConstant(facs[i])) continue;
      facs[i] = facs[last]; //swap(facs[i], facs[last]);
      mults[i] = mults[last]; //swap(mults[i], mults[last]);
      --last;
    }
    facs.resize(last+1);
    mults.resize(last+1);


    // Append any (non-constant) new factors
    for (long j=0; j < len(PthRootFacs); ++j)
    {
      if (IsConstant(PthRootFacs[j])) continue;
      facs.push_back(PthRootFacs[j]);
      mults.push_back(q*RootFacs.myMultiplicities()[j]);
    }
    return factorization<RingElem>(facs, mults, LCF);
  }


  factorization<RingElem> SqFreeFactor(ConstRefRingElem f)
  {
    const char* const FnName = "SqFreeFactor";
    if (!IsPolyRing(owner(f)))
      CoCoA_ERROR(ERR::NYI, FnName);
    const PolyRing P = owner(f);
    if (!IsTrueGCDDomain(P))
      CoCoA_ERROR(ERR::NotTrueGCDDomain, FnName);
    if (IsZero(f))
      CoCoA_ERROR(ERR::NotNonZero, FnName);

    const ring R = CoeffRing(P);
    if (IsQQ(R)) return SqFreeFactorChar0(f);
    if (IsFiniteField(R)) return SqFreeFactorCharP(f);
    if (IsPrime(characteristic(R))) return SqFreeFactorCharP(f); // BUG BUG  DODGY!!!
    return SqFreeFactor_generic(f);
  }


  namespace // anonymous
  {
    void ContentFreeFactorLoop(vector<RingElem>& ans, RingElem& f, vector<bool> skip)
    {
      CoCoA_ASSERT(!IsZero(f));
      CoCoA_ASSERT(IsSparsePolyRing(owner(f)));
      const SparsePolyRing P = owner(f);
      CoCoA_ASSERT(IsTrueGCDDomain(CoeffRing(P)) || IsField(CoeffRing(P))); // coeffs are in field or GCD domain
      const vector<long> var = indets(f);
      if (var.empty()) return;
      for (int i=0; i < len(var); ++i)
      {
        if (skip[var[i]]) continue;
        RingElem c = ContentWRT(f, indet(P,var[i]));
        if (IsOne(c)) continue;
        f /= c;
        skip[i] = true;
        ContentFreeFactorLoop(ans, c, skip);
        CoCoA_ASSERT(IsOne(c));
      }
      if (StdDeg(f) == 0) return;
      // f is non-trivial; normalize it before adjoining it to the result.
      if (IsField(CoeffRing(P)))
      {
        ans.push_back(monic(f));
        return;
      }
      // otherwise IsTrueGCDDomain(CoeffRing(P))
      {
        RingHom embed = CoeffEmbeddingHom(P);
        ans.push_back(f/embed(content(f)));
        return;
      }
    }
  } // end of anonymous namespace

  factorization<RingElem> ContentFreeFactor(RingElem f)
  {
    const char* const FnName = "ContentFreeFactor";
    if (!IsPolyRing(owner(f)))
      CoCoA_ERROR(ERR::NotPolyRing, FnName);
    const SparsePolyRing P = owner(f);
    if (!IsField(CoeffRing(P)) && !IsTrueGCDDomain(CoeffRing(P)))
      CoCoA_ERROR(ERR::NotTrueGCDDomain, FnName); //???? what error to give here???
    if (IsZero(f))
      CoCoA_ERROR(ERR::NotNonZero, FnName);
    // Idea: could remove numerical content here... (is it worth it?)
    RingElem RemainingFactor = LC(f); // BEFORE calling ContentFreeFactorLoop
    vector<RingElem> facs;
    ContentFreeFactorLoop(facs, f, vector<bool>(NumIndets(P)));
    const vector<long> exps(len(facs), 1);
    for (int i=0; i < len(facs); ++i)
      RemainingFactor /= LC(facs[i]);
    return factorization<RingElem>(facs, exps, CoeffEmbeddingHom(P)(RemainingFactor));
  }


} // end of namespace CoCoA

// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/factor.C,v 1.7 2014/07/08 15:24:23 abbott Exp $
// $Log: factor.C,v $
// Revision 1.7  2014/07/08 15:24:23  abbott
// Summary: Removed AsQuotientRing from an assertion
// Author: JAA
//
// Revision 1.6  2014/07/08 08:40:18  abbott
// Summary: Removed AsFractionField
// Author: JAA
//
// Revision 1.5  2014/07/07 14:12:35  abbott
// Summary: Corrected silly typo
// Author: JAA
//
// Revision 1.4  2014/07/07 13:24:48  abbott
// Summary: Removed AsPolyRing; Removed AsSparsePolyRing
// Author: JAA
//
// Revision 1.3  2014/05/06 13:20:41  abbott
// Summary: Changed names (my)MaxExponents into (my)Deg
// Author: JAA
//
// Revision 1.2  2014/04/30 16:26:49  abbott
// Summary: Replaced size_t by long; Replaced X.size() by len(X)
// Author: JAA
//
// Revision 1.1  2014/04/11 15:06:38  abbott
// Summary: Renamed from TmpFactor.H/C
// Author: JAA
//
// Revision 1.31  2014/04/08 16:35:30  abbott
// Summary: Changed SqfreeFactor to SqFreeFactor
// Author: JAA
//
// Revision 1.30  2014/03/24 12:09:21  abbott
// Summary: Major revision to public interface of factorization template class
// Author: JAA
//
// Revision 1.29  2014/01/16 16:15:13  abbott
// Renamed  SqfrDecomp_generic  into  SqfrFactor_generic  & cleaned the code slightly.
//
// Revision 1.28  2013/10/22 14:02:09  abbott
// Added code for SqFreeFactor.
//
// Revision 1.27  2013/09/12 12:59:52  abbott
// Minor cleaning of factorize fn when applied to a BigInt.
//
// Revision 1.26  2013/04/16 15:51:00  bigatti
// -- added "leak" warning comment
//
// Revision 1.25  2013/02/21 12:55:23  abbott
// Improved factor: now works in Fp[x,y,z] provided given poly is univariate.
//
// Revision 1.24  2013/01/28 13:05:22  abbott
// Fixed two memory leaks in GCD_DMPZ (found by valgrind).
//
// Revision 1.23  2012/10/15 12:38:12  abbott
// Corrected signature for ContentFreeFactorLoop.
//
// Revision 1.22  2012/10/05 09:32:12  abbott
// Changed myExponents into myMultiplicities.
//
// Revision 1.21  2012/05/28 09:18:20  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.20  2012/05/22 10:02:37  abbott
// Removed IsGCDDomain; substituted by IsTrueGCDDomain.
// Added IsFractionFieldOfGCDDomain.
//
// Revision 1.19  2012/05/04 20:02:20  abbott
// Added brackets as suggested by newish version of g++.
//
// Revision 1.18  2012/04/27 15:10:00  abbott
// Cleaned ConvertDUPFFToRingElem; inproved factor fn.
//
// Revision 1.17  2012/04/18 14:26:29  abbott
// Corrected an incorrect assertion (found via anna.cocoa5 test).
//
// Revision 1.16  2012/04/15 20:14:33  abbott
// Added several fns:
//   temporary hacks GCD_DUPFF and GCD_DMPZ to use the old CoCoA-4 code
//   ContentFreeFactor.
//
// Revision 1.15  2012/02/24 13:47:29  abbott
// Temporary check-in of a version which compiles (but new fns DO NOT YET WORK).
//
// Revision 1.14  2012/02/10 10:29:07  bigatti
// -- changed RingZ.H, RingQ.H --> RingZZ.H, RingQQ.H
//
// Revision 1.13  2012/02/08 15:14:30  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.12  2011/08/14 15:52:16  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.11  2011/06/27 12:56:37  bigatti
// -- just a comment Anna --> AMB
//
// Revision 1.10  2009/12/03 17:41:20  abbott
// Added include of cstdlib for malloc.
//
// Revision 1.9  2009/10/26 15:43:25  bigatti
// -- added CopyFromMPZ in ZZ ctor
//
// Revision 1.8  2009/07/24 14:22:40  abbott
// Changed interface to factorizer: new function name (now "factor") with new signature.
// Did some cleaning too.
//
// Revision 1.7  2009/07/02 16:34:32  abbott
// Minor change to keep the compiler quiet.
//
// Revision 1.6  2009/05/14 09:20:03  abbott
// Changed a comment.
//
// Revision 1.5  2008/10/07 15:46:10  abbott
// Added missing log/history lines at the ends of a few files.
//
//
