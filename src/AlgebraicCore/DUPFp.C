#include "CoCoA/DUPFp.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/error.H"
#include "CoCoA/utils.H"
// Next is for STOPGAP impl of factor
#include "TmpFactorDir/DUPFFfactor.h"

#include <algorithm>
using std::max;
#include <iostream>
#include <numeric>
using std::inner_product;
#include <vector>
using std::vector;

namespace CoCoA
{

  DUPFp::DUPFp(long maxdeg, const SmallFpImpl& arith):
      myArith(arith)
  { myCoeffs.reserve(maxdeg+1); }


  //inline void swap(DUPFp& a, DUPFp& b) { std::swap(a.myCoeffs, b.myCoeffs); }


  void AssignZero(DUPFp& f)
  {
    f.myCoeffs.clear();
  }


  void AssignOne(DUPFp& f)
  {
    f.myCoeffs.clear();
    f.myCoeffs.push_back(1);
  }


  bool IsZero(const DUPFp& f)
  {
    return f.myCoeffs.empty();
  }


  long deg(const DUPFp& f) // deg(0) = -1
  {
    return len(f.myCoeffs)-1;
  }


  void FixDeg(DUPFp& f)
  {
    long d = len(f.myCoeffs);
    while (--d >= 0 && f.myCoeffs[d] == 0) {}
    f.myCoeffs.resize(d+1);
  }


  SmallFpImpl::value_t LC(const DUPFp& f)
  {
    if (f.myCoeffs.empty()) return 0;
    return f.myCoeffs.back();
  }


  void MakeMonic(DUPFp& f)
  {
    if (IsZero(f)) CoCoA_ERROR(ERR::DivByZero, "MakeMonic");
    f /= LC(f);
  }

  DUPFp monic(const DUPFp& f)
  {
    if (IsZero(f)) CoCoA_ERROR(ERR::DivByZero, "monic");
    return f/LC(f);
  }


  DUPFp operator*(const DUPFp& f, SmallFpImpl::value_t c)
  {
    const SmallFpImpl& ModP = f.myArith;
    CoCoA_ASSERT(c == ModP.myReduce(c));
    if (c == 0) return DUPFp(0, ModP);
    if (c == 1) return f;
    const long d = deg(f);
    DUPFp ans(d, ModP);
    ans.myCoeffs.resize(d+1);
    for (long i=0; i <= d; ++i)
      ans.myCoeffs[i] = ModP.myMul(f.myCoeffs[i], c);
    return ans;
  }

  DUPFp& operator*=(DUPFp& f, SmallFpImpl::value_t c)
  {
    const SmallFpImpl& ModP = f.myArith;
    CoCoA_ASSERT(c == ModP.myReduce(c));
    if (c == 0)
    {
      f.myCoeffs.clear();
      return f;
    }
    if (c == 1) return f;
    const long N = len(f.myCoeffs);
    for (long i=0; i < N; ++i)
      f.myCoeffs[i] = ModP.myMul(f.myCoeffs[i], c);
    return f;
  }


  DUPFp operator/(const DUPFp& f, SmallFpImpl::value_t c)
  {
    const SmallFpImpl& ModP = f.myArith;
    CoCoA_ASSERT(c == ModP.myReduce(c));
    if (c == 0)
      CoCoA_ERROR(ERR::DivByZero, "op/(DUPFp,n)");
    if (c == 1) return f;
    return f * ModP.myDiv(1,c);
  }


  DUPFp& operator/=(DUPFp& f, SmallFpImpl::value_t c)
  {
    const SmallFpImpl& ModP = f.myArith;
    CoCoA_ASSERT(c == ModP.myReduce(c));
    if (c == 0)
      CoCoA_ERROR(ERR::DivByZero, "op/=(DUPFp,n)");

    if (c == 1) return f;
///???    f *= ModP.myDiv(1, c);
    const SmallFpImpl::value_t crecip = ModP.myDiv(1, c);
    const long N = len(f.myCoeffs);
    for (long i=0; i < N; ++i)
      f.myCoeffs[i] = ModP.myMul(f.myCoeffs[i], crecip);
    return f;
  }


  void add(DUPFp& lhs, const DUPFp& f, const DUPFp& g)
  {
    CoCoA_ASSERT(f.myArith == g.myArith);
    const SmallFpImpl& ModP = f.myArith;
    const long df = deg(f);
    const long dg = deg(g);
    const vector<SmallFpImpl::value_t>& F = f.myCoeffs;
    const vector<SmallFpImpl::value_t>& G = g.myCoeffs;
    // Compute deg of ans in d...
    long d = max(df, dg);
    if (df == dg)
      while (d >= 0 && ModP.myAdd(F[d], G[d]) == 0)
        --d;

    vector<SmallFpImpl::value_t>& H = lhs.myCoeffs;
    H.resize(d+1);
    if (d < 0) return;
    for (; d>df; --d) H[d] = G[d];
    for (; d>dg; --d) H[d] = F[d];
    for (; d>=0; --d)
      H[d] = ModP.myAdd(F[d], G[d]);
  }


  DUPFp operator+(const DUPFp& f, const DUPFp& g)
  {
    DUPFp ans(0, f.myArith);
    add(ans, f, g);
    return ans;
  }


  void sub(DUPFp& lhs, const DUPFp& f, const DUPFp& g)
  {
    CoCoA_ASSERT(f.myArith == g.myArith);
    const SmallFpImpl& ModP = f.myArith;
    const long df = deg(f);
    const long dg = deg(g);
    const vector<SmallFpImpl::value_t>& F = f.myCoeffs;
    const vector<SmallFpImpl::value_t>& G = g.myCoeffs;
    // Compute deg of ans in d...
    long d = max(df, dg);
    if (df == dg)
      while (d >= 0 && F[d] == G[d])
        --d;

    vector<SmallFpImpl::value_t>& H = lhs.myCoeffs;
    H.resize(d+1);
    if (d < 0) return;
    for (; d>df; --d) H[d] = ModP.myNegate(G[d]);
    for (; d>dg; --d) H[d] = F[d];
    for (; d>=0; --d)
      H[d] = ModP.mySub(F[d], G[d]);
  }


  DUPFp operator-(const DUPFp& f, const DUPFp& g)
  {
    DUPFp ans(0, f.myArith);
    sub(ans, f, g);
    return ans;
  }


  void mul(DUPFp& lhs, const DUPFp& f, const DUPFp& g)
  {
    typedef SmallFpImpl::value_t FpElem;
    CoCoA_ASSERT(&lhs != &f);
    CoCoA_ASSERT(&lhs != &g);
    CoCoA_ASSERT(f.myArith == g.myArith);
    const SmallFpImpl& ModP = f.myArith;
    const long df = deg(f);
    const long dg = deg(g);
    const vector<FpElem>& F = f.myCoeffs;
    const vector<FpElem>& G = g.myCoeffs;

    if (IsZero(f) || IsZero(g))
    {
      AssignZero(lhs);
      return;
    }
    if (df < dg) { mul(lhs, g, f); return; }
    // Compute deg of ans in d...
    const long d = df + dg;
    vector<FpElem>& H = lhs.myCoeffs;
    H.clear();
    H.resize(d+1);

    const long MaxIters = ModP.myMaxIters();

    // This impl is about 50% slower than the one after it
//     long count = 0;
//     long pause = 2*MaxIters;
//     for (long i=0; i <= dg; ++i)
//     {
//       if (G[i] == 0) continue;
//       ++count;
//       for (long j=0; j <= df; ++j)
//         H[i+j] += F[j]*G[i];
//       if (count != pause) continue;
//       pause += MaxIters;
//       const long stop = df-MaxIters;
//       for (long j=1; j <= stop; ++j)
//         H[i+j] = ModP.myReduceQ(H[i+j]);
//     }
//     for (long i=0; i <=d; ++i)
//       H[i] = ModP.myReduce(H[i]);

    const long BlockSize = (dg<=2*MaxIters) ? 0 : MaxIters;
    for (long k=0; k <= d; ++k)
    {
      const FpElem *flast, *fi, *gj;
      if (k >= df) flast = &F[df]; else flast = &F[k];
      if (k >= dg) { fi = &F[k-dg]; gj = &G[dg]; }
      else         { fi = &F[0];    gj = &G[k]; }
    
      FpElem sum = 0;
      if (BlockSize > 0)
        for (const FpElem* pause=fi+2*BlockSize; pause <= flast; pause += BlockSize)
        {
          for (; fi < pause; ++fi, --gj)
            sum += *fi * *gj;
          sum = ModP.myHalfNormalize(sum);
        }
      for (; fi <= flast; ++fi, --gj)
        sum += *fi * *gj;
      H[k] = ModP.myNormalize(sum);
    }
  }

  DUPFp operator*(const DUPFp& f, const DUPFp& g)
  {
    DUPFp ans(0, f.myArith);
    mul(ans, f, g);
    return ans;
  }


/***************************************************************************/
/* This is an "in place" squaring function.                                */
/* The multiplication routine above cannot do "in place" squaring, and     */
/* while this routine can easily be modified to do general multiplication  */
/* it is noticeably slower than the routine above.                         */


  DUPFp square(const DUPFp& f)
  {
    CoCoA_ERROR(ERR::NYI, "square for DUPFp");
return f*f;
//     if (IsZero(f)) return DUPFp(0, f.myArith);
//     const long df = deg(f);
//     DUPFp ans(2*df, f.myArith);
//     vector<FpElem>& G = ans.myCoeffs;
//     G.resize(2*df);

//     const long BlockSize = (dg<=2*MaxIters) ? 0 : MaxIters;
//     for (long k=df-1; k >= 0; --k)
//     {
//       const FpElem *fi, *gj;
//       const SmallFpImpl::value_* flast = &F[df];
//       const SmallFpImpl::value_* fi = &F[0];
//       if (2*k < df) flast = &F[2*k];
//       else fi = &F[2*k-df];

//       const FpElem* fj = flast;
    
//       FpElem sum_even = 0;
//       FpElem sum_odd = 0;
//       if (BlockSize > 0)
//         for (const FpElem* pause=fi+2*BlockSize; pause <= flast; pause += BlockSize)
//         {
//           for (; fi < pause; ++fi, --fj)
//           {
//             sum_even += *fi * *fj;
//             sum_odd += *fi * *fj;
//           }
//           sum_even = ModP.myReduceQ(sum_even);
//           sum_odd = ModP.myReduceQ(sum_odd);
//         }
//       for (; fi <= flast; ++fi, --fj)
//       {
//         sum_even += *fi * *gj;
//         sum_odd += *fi * *fj;
//       }
//       H[k] = ModP.myReduce(sum);
//     }


  }


  void power_loop(DUPFp& ans, const DUPFp& base, long exp)
  {
// BUG BUG BUG comment out temporarily to allow compilation
//     if (exp == 1) { ans = base; return; }
//     power_loop(ans, base, exp/2); // integer division!!
//     square(ans);
//     if ((exp&1) == 0) return;
//     ans *= base;
  }

  DUPFp power(const DUPFp& base, long exp)
  {
    if (exp < 0) CoCoA_ERROR(ERR::NegExp, "power(DUPFp,n)");
    if (IsZero(base))
    {
      DUPFp ans(0, base.myArith);
      if (exp == 0) AssignOne(ans);
      return ans;
    }
    DUPFp ans(exp*deg(base), base.myArith);
    if (exp == 0) { AssignOne(ans); return ans; }
    power_loop(ans, base, exp);
    return ans;
  }


  // f += c*x^exp*g
  void ShiftAdd(DUPFp& f, const DUPFp& g, SmallFpImpl::value_t c, long DegShift)
  {
    typedef SmallFpImpl::value_t FpElem;
    CoCoA_ASSERT(f.myArith == g.myArith);
    if (c == 0) return;
    const SmallFpImpl& ModP = f.myArith;
    const long df = deg(f);
    const long dg = deg(g);
    vector<FpElem>& F = f.myCoeffs;
    const vector<FpElem>& G = g.myCoeffs;
    if (dg+DegShift > df) F.resize(dg+DegShift+1); // fills with 0

//     FpElem* Fj= &F[dg+exp];
//     const FpElem* const G0 = &G[0];
//     for (const FpElem* Gi = &G[dg]; Gi >= G0; --Gi, --Fj)
//       ModP.myIsZeroAddMul(*Fj, c, *Gi);
/////    std::clog<<"ShiftAdd BEFORE LOOP DegShift="<<DegShift<<std::endl;
/////    std::clog<<"ShiftAdd f="<<f<<std::endl;
/////    std::clog<<"ShiftAdd g="<<g<<std::endl;
    for (long i=dg; i >= 0; --i)
      ModP.myIsZeroAddMul(F[i+DegShift], c, G[i]); // F[i+DegShift] = ModP.myNormalize(F[i+exp] + c*G[i]);
/////    std::clog<<"ShiftAdd f="<<f<<std::endl;

    FixDeg(f);
  }


  DUPFp operator/(const DUPFp& num, const DUPFp& den)
  {
    CoCoA_ASSERT(num.myArith == den.myArith);
    if (IsZero(den)) CoCoA_ERROR(ERR::DivByZero, "operator/ for DUPFp");
    if (IsZero(num)) return DUPFp(0, num.myArith);//??? return num;???
    const long dnum = deg(num);
    const long dden = deg(den);
    DUPFp quot(dnum-dden, num.myArith);
    DUPFp rem(dden-1, num.myArith);
    QuoRem(quot, rem, num, den);
    return quot;
  }

  DUPFp operator%(const DUPFp& num, const DUPFp& den)
  {
    CoCoA_ASSERT(num.myArith == den.myArith);
    if (IsZero(den)) CoCoA_ERROR(ERR::DivByZero, "operator/ for DUPFp");
    if (IsZero(num)) return DUPFp(0, num.myArith);
    const long dnum = deg(num);
    const long dden = deg(den);
    DUPFp quot(dnum-dden, num.myArith);
    DUPFp rem(dden-1, num.myArith);
    QuoRem(quot, rem, num, den);
    return rem;
  }


  void QuoRem(DUPFp& q, DUPFp& r, const DUPFp& num, const DUPFp& den)
  {
    if (&q == &r || &q == &num || &q == &den ||
        &r == &num || &r == &den)
      CoCoA_ERROR(ERR::BadArg, "QuoRem args alias");
    const SmallFpImpl& ModP = num.myArith;
    CoCoA_ASSERT(den.myArith == ModP);
    CoCoA_ASSERT(q.myArith == ModP);
    CoCoA_ASSERT(r.myArith == ModP);
    if (IsZero(den)) CoCoA_ERROR(ERR::DivByZero, "QuoRem for DUPFp");
    if (IsZero(num)) { AssignZero(q); AssignZero(r); return; }

    CoCoA_ERROR(ERR::NYI, "QuoRem");

// CHECK FOR ALIASING!!!
//   if (quot == rem || quot == den || rem == den)
//   {
//     JERROR(JERROR_DIV4_ARGS);
//     return;
//   }
    const long dnum = deg(num);
    const long dden = deg(den);

    if (dnum < dden) /* quotient is zero, remainder = num */
    {
      AssignZero(q);
      r = num;
      return;
    }

    typedef SmallFpImpl::value_t FpElem;
    vector<FpElem>& Q = q.myCoeffs;
    vector<FpElem>& R = r.myCoeffs;
///???    const vector<FpElem>& N = num.myCoeffs;
    const vector<FpElem>& D = den.myCoeffs;

    long dquot = dnum-dden;
    if (dquot > deg(q)) Q.resize(dquot+1);
    r = num;


    const FpElem lcdrecip = ModP.myDiv(1, D[dden]);
    const long MaxIters = ModP.myMaxIters();
    const FpElem p = ModP.myModulus();
    long count = 0;
///???  if (dden < k) k = dquot+1; /* disable "shifting" if den has low degree */

    for(long dtmp = dnum; dquot >= 0; --dquot, --dtmp)
    {
      R[dtmp] = ModP.myReduce(R[dtmp]);
      if (R[dtmp] == 0) { Q[dquot] = 0; continue; }
      Q[dquot] =  ModP.myMul(R[dtmp], lcdrecip);
      const FpElem qq = p - Q[dquot];
      for (long i=0; i < dden; ++i)
        R[i+dquot] += qq*D[i]; // NOT normalized!!
      if (++count < MaxIters) continue;
      count = 0;
      for (long i=MaxIters; i < dden; ++i)
        R[i+dquot] = ModP.myHalfNormalize(R[i+dquot]);
    }

    for(long i=0; i < dden; ++i) R[i] = ModP.myNormalize(R[i]);

    FixDeg(r);
  }


  DUPFp deriv(const DUPFp& f)
  {
    long d = deg(f);
    const SmallFpImpl& ModP = f.myArith;
    const long p = ModP.myModulus();
    while (d > 0)
    {
      if ((f.myCoeffs[d] != 0) && (d%p != 0)) break;
      --d;
    }
    if (d == 0) return DUPFp(0, ModP);
    DUPFp fprime(d-1, ModP);
    for (long i=1; i <= d; ++i)
      fprime.myCoeffs[i-1] = ModP.myMul(ModP.myReduce(i), f.myCoeffs[i]);
    return fprime;
  }


  // anon namespace?
  // OVERWRITES ARGS!!!
  DUPFp EuclidAlgm(DUPFp& f, DUPFp& g)
  {
    typedef SmallFpImpl::value_t FpElem;
    CoCoA_ASSERT(f.myArith == g.myArith);
    if (deg(f) < deg(g)) return EuclidAlgm(g, f);

    const SmallFpImpl& ModP = f.myArith;
    while (!IsZero(g))
    {
      while (deg(f) >= deg(g))
      {
        const FpElem q = ModP.myNegate(ModP.myDiv(LC(f), LC(g)));
        ShiftAdd(f, g, q, deg(f)-deg(g));
      }
      swap(f, g);
    }
    return f;
  }

  DUPFp gcd(const DUPFp& f, const DUPFp& g)
  {
    if (f.myArith != g.myArith) CoCoA_ERROR(ERR::MixedRings, "gcd for DUPFp");
    if (IsZero(f)) return monic(g);
    if (IsZero(g)) return monic(f);
    if (deg(f) == 0) return monic(f);
    if (deg(g) == 0) return monic(g);
    DUPFp fcopy = f;
    DUPFp gcopy = g;
    return monic(EuclidAlgm(fcopy, gcopy));
  }

  DUPFp ExtEuclidAlgm(DUPFp& cf, DUPFp& cg, DUPFp& f, DUPFp& g)
  {
    typedef SmallFpImpl::value_t FpElem;
    CoCoA_ASSERT(f.myArith == g.myArith);
    if (deg(f) < deg(g)) return ExtEuclidAlgm(cg, cf, g, f);

    const SmallFpImpl& ModP = f.myArith;
    DUPFp m11(deg(g), ModP); AssignOne(m11);
    DUPFp m12(deg(f), ModP);
    DUPFp m21(deg(g), ModP);
    DUPFp m22(deg(f), ModP); AssignOne(m22);
    while (!IsZero(g))
    {
      while (deg(f) >= deg(g))
      {
        const FpElem q = ModP.myNegate(ModP.myDiv(LC(f), LC(g)));
        const long DegShift = deg(f)-deg(g);
/////        std::clog<<"EEA: q="<<q<<" shift="<<DegShift<<std::endl;
        ShiftAdd(f, g, q, DegShift);
/////        std::clog<<"BEFORE UPDATE:\n";
/////        std::clog<<"m11="<<m11<<std::endl;
/////        std::clog<<"m12="<<m12<<std::endl;
/////        std::clog<<"m21="<<m21<<std::endl;
/////        std::clog<<"m22="<<m22<<std::endl;
        ShiftAdd(m11, m21, q, DegShift);
        ShiftAdd(m12, m22, q, DegShift);
/////        std::clog<<"AFTER UPDATE:\n";
/////        std::clog<<"m11="<<m11<<std::endl;
/////        std::clog<<"m12="<<m12<<std::endl;
/////        std::clog<<"m21="<<m21<<std::endl;
/////        std::clog<<"m22="<<m22<<std::endl;
      }
      swap(f, g);
      swap(m11, m21);
      swap(m12, m22);
    }
    swap(cf, m11); // equiv. cf = m11;
    swap(cg, m12); // equiv. cg = m12;
    return f;
  }

  DUPFp exgcd(DUPFp& cf, DUPFp& cg, const DUPFp& f, const DUPFp& g)
  {
    if (f.myArith != g.myArith) CoCoA_ERROR(ERR::MixedRings, "exgcd for DUPFp");
    if (cf.myArith != cg.myArith) CoCoA_ERROR(ERR::MixedRings, "exgcd for DUPFp");
    if (f.myArith != cf.myArith) CoCoA_ERROR(ERR::MixedRings, "exgcd for DUPFp");
    if (IsZero(f) && IsZero(g))  { AssignZero(cf); AssignZero(cg); return f; } // all zero!
    if (IsZero(f) || deg(g) == 0) { AssignZero(cf); AssignOne(cg); cg /= LC(g); return monic(g); }
    if (IsZero(g) || deg(f) == 0) { AssignOne(cf); AssignZero(cg); cf /= LC(f); return monic(f); }

    DUPFp fcopy(f); // work on copies of f & g!
    DUPFp gcopy(g);
    DUPFp h = ExtEuclidAlgm(cf, cg, fcopy, gcopy);
    cf /= LC(h);
    cg /= LC(h);
    return monic(h);
  }


  SmallFpImpl::value_t eval(const DUPFp& f, SmallFpImpl::value_t x)
  {
    typedef SmallFpImpl::value_t FpElem;
    if (IsZero(f)) return 0;
    if (deg(f) == 0) return LC(f);
    const vector<FpElem>& F = f.myCoeffs;
    if (x == 0) return F[0];
    const SmallFpImpl& ModP = f.myArith;
    FpElem ans = 0;
    for (long i=deg(f); i >= 0; --i)
      ans = ModP.myAdd(F[i], ModP.myMul(ans, x));
    return ans;
  }


  bool operator==(const DUPFp& f, const DUPFp& g)
  {
    if (f.myArith != g.myArith) CoCoA_ERROR(ERR::MixedRings, "op== for DUPFp");
    if (deg(f) != deg(g)) return false;
    return f.myCoeffs == g.myCoeffs;  // let std::vector do the work :-)
  }

  bool operator!=(const DUPFp& f, const DUPFp& g)
  {
    return !(f == g);
  }


  std::ostream& operator<<(std::ostream& out, const DUPFp& f)
  {
    if (IsZero(f)) { out << "0"; return out; }
    const long d = deg(f);
    const long p = f.myArith.myModulus();
    out << "DUPFp(char=" << p << ", coeffs=[" << f.myCoeffs[d];  ///???? BUG????  use myExport???
    for (long i=d-1; i >= 0; --i)
      out << ", " << f.myCoeffs[i];
    out << "])";
    return out;
  }


  DUPFp MulMod(const DUPFp& f, const DUPFp& g, const DUPFp& m)
  {
    return (f*g)%m;
  }

  DUPFp PowerMod(const DUPFp& f, long exp, const DUPFp& m)
  {
    CoCoA_ASSERT(exp >= 0);
    CoCoA_ASSERT(f.myArith == m.myArith);
    CoCoA_ASSERT(deg(m) > 0);
    CoCoA_ASSERT(deg(f) < deg(m));
    const SmallFpImpl& ModP = f.myArith;
    DUPFp ans(deg(m), ModP);
    if (exp == 0) { AssignOne(ans); return ans; }
    if (exp == 1) return f;
    ans = PowerMod(f, exp/2, m);
    ans = MulMod(ans, ans, m);
    if ((exp&1) == 0) return ans;
    return MulMod(ans, f, m);
  }


  bool IsSqfr(const DUPFp& f)
  {
    if (IsZero(f)) return false;
    if (deg(f) <= 1) return true;
    return (deg(gcd(f, deriv(f))) == 0);
  }


  DUPFp PthRoot(const DUPFp& f)
  {
    typedef SmallFpImpl::value_t FpElem;
    if (IsZero(f) || deg(f) == 0) return f;
    const SmallFpImpl& ModP = f.myArith;
    const long p = ModP.myModulus();
    const long degf = deg(f);
    CoCoA_ASSERT(degf%p == 0);
    const vector<FpElem>& F = f.myCoeffs;
#ifdef CoCoA_DEBUG
    for (long i=0; i <= degf; ++i)
      CoCoA_ASSERT(i%p == 0 || F[i] == 0);
#endif
    DUPFp ans(degf/p, ModP);
    vector<FpElem>& A = ans.myCoeffs;
    A.resize(1+degf/p);
    for (long i=0; p*i <= degf; ++i)
      A[i] = F[p*i];
    return ans;
  }


  factorization<DUPFp> SqfrDecomp(DUPFp f)
  {
    if (IsZero(f)) CoCoA_ERROR(ERR::NotNonZero, "SqfrDecomp");
    CoCoA_ASSERT(IsSqfr(f));
    vector<DUPFp> factors;
    vector<long> multiplicities;
    const SmallFpImpl& ModP = f.myArith;
    const long p = ModP.myModulus();
    long q = 1; // q will always be a power of p.

    while (deg(f) > 0)
    {
      DUPFp g = gcd(f, deriv(f));
      DUPFp radical = f/g;
      swap(f, g);
      long pwr = 0;
      while (deg(radical) > 0)
      {
        ++pwr;
        if (pwr%p == 0) { f = f/radical; continue; }
        DUPFp NewRadical = gcd(f, radical);
        if (deg(NewRadical) < deg(radical))
        {
          factors.push_back(radical/NewRadical);
          multiplicities.push_back(pwr*q);
          swap(radical, NewRadical);
        }
        f = f/radical;
      }
      q *= p;
      f = PthRoot(f);
    }
    DUPFp RemainingFactor(0, ModP); AssignOne(RemainingFactor); RemainingFactor *= LC(f);
    return factorization<DUPFp>(factors, multiplicities, RemainingFactor);
  }

  factorization<DUPFp> DistinctDegFactor(DUPFp f)
  {
    if (IsZero(f)) CoCoA_ERROR(ERR::NotNonZero, "DistinctDegFactor");
    CoCoA_ASSERT(IsSqfr(f));
    vector<DUPFp> factors;
    const SmallFpImpl& ModP = f.myArith;
    const long p = ModP.myModulus();
    DUPFp x(1, ModP); x.myCoeffs.resize(1); x.myCoeffs[1]=1;
    DUPFp xpower(deg(f)-1, ModP);
    xpower = x;
    for (long d=1; d <= deg(f)/2; ++d)
    {
      xpower = PowerMod(xpower, p, f);
///???      xpower -= x;
      const DUPFp g = gcd(xpower-x, f);
///???      xpower += x;
      factors.push_back(g);
      if (deg(g) == 0) continue;
      f = f/g;
      xpower = xpower%f;
    }
    if (deg(f) > 0) factors.push_back(monic(f));
    DUPFp RemainingFactor(0, ModP); AssignOne(RemainingFactor); RemainingFactor *= LC(f);
    const vector<long> mult(len(factors), 1);
    return factorization<DUPFp>(factors, mult, RemainingFactor);
  }


  DUPFp ConvertToDUPFp(const SmallFpImpl& ModP, ConstRefRingElem f)
  {
    if (!IsSparsePolyRing(owner(f))) CoCoA_ERROR(ERR::NotSparsePolyRing, "ConvertToDUPFp");
    CoCoA_ASSERT(UnivariateIndetIndex(f) >= 0);
//????    if (UnivariateIndetIndex(f) < 0) CoCoA_ERROR("not univariate","ConvertToDUPFp");
    const long degf = StdDeg(f);
    DUPFp ans(degf, ModP);
    ans.myCoeffs.resize(degf+1);
    for (SparsePolyIter it = BeginIter(f); !IsEnded(it); ++it)
    {
      ans.myCoeffs[StdDeg(PP(it))] = ModP.myReduce(ConvertTo<BigInt>(coeff(it)));
    }
    return ans;
  }

  RingElem ConvertFromDUPFp(ConstRefRingElem x, const DUPFp& f)
  {
    if (!IsPolyRing(owner(x))) CoCoA_ERROR(ERR::NotPolyRing, "ConvertFromDUPFp");
    if (!IsIndet(x)) CoCoA_ERROR("x not indet", "ConvertFromDUPFp");
    const PolyRing P = owner(x);
    RingElem ans(P);
    if (IsZero(f)) return ans;
    const long degf = deg(f);
    for (long i=0; i <= degf; ++i)
      if (f.myCoeffs[i] != 0)
        ans += f.myCoeffs[i]*power(x,i);
    return ans;
  }


  // ******* STOPGAP IMPL *********
  factorization<DUPFp> factor(const DUPFp& f)
  {
    const SmallFpImpl& ModP = f.myArith;
    const long p = ModP.myModulus();
    const long d = deg(f);

    FF Fp = FFctor(p);
    FFselect(Fp);
    DUPFF F = DUPFFnew(d);
    { FFelem* CoeffVec = F->coeffs; for (long i=0; i<=d; ++i) CoeffVec[i] = f.myCoeffs[i]; }
    F->deg = d;
    DUPFFlist facs = DUPFFfactor(F);
    DUPFFfree(F);
    const long n = DUPFFlist_length(facs);
//std::clog<<"factor: n="<<n<<std::endl;
    vector<DUPFp> irreds; irreds.reserve(n);
    vector<long> multiplicity; multiplicity.reserve(n);
    for (DUPFFlist curr=facs; curr != 0/*nullptr*/; curr = curr->next)
    {
      const long DegFac = DUPFFdeg(curr->poly);
      DUPFp fac(DegFac, ModP); fac.myCoeffs.resize(DegFac+1);
      for (long i=0; i <= DegFac; ++i) fac.myCoeffs[i] = ModP.myReduce(curr->poly->coeffs[i]);
//std::clog<<"factor: fac="<<fac<<std::endl;
      irreds.push_back(monic(fac));
      multiplicity.push_back(curr->deg);
    }
    DUPFFlist_dtor(facs);
///???    FFselect(0);
    FFdtor(Fp);
    DUPFp scalar(0, ModP); scalar.myCoeffs.resize(1); scalar.myCoeffs[0] = LC(f);
    return factorization<DUPFp>(irreds, multiplicity, scalar);
  }

} // end of namespace CoCoA
