//   Copyright (c) 2010-2014 Anna Maria Bigatti, 2010 Giovanni Lagorio (lagorio@disi.unige.it)
//
//   This file is part of the source of CoCoALib, the CoCoA Library.
//
//   CoCoALib is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   CoCoALib is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with CoCoALib.  If not, see <http://www.gnu.org/licenses/>.

#include "CoCoALibSupplement.H"
#include "CoCoA/library.H"

#include <sstream>  // for ErrorMessage
#include <string>
//#include <vector> // already included by library.H->io.H
using std::vector;

namespace CoCoA {


  namespace // anonymous namespace
  {
    // for GetRow_forC5 and NBM_forC5  -- indexed from 0 as in CoCoALib
    std::vector<RingElem> GetRow_long(ConstMatrixView M, long i)
    {
      vector<RingElem> v;
      for (long j=0; j<NumCols(M); ++j)  v.push_back(M(i,j));
      return v;
    }

    // for NBM_forC5
    std::vector<std::vector<RingElem> > GetRows(ConstMatrixView M)
    {
      vector<vector<RingElem> > v;
      for (long i=0; i<NumRows(M); ++i)  v.push_back(GetRow_long(M,i));
      return v;
    }

    // for NBM_forC5 and QuotientBasis_forC5
    std::vector<RingElem> VectorRingElem(const SparsePolyRing& P, std::vector<PPMonoidElem> v)
    {
      const long N = len(v);
      vector<RingElem> ans; ans.reserve(N);
      for (long i=0; i < N; ++i)
        ans.push_back(monomial(P, 1, v[i]));
      return ans;
    }
  } // anonymous namespace
  

  std::vector<long> VectorLong(const std::vector<BigInt>& BigIntVec, const std::string& FuncName)
  {
    vector<long> v;
    long tmp;
    for (long i=0; i<len(BigIntVec); ++i)
      if (IsConvertible(tmp, BigIntVec[i])) v.push_back(tmp);
      else CoCoA_ERROR(ERR::BadSymbolRange, FuncName);
    return v;
  }
  

  std::vector<long> VectorLongDecr1(const std::vector<BigInt>& BigIntVec, const ERR::ID& ErrID, const std::string& FuncName)
  {
    vector<long> v;
    long tmp;
    for (long i=0; i<len(BigIntVec); ++i)
      if (IsConvertible(tmp, BigIntVec[i])) v.push_back(tmp-1);
      else CoCoA_ERROR(ErrID, FuncName);
    return v;
  }
  

  bool IsTerm_forC5(ConstRefRingElem f)
  {
    return IsMonomial(f)&&IsOne(LC(f));
  }

  
  RingElem LT_forC5(ConstRefRingElem f)
  {
    return LPP_forC5(f);
  }  ///< NB result belongs to ring (owner(f))


  ModuleElem LT_forC5(const ModuleElem& f)
  {
    //    return LPP(f)*gens(owner(f))[LPos(f)];
    return LPP_forC5(f)*gens(owner(f))[LPosn(f)];
  }  ///< NB result belongs to module


  long LPosn_forC5(const ModuleElem& f)
  { return LPosn(f)+1; }


  long FirstNonZeroPosn_forC5(const ModuleElem& f)
  { return FirstNonZeroPosn(f)+1; }
  

  RingElem LPP_forC5(ConstRefRingElem f)
  {
    ring R = owner(f);
    if (!IsPolyRing(R)) CoCoA_ERROR(ERR::NotElemPolyRing, "LPP_forC5");
    if (IsSparsePolyRing(R)) return monomial(R, 1, LPP(f));
    return IndetPower(R, 0, deg(f)); // univariate case
  }


  RingElem LPP_forC5(const ModuleElem& f)
  {
    const ring R = RingOf(owner(f));
    if (!IsPolyRing(R)) CoCoA_ERROR(ERR::NotElemPolyRing, "LPP_forC5");
    if (IsSparsePolyRing(R)) return monomial(R, 1, LPP(f));
    //    PolyRing Rx = AsPolyRing(R); // univariate case
    //    return IndetPower(Rx, 0, deg(f));
    CoCoA_ERROR(ERR::NYI, "LPP_forC5 univariate module");
    return IndetPower(R, 0, 1);
  }


  RingElem LM_forC5(ConstRefRingElem f)
  {
    ring R = owner(f);
    if (!IsPolyRing(R)) CoCoA_ERROR(ERR::NotElemPolyRing, "LM_forC5");
    if (IsSparsePolyRing(R)) return monomial(R, LC(f), LPP(f));
    return CoeffEmbeddingHom(R)(LC(f))*IndetPower(R, 0, deg(f));    
  }


  ModuleElem LM_forC5(const ModuleElem& f)
  {
    long LPosnf = LPosn(f);
    return LM_forC5(f[LPosnf])*gens(owner(f))[LPosnf];
  }  ///< NB result belongs to module


  namespace
  {
    RingElem CoeffOfTermSparse(ConstRefRingElem f, ConstRefPPMonoidElem pp)
    {
      for (SparsePolyIter itf=BeginIter(f); !IsEnded(itf); ++itf)
        if (PP(itf) == pp) return coeff(itf);
      return zero(CoeffRing(owner(f)));
    }
    
    std::vector<RingElem> CoefficientsDense(ConstRefRingElem f)
    {
      vector<RingElem> v;
      const DenseUPolyRing P = owner(f);
      for (long i=P->myDegPlus1(raw(f))-1; i>=0; --i)
        if (!IsZero(coeff(f,i))) v.push_back(coeff(f,i));
      return v;
    }
    
    std::vector<RingElem> CoefficientsSparse(ConstRefRingElem f)
    {
      vector<RingElem> v;
      for (SparsePolyIter itf=BeginIter(f); !IsEnded(itf); ++itf)
        v.push_back(coeff(itf));
      return v;
    }
    
    std::vector<RingElem> MonomialsDense(ConstRefRingElem f)
    {
      vector<RingElem> v;
      const DenseUPolyRing P = owner(f);
      for (long i=P->myDegPlus1(raw(f))-1; i>=0; --i)
        if (!IsZero(coeff(f,i))) v.push_back(monomial(P, coeff(f,i), i));
      return v;
    }
    
    std::vector<RingElem> MonomialsSparse(ConstRefRingElem f)
    {
      vector<RingElem> v;
      const SparsePolyRing P = owner(f);
      for (SparsePolyIter itf=BeginIter(f); !IsEnded(itf); ++itf)
        v.push_back(monomial(P, coeff(itf), PP(itf)));
      return v;
    }
    
    std::vector<RingElem> SupportDense(ConstRefRingElem f)
    {
      CoCoA_ERROR("SupportDense is disabled", "SupportDense");
      vector<RingElem> v;
      const DenseUPolyRing P = owner(f);
      for (long i=P->myDegPlus1(raw(f))-1; i>=0; --i)
        if (!IsZero(coeff(f,i))) v.push_back(monomial(P, one(CoeffRing(P)), i));
      return v;
    }

    std::vector<RingElem> SupportSparse(ConstRefRingElem f)
    {
      vector<RingElem> v;
      const SparsePolyRing P = owner(f);
      for (SparsePolyIter itf=BeginIter(f); !IsEnded(itf); ++itf)
        v.push_back(monomial(P, one(CoeffRing(P)), PP(itf)));
      return v;
    }  
  }
  

  RingElem CoeffOfTerm_forC5(ConstRefRingElem f, ConstRefRingElem t)
  {
    ring R = owner(f);
    if (!IsPolyRing(R)) CoCoA_ERROR(ERR::NotElemPolyRing, "CoeffOfTerm_forC5");
    if (!IsMonomial(t)) CoCoA_ERROR(ERR::BadArg, "CoeffOfTerm_forC5");
    if (!IsOne(LC(t))) CoCoA_ERROR(ERR::BadArg, "CoeffOfTerm_forC5");
    if (IsSparsePolyRing(R))
      return CoeffOfTermSparse(f, LPP(t));
    return coeff(f, deg(t));
  }
  
  
  std::vector<RingElem> coefficients_forC5(ConstRefRingElem f)
  {
    ring R = owner(f);
    if (!IsPolyRing(R)) CoCoA_ERROR(ERR::NotElemPolyRing, "coefficients_forC5");
    if (IsSparsePolyRing(R))  return CoefficientsSparse(f);
    return CoefficientsDense(f);
  }
  
  
  std::vector<RingElem> monomials_forC5(ConstRefRingElem f)
  {
    ring R = owner(f);
    if (!IsPolyRing(R)) CoCoA_ERROR(ERR::NotElemPolyRing, "monomials_forC5");
    if (IsSparsePolyRing(R))  return MonomialsSparse(f);
    return MonomialsDense(f);
  }
  
  
  std::vector<RingElem> support_forC5(ConstRefRingElem f)
  {
    ring R = owner(f);
    if (!IsPolyRing(R)) CoCoA_ERROR(ERR::NotElemPolyRing, "support_forC5");
    if (IsSparsePolyRing(R))  return SupportSparse(f);
    return SupportDense(f);
  }


  namespace
  {
    // needs a decent implementation, this one isn't
    RingElem DensePolyRec(const PolyRing& P, long n, long d)
    {
      if (n==1) return IndetPower(P,0,d);
      if (d==1)
      {
        RingElem s=indet(P,0);
        for (long i=1; i<n; ++i)  s+=indet(P,i);
        return s;
      }
      return indet(P,n-1)*DensePolyRec(P,n,d-1) + DensePolyRec(P,n-1,d);
    }    
  }


  RingElem DensePoly_forC5(const ring& R, const BigInt& D)
  {
    CoCoA_ASSERT(IsPolyRing(R));
    long d;
    if (!IsConvertible(d,D))  CoCoA_ERROR(ERR::ArgTooBig, "DensePoly_forC5");
    if (d<0)  CoCoA_ERROR(ERR::NotNonNegative, "DensePoly_forC5");
    if (d==0) return one(R);
    return DensePolyRec(R, NumIndets(R), d);
  }


  std::vector<BigInt> DegreeToVec(degree d)
  {
    vector<BigInt> v;
    for (long i=0; i<GradingDim(d); ++i)  v.push_back(d[i]);
    return v;
  }


  std::vector<BigInt> wdeg_forC5(const ModuleElem& f)
  {
    degree d = wdeg(f);
    vector<BigInt> v;
    for (long i=0; i<GradingDim(d); ++i)  v.push_back(d[i]);
    return v;
  }


  BigInt NextPrime_forC5(BigInt N)
  {
    if (N < 0) CoCoA_ERROR(ERR::NotNonNegative, "NextPrime_forC5");
    long n;
    if (!IsConvertible(n, N)) CoCoA_ERROR(ERR::ArgTooBig, "NextPrime_forC5");
    long ans = NextPrime(n);
    if (ans == 0) CoCoA_ERROR(ERR::ArgTooBig, "NextPrime_forC5");
    return BigInt(ans);
  }

  factorization<BigInt> SmoothFactor_forC5(const BigInt& N, const BigInt& TrialLimit)
  {
    if (N <= 0) CoCoA_ERROR(ERR::NotPositive, "SmoothFactor_forC5");
    long tl;
    if (!IsConvertible(tl, TrialLimit)) CoCoA_ERROR(ERR::ArgTooBig, "SmoothFactor_forC5");
    return SmoothFactor(N, tl);
  }


  BigRat CpuTime_forC5()
  {
    return ConvertTo<BigRat>(CpuTime()); 
//     BigRat res;
//     IsConvertible(res, CpuTime());
//     return res;
  }
  

  //  std::string date_forC5()
  BigInt date_forC5()
  {
//     const static std::string d(__DATE__);
//     const static std::string t(__TIME__);
    // fix day of the week
    //    return "    " + d.substr(0,7) + t + d.substr(6,5);
    long date, time;
    DateTime(date, time);
    return BigInt(date);  //*1000000 + time;
  }
  

  BigInt TimeOfDay_forC5()
  {
    long date, time;
    DateTime(date, time);
    return BigInt(time);
  }
  

  const BigInt random_forC5(const BigInt& lo, const BigInt& hi)
  { return RandomBigInt(GlobalRandomSource(), lo, hi); }
  

  BigInt lcm_forC5(const std::vector<BigInt>& v)
  {  // v is not empty
    BigInt res = v[0];
    for (long i=1; i<len(v); ++i)  res = lcm(res, v[i]);
    return res;
  }

  
  BigInt gcd_forC5(const std::vector<BigInt>& v)
  {  // v is not empty
    BigInt res = v[0];
    for (long i=1; i<len(v); ++i)  res = gcd(res, v[i]);
    return res;
  }
  

  RingElem lcm_forC5(const std::vector<RingElem>& v)
  {  // v is not empty
    RingElem res = v[0];
    for (long i=1; i<len(v); ++i)  res = lcm(res, v[i]);
    return res;
  }
  

  RingElem gcd_forC5(const std::vector<RingElem>& v)
  {  // v is not empty
    RingElem res = v[0];
    for (long i=1; i<len(v); ++i)  res = gcd(res, v[i]);
    return res;
  }
  

  RingElem ContentWRT_forC5(ConstRefRingElem f, const std::vector<RingElem>& v)
  {
    vector<long> indices(len(v));
    for (long i=0; i<len(v); ++i)
      if (!IsIndet(indices[i], v[i]))
        CoCoA_ERROR(ERR::NotIndet, "ContentWRT_forC5");
    return ContentWRT(f, indices);
  }
  

  std::vector<CoeffPP> CoefficientsWRT_forC5(ConstRefRingElem f, const std::vector<RingElem>& v)
  {
    vector<long> indices(len(v));
    for (long i=0; i<len(v); ++i)
      if (!IsIndet(indices[i], v[i]))
        CoCoA_ERROR(ERR::NotIndet, "ContentWRT_forC5");
    return CoefficientsWRT(f, indices);
  }
  

  std::vector<RingElem> CoeffListWRT_forC5(ConstRefRingElem f, ConstRefRingElem x)
  {
    return CoeffVecWRT(f, x);
  }


  std::vector<RingElem> QuotientBasis_forC5(const ideal& I)
  {
    return VectorRingElem(RingOf(I), QuotientBasis(I));
  }


  void NBM_forC5(std::vector<RingElem>& QB, std::vector<RingElem>& BB, std::vector<RingElem>& AV,
                 const SparsePolyRing& P, ConstMatrixView OrigPts, ConstMatrixView OrigTolerance)
  {
    std::vector<vector<RingElem> > pts(GetRows(apply(CanonicalHom(RingQQ(), CoeffRing(P)), OrigPts)));
    if (NumRows(OrigTolerance) != 1) CoCoA_ERROR(ERR::BadMatrixSize, "NBM_forC5");
    vector<RingElem> tol(GetRow_long(OrigTolerance, 0));
    vector<PPMonoidElem> tmpQB0;
    vector<RingElem> BB0;
    vector<RingElem> AV0;
    ApproxPts::NBM(tmpQB0, BB0, AV0, P, pts, tol);
    vector<RingElem> QB0(VectorRingElem(P, tmpQB0));
    swap(QB, QB0);
    swap(BB, BB0);
    swap(AV, AV0);
  }


//   void SOI_forC5(std::vector<RingElem>& QB, std::vector<RingElem>& BB, std::vector<RingElem>& AV,
//                  const SparsePolyRing& P, ConstMatrixView OrigPts, ConstMatrixView OrigTolerance)
//   {
//     std::vector<vector<RingElem> > pts(GetRows(apply(CanonicalHom(RingQQ(), CoeffRing(P)), OrigPts)));
//     if (NumRows(OrigTolerance) != 1) CoCoA_ERROR(ERR::BadMatrixSize, "NBM_forC5");
//     vector<RingElem> tol(GetRow_long(OrigTolerance, 0));
//     vector<PPMonoidElem> tmpQB0;
//     vector<RingElem> BB0;
//     vector<RingElem> AV0;
//     ApproxPts::SOI(tmpQB0, BB0, AV0, P, pts, tol);
//     vector<RingElem> QB0(VectorRingElem(P, tmpQB0));
//     swap(QB, QB0);
//     swap(BB, BB0);
//     swap(AV, AV0);
//   }


  RingElem ClosePassingPoly_forC5(const ring& R, ConstMatrixView OrigPts, ConstMatrixView OrigTolerance)
  {
    std::vector<vector<RingElem> > pts(GetRows(apply(CanonicalHom(RingQQ(), CoeffRing(R)), OrigPts)));
    if (NumRows(OrigTolerance) != 1) CoCoA_ERROR(ERR::BadMatrixSize, "AlmostVanishing_forC5");
    vector<RingElem> tol(GetRow_long(OrigTolerance, 0));
    vector<RingElem> poly0;
    vector<ApproxPts::PointR> NewPointsVec0;
    RingElem MaxTol=zero(RingQQ());
    for (long i=0; i<len(tol); ++i)
      if (MaxTol<tol[i]) MaxTol = tol[i];
    ApproxPts::VanishPoly(NewPointsVec0, poly0, R, pts, tol, MaxTol);
    return poly0[0];
  }
  

  void PreprocessPts_forC5(const std::string& WhichAlgm,
                           std::vector< std::vector<RingElem> >& NewPts,
                           std::vector<long>& weights,
                           ConstMatrixView OrigPts,
                           ConstMatrixView epsilon)
  {
    if (!IsOrderedDomain(RingOf(OrigPts))) CoCoA_ERROR(ERR::NotOrdDom, "PreprocessGridAlgm_forC5");
    if (RingOf(epsilon) != RingOf(OrigPts)) CoCoA_ERROR(ERR::MixedRings, "PreprocessGridAlgm_forC5");
    const long NumPts = NumRows(OrigPts);
    const long dim = NumCols(OrigPts);
    if (NumCols(epsilon) != dim || NumRows(epsilon) != 1) CoCoA_ERROR(ERR::BadMatrixSize, "PreprocessGridAlgm_forC5");
    for (long j=0; j < dim; ++j)
      if (epsilon(0,j) <= 0) CoCoA_ERROR(ERR::BadArg, "PreprocessGridAlgm epsilon");
    vector< vector<RingElem> > pts(NumPts, vector<RingElem>(dim));
    for (long i=0; i < NumPts; ++i)
      for (long j=0; j < dim; ++j)
        pts[i][j] = OrigPts(i,j);
    vector<RingElem>  eps(dim);
    for (long j=0; j < dim; ++j)
      eps[j] = epsilon(0,j);

    if (WhichAlgm == "grid")
      PreprocessPtsGrid(NewPts, weights, pts, eps);
    else if (WhichAlgm == "aggr")
      PreprocessPtsAggr(NewPts, weights, pts, eps);
    else if (WhichAlgm == "subdiv")
      PreprocessPtsSubdiv(NewPts, weights, pts, eps);
    else if (WhichAlgm == "auto")
      PreprocessPts(NewPts, weights, pts, eps);
    else CoCoA_ERROR(ERR::SERIOUS, "Unknown preprocessing algm");
  }


  RingElem HilbertNumQuot_forC5(const ideal& I)
  {
    PPOrdering PPO(ordering(PPM(RingOf(I))));
    if (GradingDim(PPO) ==0)  // aff -- see hp.cpkg5
      CoCoA_ERROR(ERR::NYI, "HilbertNumQuot_forC5");
    //    if (GradingDim(PPO) >1)  // multigraded
    if (GradingDim(PPO) ==1)  // standard case
    {
      if (IsStdGraded(PPO))
        return HilbertNumQuot(I);
      else
        return MGHilbertNumQuot(I);
    }
    return MGHilbertNumQuot(I);
  }


  std::vector<BigInt> ContFrac_forC5(const BigRat& q)
  {
    vector<BigInt> ans;
    if (IsZero(q)) { ans.push_back(BigInt(0)); return ans; }
    for (ContFracIter iter(q); !IsEnded(iter); ++iter)
      ans.push_back(*iter);
    return ans;
  }


  std::vector<BigRat> CFApproximants_forC5(const BigRat& q)
  {
    vector<BigRat> ans;
    if (IsZero(q)) { ans.push_back(BigRat(0)); return ans; }
    for (CFApproximantsIter iter(q); !IsEnded(iter); ++iter)
      ans.push_back(*iter);
    return ans;
  }

  std::vector<BigInt> BinomialRepr_forC5(const BigInt& N, const BigInt& r)
  {
    long SmallR;
    if (!IsConvertible(SmallR, r))
      CoCoA_ERROR(ERR::ArgTooBig, "BinomialRepr");
    std::vector<BigInt> ans = BinomialRepr(N, SmallR);
    ans.erase(ans.begin()); // erase the first elem (which is always 0)
    return ans;
  }
  
  BigInt BinomialReprShift_forC5(const BigInt& N, const BigInt& r, const BigInt& shift1, const BigInt& shift2)
  {
    long SmallR, SmallShift1, SmallShift2;
    if (!IsConvertible(SmallR, r) ||
        !IsConvertible(SmallShift1, shift1) ||
        !IsConvertible(SmallShift2, shift2))
      CoCoA_ERROR(ERR::ArgTooBig, "BinomialRepr");
    return BinomialReprShift(N, SmallR, SmallShift1, SmallShift2);
  }
  

  BigInt NumPartitions_forC5(const BigInt& N)
  {
    long n;
    if (!IsConvertible(n, N))
      CoCoA_ERROR(ERR::ArgTooBig, "NumPartitions");
    return NumPartitions(n);
  }
  
  
  ring RingQQt_forC5(const BigInt& NumIndets)
  {
    //    long n;
    //    if (!IsConvertible(n, NumIndets)) CoCoA_ERROR(ERR::ArgTooBig, "RingQQt");
    //    return RingQQt(n);
    const ErrorInfo ErrMesg(ERR::ArgTooBig, "RingQQt");
    return RingQQt(ConvertTo<long>(NumIndets, ErrMesg));
  }


// std::vector<symbol> NewSymbols(const BigInt& n)
// {
// return NewSymbols(ConvertTo<long>(n,ErrorInfo("Arg too big","NewSymbols")));
// }
  

  //----- matrix --------------------------------------------------

  matrix NewDenseMat_forC5(const ring& R, const BigInt& NR, const BigInt& NC)
  {
    long r,c;
    if (!IsConvertible(r,NR)) CoCoA_ERROR(ERR::BadRowIndex,"NewDenseMat_forC5");
    if (!IsConvertible(c,NC)) CoCoA_ERROR(ERR::BadColIndex,"NewDenseMat_forC5");
    return NewDenseMat(R, r, c);
  }

// { // AMB 2014-01: same as (but more readable?)
//   return NewDenseMat(R,
//         ConvertTo<long>(NR,ErrorInfo(ERR::BadRowIndex,"NewDenseMat_forC5")),
//         ConvertTo<long>(NC,ErrorInfo(ERR::BadColIndex,"NewDenseMat_forC5")));
// }
  

  matrix HomogElimMat_forC5(ConstMatrixView M, const std::vector<BigInt>& ElimInds)
  {
    std::vector<long> v = VectorLong(ElimInds, "HomogElimMat_forC5");
    for (long i=0; i<len(v); ++i)  --v[i];
    return HomogElimMat(M, v);
  }
  

  matrix ElimMat_forC5(ConstMatrixView M, const std::vector<BigInt>& ElimInds)
  {
    std::vector<long> v = VectorLong(ElimInds, "ElimMat_forC5");
    for (long i=0; i<len(v); ++i)  --v[i];
    return ElimMat(M, v);
  }
  

  matrix ElimMat_forC5(const BigInt& n, const std::vector<BigInt>& ElimInds)
  {
    std::vector<long> v = VectorLong(ElimInds, "ElimMat_forC5");
    for (long i=0; i<len(v); ++i)  --v[i];
    long tmp;
    if (!IsConvertible(tmp, n)) CoCoA_ERROR(ERR::BadMatrixSize,"ElimMat_forC5");
    return ElimMat(tmp, v);
  }
  

  matrix RevLexMat_forC5(const BigInt& N)
  {
    long n;
    if (!IsConvertible(n, N)) CoCoA_ERROR(ERR::BadRowIndex, "RevLexMat_forC5");
    return NewDenseMatRevLex(n);
  }
  

  matrix StdDegLexMat_forC5(const BigInt& N)
  {
    long n;
    if (!IsConvertible(n, N)) CoCoA_ERROR(ERR::BadRowIndex, "StdDegLexMat_forC5");
    return NewDenseMatStdDegLex(n);
  }
  

  matrix StdDegRevLexMat_forC5(const BigInt& N)
  {
    long n;
    if (!IsConvertible(n, N)) CoCoA_ERROR(ERR::BadRowIndex, "StdDegRevLexMat_forC5");
    return NewDenseMatStdDegRevLex(n);
  }
  

  matrix XelMat_forC5(const BigInt& N)
  {
    long n;
    if (!IsConvertible(n, N)) CoCoA_ERROR(ERR::BadRowIndex, "XelMat_forC5");
    return NewDenseMatXel(n);
  }
  

  matrix LexMat_forC5(const BigInt& N)
  {
    long n;
    if (!IsConvertible(n, N)) CoCoA_ERROR(ERR::BadRowIndex, "LexMat_forC5");
    return NewDenseMat(IdentityMat(RingQQ(), n));
  }
  

  matrix ExtractOrdMat_forC5(const matrix& M)
  {
    return NewMatMinimize(M);
  }
  

  void SetEntry_forC5(MatrixView& M, const BigInt& I, const BigInt& J, ConstRefRingElem x)
  {
    long i,j;
    if (!IsConvertible(i, I)) CoCoA_ERROR(ERR::BadRowIndex, "SetEntry_forC5");
    if (!IsConvertible(j, J)) CoCoA_ERROR(ERR::BadColIndex, "SetEntry_forC5");
    SetEntry(M, i-1, j-1, x);
  }
  

  std::vector<RingElem> GetRow_forC5(ConstMatrixView M, const BigInt& I)
  {
    long i;
    if (!IsConvertible(i, I)) CoCoA_ERROR(ERR::BadRowIndex, "GetRow_forC5");
    return GetRow_long(M, i-1);
  }
  

  std::vector<RingElem> GetCol_forC5(ConstMatrixView M, const BigInt& I)
  {
    long i;
    if (!IsConvertible(i, I)) CoCoA_ERROR(ERR::BadColIndex, "GetCol_forC5");
    vector<RingElem> v;
    --i;
    for (long j=0; j<NumRows(M); ++j)  v.push_back(M(j,i));
    return v;
  }
  

  void SetRow_forC5(MatrixView& M, const BigInt& I, const std::vector<RingElem>& v)
  {
    long i;
    if (!IsConvertible(i, I)) CoCoA_ERROR(ERR::BadRowIndex, "SetRow_forC5");
    --i;
    for (long j=0; j<NumCols(M); ++j)  SetEntry(M, i, j, v[j]);
  }
  

  void SetCol_forC5(MatrixView& M, const BigInt& I, const std::vector<RingElem>& v)
  {
    long i;
    if (!IsConvertible(i, I)) CoCoA_ERROR(ERR::BadColIndex, "SetCol_forC5");
    --i;
    for (long j=0; j<NumRows(M); ++j)  SetEntry(M, j, i, v[j]);
  }
  

  void SwapRows_forC5(matrix& M, BigInt i1, BigInt i2)
  {
    long j1, j2;
    if (!IsConvertible(j1, i1)) CoCoA_ERROR(ERR::BadRowIndex, "SwapRows_forC5");
    if (!IsConvertible(j2, i2)) CoCoA_ERROR(ERR::BadRowIndex, "SwapRows_forC5");
    SwapRows(M, --j1, --j2);
  }
  

  void SwapCols_forC5(matrix& M, BigInt i1, BigInt i2)
  {
    long j1, j2;
    if (!IsConvertible(j1, i1)) CoCoA_ERROR(ERR::BadColIndex, "SwapCols_forC5");
    if (!IsConvertible(j2, i2)) CoCoA_ERROR(ERR::BadColIndex, "SwapCols_forC5");
    SwapCols(M, --j1, --j2);
  }


  MatrixView IdentityMat_forC5(const ring& R, const BigInt& N)
  {
    long n;
    if (!IsConvertible(n, N)) CoCoA_ERROR(ERR::BadMatrixSize, "IdentityMat_forC5");
    return IdentityMat(R, n);
  }
  

  MatrixView ZeroMat_forC5(const ring& R, const BigInt& NRows, const BigInt& NCols)
  {
    long nrows, ncols;
    if (!IsConvertible(nrows, NRows)) CoCoA_ERROR(ERR::BadMatrixSize, "ZeroMat_forC5");
    if (!IsConvertible(ncols, NCols)) CoCoA_ERROR(ERR::BadMatrixSize, "ZeroMat_forC5");
    return ZeroMat(R, nrows, ncols);
  }
  

  ConstMatrixView transposed_forC5(ConstMatrixView M)
  { return transpose(M); }


//   std::vector<RingElem> minors_forC5(const BigInt& N, ConstMatrixView M)
//   {
//     vector<RingElem> v;
//     size_t n;
//     if (!IsConvertible(n, N))
//       CoCoA_ERROR(ERR::BadRowIndex, "minors_forC5");
//     if (n>NumRows(M))  CoCoA_ERROR(ERR::BadRowIndex, "minors_forC5");
//     if (n>NumCols(M))  CoCoA_ERROR(ERR::BadColIndex, "minors_forC5");
//     vector<size_t> rows(n);
//     vector<size_t> cols(n);
// ---- we need a tuple-generator ----
//         v.push_back(det(submat(M,rows,cols)));
//     return v;
//   }
//   }




  //-------- ideal --------------------------------------------------



  //-------- modules --------------------------------------------------

  module NewFreeModule_forC5(const ring& R, ConstMatrixView M)
  {
    std::vector<degree> shifts;
    long n = NumCols(M);
    degree d(n);
    for (long i=0; i<NumRows(M); ++i) 
    {
      for (long j=0; j<n; ++j)
        SetComponent(d, j, ConvertTo<BigInt>(M(i,j)));
      shifts.push_back(d);
    }
    return NewFreeModule(R, shifts);
  }
 
  long NumCompts_forC5(const module& M)
  {
    if (!IsFGModule(M))
      CoCoA_ERROR("Expected FGModule", "NumCompts");
    return NumCompts(M);
  }
 
 ModuleElem NewFreeModuleElem(const module& M, const vector<RingElem>& v)
  {
    if (!IsFreeModule(M))
      CoCoA_ERROR("Expected FreeModule", "NewFreeModuleElem");
    const vector<ModuleElem>& e = gens(M);
    if (len(e) != len(v))
      CoCoA_ERROR("incompatible length", "NewFreeModuleElem");
    ModuleElem res(M);
    for (long i=0; i<len(e); ++i)  res += v[i]*e[i];
    return res;
  }
  
  FGModule SubmoduleCols_forC5(const module& F, ConstMatrixView M)
  {
    if (!IsFreeModule(F))
      CoCoA_ERROR("Expected FreeModule", "SubmoduleCols");
    return SubmoduleCols(F, M);
  }
  
  FGModule SubmoduleRows_forC5(const module& F, ConstMatrixView M)
  {
    if (!IsFreeModule(F))
      CoCoA_ERROR("Expected FreeModule", "SubmoduleRows");
    return SubmoduleRows(F, M);
  }
  


  //-------- ExternalLibs ----------------------------------------------

  //-------- FROBBY --------------------------------------------------
  
#ifdef CoCoA_WITH_FROBBY
  std::vector<ideal> FrbPrimaryDecomposition_forC5(const ideal& I)
  {
    std::vector<ideal> v;
    FrbPrimaryDecomposition(v,I);
    return v;
  }

  std::vector<ideal> FrbIrreducibleDecomposition_forC5(const ideal& I)
  {
    std::vector<ideal> v;
    FrbIrreducibleDecomposition(v,I);
    return v;
  }

  std::vector<ideal> FrbAssociatedPrimes_forC5(const ideal& I)
  {
    std::vector<ideal> v;
    FrbAssociatedPrimes(v,I);
    return v;
  }

  ideal FrbAlexanderDual_forC5(const ideal& I, ConstRefRingElem t)
  {
    const ring& R = owner(t);
    if (!IsPolyRing(R)) CoCoA_ERROR(ERR::NotElemPolyRing, "FrbAlexanderDual_forC5");
    if (!IsMonomial(t)) CoCoA_ERROR(ERR::BadArg, "FrbAlexanderDual_forC5");
    if (!IsOne(LC(t))) CoCoA_ERROR(ERR::BadArg, "FrbAlexanderDual_forC5");
    return FrbAlexanderDual(I, LPP(t));
  }

  
#endif // CoCoA_WITH_FROBBY

  //-------- NORMALIZ --------------------------------------------------
  
#ifdef CoCoA_WITH_NORMALIZ
  matrix NmzHilbertBasis_forC5(ConstMatrixView M)
  {
    return NewDenseMat(RingOf(M), Normaliz::HilbertBasis(Normaliz::cone(M, libnormaliz::Type::normalization)));
  }

  std::vector<RingElem> NmzNormalToricRing_forC5(const std::vector<RingElem>& v)
  {
    //convert it to a PPVector, run Normaliz and convert back to vector<RingElem>
    std::vector<RingElem> res;
    convert(res,owner(v[0]),Normaliz::NormalToricRing(Normaliz::MonomialsToPPV(v)));
    return res;
  }

  std::vector<RingElem> NmzIntClosureToricRing_forC5(const std::vector<RingElem>& v)
  {
    //convert it to a PPVector
    PPVector ppv = Normaliz::MonomialsToPPV(v);
    //run Normaliz and convert back to vector<RingElem>
    std::vector<RingElem> res;
    convert(res,owner(v[0]),Normaliz::IntClosureToricRing(ppv));
    return res;
  }


  std::vector<RingElem> NmzIntClosureMonIdeal_forC5(const std::vector<RingElem>& v)
  {
    //convert it to a PPVector
    PPVector ppv = Normaliz::MonomialsToPPV(v);
    //run Normaliz and convert back to vector<RingElem>
    std::vector<RingElem> res;
    convert(res,owner(v[0]),Normaliz::IntClosureMonIdeal(ppv));
    return res;
  }

#endif // CoCoA_WITH_NORMALIZ

}

