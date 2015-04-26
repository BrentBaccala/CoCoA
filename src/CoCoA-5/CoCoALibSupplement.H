#ifndef CoCoALibSupplement_H
#define CoCoALibSupplement_H

//   Copyright (c)  2011  Anna Bigatti

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

#include <vector>
#include "CoCoA/factorization.H"
#include "CoCoA/error.H"

namespace CoCoA
{
// fwd decl -- defined in library.H
	class BigInt;
	class BigRat;
  struct CoeffPP;
	class ConstMatrixView;
	class RingElemAlias;
  typedef const RingElemAlias& ConstRefRingElem;
	class MatrixView;
	class RingElem;
	class RingHom;
	class degree;
	class ideal;
	class matrix;
	class FGModule;
	class module;
	class ModuleElem;
	class ring;
        class SparsePolyRing;

  std::vector<long> VectorLong(const std::vector<BigInt>& BigIntVec, const std::string& FuncName);
  std::vector<long> VectorLongDecr1(const std::vector<BigInt>& BigIntVec, const ERR::ID& ErrID, const std::string& FuncName);
  
  std::vector<BigInt> DegreeToVec(degree d);  // for Value::from

  RingElem LPP_forC5(ConstRefRingElem f); ///< NB result belongs to owner(f)
  RingElem LPP_forC5(const ModuleElem& f);
  RingElem LM_forC5(ConstRefRingElem f);  ///< NB result belongs to owner(f)
  ModuleElem LM_forC5(const ModuleElem& f);
  bool IsTerm_forC5(ConstRefRingElem f);
  RingElem LT_forC5(ConstRefRingElem f); ///< NB result belongs to owner(f)
  ModuleElem LT_forC5(const ModuleElem& f);
  long LPosn_forC5(const ModuleElem& f);
  long FirstNonZeroPosn_forC5(const ModuleElem& f);
  RingElem CoeffOfTerm_forC5(ConstRefRingElem f, ConstRefRingElem t); ///< NB result belongs to CoeffRing(owner(f))
  std::vector<RingElem> coefficients_forC5(ConstRefRingElem f);  ///< NB result belongs to CoeffRing(owner(f))
  std::vector<RingElem> monomials_forC5(ConstRefRingElem f);  ///< NB result belongs to owner(f)
  std::vector<RingElem> support_forC5(ConstRefRingElem f);  ///< NB result belongs to owner(f)

  RingElem DensePoly_forC5(const ring& P, const BigInt& D);

  BigInt NextPrime_forC5(BigInt N);
  factorization<BigInt> SmoothFactor_forC5(const BigInt& N, const BigInt& TrialLimit);
  BigRat CpuTime_forC5();
  //  std::string date_forC5();
  BigInt date_forC5();
  BigInt TimeOfDay_forC5();
  const BigInt random_forC5(const BigInt& lo, const BigInt& hi);
  
  BigInt lcm_forC5(const std::vector<BigInt>& v);
  BigInt gcd_forC5(const std::vector<BigInt>& v);
  RingElem lcm_forC5(const std::vector<RingElem>& v);
  RingElem gcd_forC5(const std::vector<RingElem>& v);  
  RingElem ContentWRT_forC5(ConstRefRingElem f, const std::vector<RingElem>& v);
  std::vector<CoeffPP> CoefficientsWRT_forC5(ConstRefRingElem f, const std::vector<RingElem>& v);
  std::vector<RingElem> CoeffListWRT_forC5(ConstRefRingElem f, ConstRefRingElem x);

  std::vector<RingElem> QuotientBasis_forC5(const ideal& I);
  RingElem HilbertNumQuot_forC5(const ideal& I);
  std::vector<BigInt> ContFrac_forC5(const BigRat& q);
  std::vector<BigRat> CFApproximants_forC5(const BigRat& q);
  std::vector<BigInt> BinomialRepr_forC5(const BigInt& N, const BigInt& r);
  BigInt BinomialReprShift_forC5(const BigInt& N, const BigInt& r, const BigInt& shift1, const BigInt& shift2);
  BigInt NumPartitions_forC5(const BigInt& N);

  ring RingQQt_forC5(const BigInt& NumIndets);

  //----- matrix --------------------------------------------------
  matrix NewDenseMat_forC5(const ring& R, const BigInt& NR, const BigInt& NC);
  matrix HomogElimMat_forC5(ConstMatrixView M, const std::vector<BigInt>& ElimInds);
  matrix ElimMat_forC5(ConstMatrixView M, const std::vector<BigInt>& ElimInds);
  matrix ElimMat_forC5(const BigInt& n, const std::vector<BigInt>& ElimInds);
  matrix LexMat_forC5(const BigInt& N);  
  matrix RevLexMat_forC5(const BigInt& N);
  matrix StdDegLexMat_forC5(const BigInt& N);
  matrix StdDegRevLexMat_forC5(const BigInt& N);
  matrix XelMat_forC5(const BigInt& N);  
  matrix ExtractOrdMat_forC5(const matrix& M);

  void SetEntry_forC5(MatrixView& M, const BigInt& I, const BigInt& J, ConstRefRingElem x);
  std::vector<RingElem> GetRow_forC5(ConstMatrixView M, const BigInt& I);
  std::vector<RingElem> GetCol_forC5(ConstMatrixView M, const BigInt& I);
  void SetRow_forC5(MatrixView& M, const BigInt& I, const std::vector<RingElem>& v);
  void SetCol_forC5(MatrixView& M, const BigInt& I, const std::vector<RingElem>& v);
  void SwapRows_forC5(matrix& M, BigInt i1, BigInt i2);
  void SwapCols_forC5(matrix& M, BigInt i1, BigInt i2);
  MatrixView IdentityMat_forC5(const ring& R, const BigInt& N);
  MatrixView ZeroMat_forC5(const ring& R, const BigInt& NRows, const BigInt& NCols);
  ConstMatrixView transposed_forC5(ConstMatrixView M);

  //---- points ----
  void NBM_forC5(std::vector<RingElem>& QB, std::vector<RingElem>& BB, std::vector<RingElem>& AV, const SparsePolyRing& P, ConstMatrixView OrigPts, ConstMatrixView OrigTolerance);
  //  void SOI_forC5(std::vector<RingElem>& QB, std::vector<RingElem>& BB, std::vector<RingElem>& AV, const SparsePolyRing& P, ConstMatrixView OrigPts, ConstMatrixView OrigTolerance);
  RingElem ClosePassingPoly_forC5(const ring& P, ConstMatrixView OrigPts, ConstMatrixView OrigTolerance);
  void PreprocessPts_forC5(const std::string& WhichAlgm, std::vector< std::vector<RingElem> >& NewPts, std::vector<long>& weights, ConstMatrixView OrigPts, ConstMatrixView epsilon);

  //---- module ----
  module NewFreeModule_forC5(const ring& R, ConstMatrixView M);
  long NumCompts_forC5(const module& M);
  ModuleElem NewFreeModuleElem(const module& M, const std::vector<RingElem>& v);
  FGModule SubmoduleCols_forC5(const module& F, ConstMatrixView M);
  FGModule SubmoduleRows_forC5(const module& F, ConstMatrixView M);

  //---- FROBBY ----  
#ifdef CoCoA_WITH_FROBBY
  std::vector<ideal> FrbAssociatedPrimes_forC5(const ideal& I);
  std::vector<ideal> FrbIrreducibleDecomposition_forC5(const ideal& I);
  std::vector<ideal> FrbPrimaryDecomposition_forC5(const ideal& I);
  ideal FrbAlexanderDual_forC5(const ideal& I, ConstRefRingElem f);
#endif

  //---- Normaliz ----  
#ifdef CoCoA_WITH_NORMALIZ
  matrix NmzHilbertBasis_forC5(ConstMatrixView M);
  std::vector<RingElem> NmzNormalToricRing_forC5(const std::vector<RingElem>& ppv);
  std::vector<RingElem> NmzIntClosureToricRing_forC5(const std::vector<RingElem>& ppv);
  std::vector<RingElem> NmzIntClosureMonIdeal_forC5(const std::vector<RingElem>& ppv);
#endif // CoCoA_WITH_NORMALIZ

}
#endif // CoCoALibSupplement_H

