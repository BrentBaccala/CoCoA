//   Copyright (c)  2010  John Abbott

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


#include "CoCoA/BuildInfo.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/TmpJBEnv.H"
//#include "CoCoA/assert.H"  // for using TEST_ASSERT
#include "CoCoA/error.H"

using namespace CoCoA;
using namespace std;

// Handy macro for making assertions.
#define TEST_ASSERT(cond)   do { if (!(cond)) CoCoA::AssertionFailed(#cond, __FILE__, __LINE__); } while(0)

void test_one(const std::vector<RingElem>& res)
{
  TEST_ASSERT(res.size() == 1);
  TEST_ASSERT(IsOne(res[0]));
}

void test_principal_ideal(const std::vector<RingElem>& to_test, RingElem result )
{  
  TEST_ASSERT(to_test.size() == 1);
  TEST_ASSERT(to_test[0] == result);
}

void program()
{
  GlobalManager CoCoAFoundations;
  ring Q = RingQQ();


  SparsePolyRing polyRing = NewPolyRing(Q, 1);
  const vector<RingElem> x = indets(polyRing);
  vector<RingElem> pols;
  vector<RingElem> result;

  // test empty input
  try 
  {
    JanetBasis(pols); TEST_ASSERT(!"NEVER GET HERE!");
  }
  catch(const CoCoA::ErrorInfo& err)
  {
    TEST_ASSERT(err == ERR::nonstandard);
  }

  // test zero input
  
  pols.push_back(zero(polyRing));

  try 
  {
    JanetBasis(pols); TEST_ASSERT(!"NEVER GET HERE!");
  }
  catch(const CoCoA::ErrorInfo& err)
  {
    TEST_ASSERT(err == ERR::nonstandard);
  }


  // test input only one
  pols.clear();
  pols.push_back(one(polyRing));

  result = JanetBasis(pols, TQDegree);
  test_one(result);

  result = JanetBasis(pols, TQBlockLow);
  test_one(result);

  result = JanetBasis(pols, TQBlockHigh);
  test_one(result);


  // test input only x
  pols.clear();
  pols.push_back(x[0]);

  result = JanetBasis(pols, TQDegree);
  test_principal_ideal(result, x[0]);

  result = JanetBasis(pols, TQBlockLow);
  test_principal_ideal(result, x[0]);

  result = JanetBasis(pols, TQBlockHigh);
  test_principal_ideal(result, x[0]);


  // test input only x, x^2, x^2 -> result must be x
  pols.clear();
  pols.push_back(x[0]);
  pols.push_back(x[0] * x[0]);
  pols.push_back(x[0] * x[0]);

  result = JanetBasis(pols, TQDegree);
  test_principal_ideal(result, x[0]);

  result = JanetBasis(pols, TQBlockLow);
  test_principal_ideal(result, x[0]);

  result = JanetBasis(pols, TQBlockHigh);
  test_principal_ideal(result, x[0]);




  pols.clear();
  pols.push_back(x[0]);
  JBMill mill = ExtendedJanetBasis(pols);
  
  // Test -- JBNonMultVar;
  std::map<PPMonoidElem, std::vector<bool> > ResJBNonMultVar = JBNonMultVar(mill);
  TEST_ASSERT(1 == ResJBNonMultVar.size());
  std::vector<bool> NonMultVarVector = (*ResJBNonMultVar.begin()).second;
  TEST_ASSERT(1 == NonMultVarVector.size());
  TEST_ASSERT(NonMultVarVector[0] == false);


  // Test -- JBMultVar
  std::map<PPMonoidElem, std::vector<bool> > ResJBMultVar = JBMultVar(mill);
  TEST_ASSERT(1 == ResJBMultVar.size());
  std::vector<bool> MultVarVector = (*ResJBMultVar.begin()).second;
  TEST_ASSERT(1 == MultVarVector.size());
  TEST_ASSERT(MultVarVector[0] == true);


  //Test -- JBIsPommaretBasis
  TEST_ASSERT(JBIsPommaretBasis(mill) == true);

  //Test -- JBIsHomogenous
  TEST_ASSERT(JBIsHomogenous(mill) == true);

  //Test -- JBIsMonomialIdeal
  TEST_ASSERT(JBIsMonomialIdeal(mill) == true);

  //Test -- JBHilbertPol
  SparsePolyRing HilbPolPolyRing = NewPolyRing(Q, SymbolRange("s",0,0), NewStdDegLexOrdering(1));
  const std::vector<RingElem> s = indets(HilbPolPolyRing);

  TEST_ASSERT(IsZero(JBHilbertPol(mill, s[0])));

  //Test -- JBHilbertFunc
  TEST_ASSERT(JBHilbertFunc(mill, BigInt(0)) == BigInt(1));
  TEST_ASSERT(JBHilbertFunc(mill, BigInt(1)) == BigInt(0));
  TEST_ASSERT(JBHilbertFunc(mill, BigInt(2)) == BigInt(0));
  TEST_ASSERT(JBHilbertFunc(mill, BigInt(3)) == BigInt(0));

  //Test -- JBHilbertSeries
  FractionField QR = NewFractionField(polyRing);
  RingHom embeddedHom(EmbeddingHom(QR));
  RingElem qx = embeddedHom(x[0]);
  
  TEST_ASSERT(JBHilbertSeries(mill, qx) == qx/(1-qx));

  //Test -- JBStandardRepresentation
  RingElem elem(x[0]*x[0] + x[0] + 1);
  std::pair<std::map<PPMonoidElem, RingElem>, RingElem> res(JBStandardRepresentation(mill, elem));
  TEST_ASSERT(res.second == one(polyRing));
  TEST_ASSERT((*res.first.begin()).second == x[0] + 1);

  RingElem JBSRZero(zero(polyRing));
  res = JBStandardRepresentation(mill, JBSRZero);
  TEST_ASSERT(res.second == zero(polyRing));
  TEST_ASSERT((*res.first.begin()).second == zero(polyRing));

  //Test -- JBNormalForm
  // RingElem elem(x[0]*x[0] + x[0] + 1); --> defined above
  TEST_ASSERT(one(polyRing) == JBNormalForm(mill, elem));

  // RingElem JBSRZero(zero(polyRing)); --> defined above
  TEST_ASSERT(zero(polyRing) == JBNormalForm(mill, JBSRZero));

  //Test -- JBSyzygy
  FGModule firstSyz = JBSyzygy(mill);
  vector<ModuleElem> gensFirstSyz = gens(firstSyz);
  TEST_ASSERT(0 == gensFirstSyz.size());

  //Test -- JBDim
  TEST_ASSERT(0 == JBDim(mill));

  //Test -- JBProjDim
  TEST_ASSERT(1 - 1 == JBProjDim(mill));

  //Test -- JBDepth
  TEST_ASSERT(1 == JBDepth(mill));

  //Test -- JBCls ConstRefRingElem
  elem = x[0] * x[0] + 1;
  TEST_ASSERT(1 == JBCls(mill, elem));

  //Test -- JBCls PPMonoidElem
  PPMonoidElem ElemPPM(LPP(elem)); // == x^2
  TEST_ASSERT(1 == JBCls(mill, ElemPPM));

  //Test -- JBCls for 1
  elem = one(polyRing);
  TEST_ASSERT(1 == JBCls(mill, elem));

  //Test -- JBMinCls 
  TEST_ASSERT(1 == JBMinCls(mill));

  //Test -- JBElementsWithClass
  std::vector<RingElem> ClassZero(JBElementsWithClass(mill, BigInt(0)));
  std::vector<RingElem> ClassOne(JBElementsWithClass(mill, BigInt(1)));
  std::vector<RingElem> ClassTwo(JBElementsWithClass(mill, BigInt(2)));
  TEST_ASSERT(ClassZero.empty());

  TEST_ASSERT(ClassOne.size() == 1);
  TEST_ASSERT(ClassOne[0] == x[0]);

  TEST_ASSERT(ClassTwo.empty());

  //Test -- JBExtremalBettiNumbers
  std::map<std::pair<long, long>, long> bettis(JBExtremalBettiNumbers(mill));
  /* Expected Result:
      Deg
       1   1  --- x
       |   |      |
       0   0  --- 0
       |   |      |
    Module P0 --- I
    
   */

  TEST_ASSERT(bettis.size() == 1);
  TEST_ASSERT((*bettis.begin()).first.first == 0);
  TEST_ASSERT((*bettis.begin()).first.second == 1);
  TEST_ASSERT((*bettis.begin()).second == 1);

  //Test -- JBRegularSequence
  vector<RingElem> RegSeq(JBRegularSequence(mill));
  TEST_ASSERT(RegSeq.size() == 1);
  TEST_ASSERT(RegSeq[0] == x[0]);

  //Test -- JBMaxStronglyIndependentSet
  vector<RingElem> MaxIndepent(JBMaxStronglyIndependentSet(mill));
  TEST_ASSERT(MaxIndepent.size() == 0);

  //Test -- JBIsCohenMacaulay
  TEST_ASSERT(false == JBIsCohenMacaulay(mill));

  //Test -- JBRegularity && JBCastelnuovoMumfordRegularity
  TEST_ASSERT(1 == JBRegularity(mill));
  TEST_ASSERT(1 == JBCastelnuovoMumfordRegularity(mill));

  //Test -- JBDegPommaretClass
  TEST_ASSERT(-1 == JBDegPommaretClass(mill, BigInt(0)));
  TEST_ASSERT(-1 == JBDegPommaretClass(mill, BigInt(-1)));
  TEST_ASSERT(1 == JBDegPommaretClass(mill, BigInt(1)));
  TEST_ASSERT(-1 == JBDegPommaretClass(mill, BigInt(2)));

  //Test -- JBSatiety
  TEST_ASSERT(1 == JBSatiety(mill));

  //Test -- JBSaturation
  vector<RingElem> saturation(JBSaturation(mill));
  TEST_ASSERT(1 == saturation.size());
  TEST_ASSERT(IsOne(saturation[0]));

  //Test -- JBComplementaryDecomposition
  vector< std::pair<PPMonoidElem, std::vector<bool> > > CompDecomp(JBComplementaryDecomposition(mill));
  TEST_ASSERT(1 == CompDecomp.size());
  TEST_ASSERT(IsOne(CompDecomp[0].first));
  TEST_ASSERT(1 == CompDecomp[0].second.size());
  TEST_ASSERT(false == CompDecomp[0].second[0]);

  //Test -- JBStandardPairs
  vector< std::pair<PPMonoidElem, std::vector<bool> > > StdPairs(JBStandardPairs(mill));
  TEST_ASSERT(1 == StdPairs.size());
  TEST_ASSERT(IsOne(StdPairs[0].first));
  TEST_ASSERT(1 == StdPairs[0].second.size());
  TEST_ASSERT(false == StdPairs[0].second[0]);
  
// std::vector<RingElem> JBSocle(const JBMill& mill)

  pols.clear();
  pols.push_back(one(polyRing));
  JBMill OneMill = ExtendedJanetBasis(pols);

  // Test -- Correct Janet Basis
  std::vector<RingElem> ReturnedJB(JBReturnJB(OneMill));
  TEST_ASSERT(1 == ReturnedJB.size());
  TEST_ASSERT(IsOne(ReturnedJB[0]));

  ResJBNonMultVar = JBNonMultVar(OneMill);
  TEST_ASSERT(1 == ResJBNonMultVar.size());
  NonMultVarVector = (*ResJBNonMultVar.begin()).second;
  TEST_ASSERT(1 == NonMultVarVector.size());
  TEST_ASSERT(NonMultVarVector[0] == false);


  // Test -- JBMultVar
  ResJBMultVar = JBMultVar(OneMill);
  TEST_ASSERT(1 == ResJBMultVar.size());
  MultVarVector = (*ResJBMultVar.begin()).second;
  TEST_ASSERT(1 == MultVarVector.size());
  TEST_ASSERT(MultVarVector[0] == true);


  //Test -- JBIsPommaretBasis
  TEST_ASSERT(JBIsPommaretBasis(OneMill) == true);

  //Test -- JBIsHomogenous
  TEST_ASSERT(JBIsHomogenous(OneMill) == true);

  //Test -- JBIsMonomialIdeal
  TEST_ASSERT(JBIsMonomialIdeal(OneMill) == true);

  //Test -- JBHilbertPol
  TEST_ASSERT(IsZero(JBHilbertPol(OneMill, s[0])));

  //Test -- JBHilbertFunc
  JBHilbertFunc(OneMill, BigInt(0));
  TEST_ASSERT(0 == JBHilbertFunc(OneMill, BigInt(0)));
  TEST_ASSERT(0 == JBHilbertFunc(OneMill, BigInt(1)));
  TEST_ASSERT(0 == JBHilbertFunc(OneMill, BigInt(2)));
  TEST_ASSERT(0 == JBHilbertFunc(OneMill, BigInt(3)));
  
  //Test -- JBHilbertSeries
  JBHilbertSeries(OneMill, qx);
  TEST_ASSERT(JBHilbertSeries(OneMill, qx) == 1/(1-qx));

  //Test -- JBStandardRepresentation
  elem = x[0]*x[0] + x[0] + 1;
  res = JBStandardRepresentation(OneMill, elem);
  TEST_ASSERT(IsZero(res.second));
  TEST_ASSERT((*res.first.begin()).second == elem);

  //Test -- JBNormalForm
  TEST_ASSERT(IsZero(JBNormalForm(OneMill, elem)));

  //Test -- JBSyzygy
  firstSyz = JBSyzygy(OneMill);
  gensFirstSyz = gens(firstSyz);
  TEST_ASSERT(0 == gensFirstSyz.size());

  //Test -- JBDim
  TEST_ASSERT(0 == JBDim(OneMill));

  //Test -- JBProjDim
  TEST_ASSERT(0 == JBProjDim(OneMill));

  //Test -- JBDepth
  TEST_ASSERT(1 == JBDepth(OneMill));

  //Test -- JBMinCls 
  TEST_ASSERT(1 == JBMinCls(OneMill));

  //Test -- JBElementsWithClass
  TEST_ASSERT(JBElementsWithClass(OneMill, BigInt(0)).empty());
  TEST_ASSERT(1 == JBElementsWithClass(OneMill, BigInt(1)).size());
  TEST_ASSERT(JBElementsWithClass(OneMill, BigInt(2)).empty());

  //Test -- JBDegPommaretClass
  TEST_ASSERT(-1 == JBDegPommaretClass(OneMill, BigInt(-1)));
  TEST_ASSERT(-1 == JBDegPommaretClass(OneMill, BigInt(0)));
  TEST_ASSERT(0 == JBDegPommaretClass(OneMill, BigInt(1)));
  TEST_ASSERT(-1 == JBDegPommaretClass(OneMill, BigInt(2)));

  //Test -- JBExtremalBettiNumbers
  std::map<std::pair<long, long>, long> OneBettis(JBExtremalBettiNumbers(OneMill));
  /* Expected Result:
      Deg
       1   0  --- 0
       |   |      |
       0   1  --- 1
       |   |      |
    Module P0 --- I
    
   */

  TEST_ASSERT(OneBettis.size() == 1);
  TEST_ASSERT((*OneBettis.begin()).first.first == 0);
  TEST_ASSERT((*OneBettis.begin()).first.second == 0);
  TEST_ASSERT((*OneBettis.begin()).second == 1);
  
  //Test -- JBRegularSequence
  RegSeq = JBRegularSequence(OneMill);
  TEST_ASSERT(RegSeq.size() == 1);
  TEST_ASSERT(RegSeq[0] == x[0]);
  
  //Test -- JBMaxStronglyIndepenentSet
  TEST_ASSERT(0 == JBMaxStronglyIndependentSet(OneMill).size());

  //Test -- JBIsCohenMacaulay
  TEST_ASSERT(false == JBIsCohenMacaulay(OneMill));

  //Test -- JBRegularity && JBCastelnuovoMumfordRegularity
  TEST_ASSERT(0 == JBRegularity(OneMill));
  TEST_ASSERT(0 == JBCastelnuovoMumfordRegularity(OneMill));

  //Test -- JBSatiety
  TEST_ASSERT(0 == JBSatiety(OneMill));

  //Test -- JBSaturation
  saturation = JBSaturation(OneMill);
  TEST_ASSERT(1 == saturation.size());
  TEST_ASSERT(IsOne(saturation[0]));

  //Test -- JBComplementaryDecomposition
  CompDecomp = JBComplementaryDecomposition(OneMill);
  TEST_ASSERT(0 == CompDecomp.size());
  // TEST_ASSERT(IsOne(CompDecomp[0].first));
  // TEST_ASSERT(1 == CompDecomp[0].second.size());
  // TEST_ASSERT(false == CompDecomp[0].second[0]);

  //Test -- JBStandardPairs
  StdPairs = JBStandardPairs(OneMill);
  TEST_ASSERT(0 == StdPairs.size());
  // TEST_ASSERT(IsOne(StdPairs[0].first));
  // TEST_ASSERT(1 == StdPairs[0].second.size());
  // TEST_ASSERT(false == StdPairs[0].second[0]);
  
}


// Use main() to handle any uncaught exceptions and warn the user about them.
int main()
{
  try
  {
    program();
    return 0;
  }
  catch (const CoCoA::ErrorInfo& err)
  {
    cerr << "***ERROR***  UNCAUGHT CoCoA Error";
    ANNOUNCE(cerr, err);
  }
  catch (const std::exception& exc)
  {
    cerr << "***ERROR***  UNCAUGHT std::exception: " << exc.what() << endl;
  }
  catch(...)
  {
    cerr << "***ERROR***  UNCAUGHT UNKNOWN EXCEPTION" << endl;
  }

  BuildInfo::PrintAll(cerr);
  return 1;
}
