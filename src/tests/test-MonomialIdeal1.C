//   Copyright (c)  2008  Anna Bigatti

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
#include "CoCoA/FractionField.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingDistrMPolyClean.H"
#include "CoCoA/RingDistrMPolyInlFpPP.H"
#include "CoCoA/RingDistrMPolyInlPP.H"
#include "CoCoA/RingFp.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/symbol.H"
// only for testing:
#include "CoCoA/time.H"
#include "CoCoA/VectorOperations.H"

using namespace CoCoA;

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <vector>
using std::vector;

//----------------------------------------------------------------------
// Test for monomial ideal operations
// functions: GBasis, intersect // todo NF, IsElem, IsContained
// environments: DMP (Fp, Q), DMPI (Q), DMPII
//----------------------------------------------------------------------

#define TEST_ASSERT(cond)   do { if (!(cond)) CoCoA::AssertionFailed(#cond, __FILE__, __LINE__); } while(0)


void TestSparsePolyRing(SparsePolyRing P)
{
  //  cout << "TESTING: " << P << endl << endl;

  RingElem x = indet(P,0);
  RingElem y = indet(P,1);
  RingElem z = indet(P,2);

  ideal I = ideal(x*y, y*z);
  ideal J = ideal(x*x, x*y*y);
  //  cout << I << endl;
  //  cout << J << endl;
  J = intersect(I, J);
  //  cout << J << endl;
  //  cout << TidyGens(J) << endl;
  TEST_ASSERT(J == ideal(x*y*y, x*x*y));
  const vector<RingElem>& g = TidyGens(J);
  TEST_ASSERT(g.size() == 2);
}


void program()
{
  GlobalManager CoCoAFoundations;

  const QuotientRing Fp = NewZZmod(101);
  const SparsePolyRing Qxyz = NewPolyRing(RingQQ(), symbols("x", "y", "z"));
  const SparsePolyRing Qx = NewPolyRing_DMPI(RingQQ(), 4);
  const SparsePolyRing Fpy = NewPolyRing_DMPI(Fp, 3);
  const SparsePolyRing FpyII = NewPolyRing_DMPII(Fp, 4);

  TestSparsePolyRing(Qxyz);
  TestSparsePolyRing(Fpy);
  TestSparsePolyRing(FpyII);
  TestSparsePolyRing(Qx);

  {

  ring Fp = NewZZmod(32003);          // coefficient ring
  SparsePolyRing Fpx = NewPolyRing(Fp, 8); // Fp[x[0..7]]
  SparsePolyRing P = Fpx;
  double t0;  // for CpuTime  
  bool IsPrintingMode = false;
  //  IsPrintingMode = true;

  const vector<RingElem>& x = indets(P);
  vector<RingElem> g;
  back_inserter(g) = power(x[2],2) * power(x[5],4);
  back_inserter(g) = power(x[1],3) * power(x[4],4);
  back_inserter(g) = power(x[1],3) * power(x[5],4);
  back_inserter(g) = power(x[3],3) * power(x[6],4);
  back_inserter(g) = power(x[3],4) * power(x[6],3);

  ideal J1(g);
  ideal J2(x[1]*x[2], x[2]*x[3]);
  if (IsPrintingMode) cout << "J1  = " << J1 << endl; 
  if (IsPrintingMode) cout << "J2  = " << J2 << endl << endl;

  TEST_ASSERT(AreGensMonomial(J1));
  TEST_ASSERT(!AreGensSquareFreeMonomial(J1));

  TEST_ASSERT(AreGensMonomial(J2));
  TEST_ASSERT(AreGensSquareFreeMonomial(J2));

  TEST_ASSERT(!AreGensMonomial(ideal(x[0],x[1]-1)));
  TEST_ASSERT(!AreGensSquareFreeMonomial(ideal(x[0],x[1]-1)));

  TEST_ASSERT(AreGensMonomial(ideal(x[0],x[1]*x[1])));
  TEST_ASSERT(!AreGensSquareFreeMonomial(ideal(x[0],x[1]*x[1])));

  t0 = CpuTime();
  ideal I = intersect(J1, J2);
  if (IsPrintingMode) cout << "Cpu Time = " << CpuTime()-t0 << endl;
  if (IsPrintingMode) cout << "intersect(J1, J2) = " << I << endl;
  if (IsPrintingMode) cout << endl;
  TEST_ASSERT(I ==
              ideal(power(x[2],2)*x[3]*power(x[5],4), x[1]*power(x[2],2)*power(x[5],4),
                    x[2]*power(x[3],3)*power(x[6],4), x[2]*power(x[3],4)*power(x[6],3)) +
              ideal(power(x[1],3)*x[2]*power(x[5],4), power(x[1],3)*x[2]*power(x[4],4)));
  

  t0 = CpuTime();
  vector<ideal> PrimDec = PrimaryDecomposition(J2);
  if (IsPrintingMode) cout << "Cpu Time = " << CpuTime()-t0 << endl;
  if (IsPrintingMode) cout << "PrimaryDecomposition(J2) = " << PrimDec << endl;
  TEST_ASSERT(PrimDec[0] == ideal(x[2]));
  TEST_ASSERT(PrimDec[1] == ideal(x[1], x[3]));

  t0 = CpuTime();
  I = J1 * J2;
  if (IsPrintingMode) cout << "Cpu Time = " << CpuTime()-t0 << endl;
  if (IsPrintingMode) cout << "J1 * J2 = " << I << endl;
  if (IsPrintingMode) cout << endl;

  t0 = CpuTime();
  I = colon(J1, J2);
  if (IsPrintingMode) cout << "Cpu Time = " << CpuTime()-t0 << endl;
  if (IsPrintingMode) cout << "colon(J1, J2) = " << I << endl;
  if (IsPrintingMode) cout << endl;
  TEST_ASSERT(I ==
              ideal(x[2]*power(x[5],4), power(x[3],3)*power(x[6],4),
                    power(x[3],4)*power(x[6],3), power(x[1],3)*power(x[5],4))+
              ideal(power(x[1],3)*power(x[4],4),
                    power(x[1],2)*power(x[3],2)*power(x[5],4)*power(x[6],4),
                    power(x[1],2)*power(x[3],2)*power(x[4],4)*power(x[6],4),
                    power(x[1],2)*power(x[3],3)*power(x[5],4)*power(x[6],3))+
              ideal(power(x[1],2)*power(x[3],3)*power(x[4],4)*power(x[6],3)));

  t0 = CpuTime();
  std::vector<RingElem> ElimInds;
  ElimInds.push_back(x[1]);
  ElimInds.push_back(x[2]);  
  I = J1;
  MakeUnique(I)->myElim(ElimInds);
  if (IsPrintingMode) cout << "Cpu Time = " << CpuTime()-t0 << endl;
  if (IsPrintingMode) cout << "elim(ElimInds, J2) = " << I << endl;
  if (IsPrintingMode) cout << endl;
  TEST_ASSERT(I ==
              ideal(power(x[3],3)*power(x[6],4),
                    power(x[3],4)*power(x[6],3)));
  }

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
