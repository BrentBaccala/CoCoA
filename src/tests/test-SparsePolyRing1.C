//   Copyright (c)  2007  John Abbott

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
#include "CoCoA/IntOperations.H"
#include "CoCoA/matrix.H"
#include "CoCoA/MatrixForOrdering.H"
#include "CoCoA/MatrixView.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingDistrMPolyClean.H"
#include "CoCoA/RingDistrMPolyInlFpPP.H"
#include "CoCoA/RingDistrMPolyInlPP.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/symbol.H"
#include "CoCoA/degree.H"

using namespace CoCoA;

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <vector>
using std::vector;

//----------------------------------------------------------------------
// Test for SparsePolyRing functions on RingElem
// functions: indets, *, wdeg, StdDeg, IsHomog, homog, gcd, CmpWDeg
// environments: DMP (Fp, ZZ), DMPI (Q), DMPII
//----------------------------------------------------------------------

#define TEST_ASSERT(cond)   do { if (!(cond)) CoCoA::AssertionFailed(#cond, __FILE__, __LINE__); } while(0)


void TestSparsePolyRing(SparsePolyRing P)
{
  cout << "TESTING: " << P << endl << endl;
  cout << std::boolalpha; // prints true/false for bool

  
  const vector<RingElem>& x = indets(P);
  RingElem f1 = x[0]*x[0] + 3*x[1] -2,
    f2 = x[1]*x[1] - 5*x[0],
    f3 = x[0]-x[1];
  ConstRefRingElem h = x[NumIndets(P)-1];

  cout << "Given f1 = " << f1 << endl;
  cout << "  and f2 = " << f2 << endl << endl;
 
  cout << "f1*f2   gives  " << f1*f2 << endl;
  //  cout << "gcd(f1,f2)   gives  " << gcd(f1,f2) << endl;
  cout << "wdeg(f1)   gives  " << wdeg(f1) << endl;
  cout << "wdeg(f2)   gives  " << wdeg(f2) << endl;
  cout << "StdDeg(f1)   gives  " << StdDeg(f1) << endl;
  cout << "StdDeg(f2)   gives  " << StdDeg(f2) << endl;
  cout << "IsHomog(f1)   gives  " << IsHomog(f1) << endl;
  cout << "IsHomog(f2)   gives  " << IsHomog(f2) << endl;
  cout << "IsHomog(f3)   gives  " << IsHomog(f3) << endl;
  cout << "CmpWDeg(f1,f2)        gives  " << sign(CmpWDeg(f1,f2)) << endl; // CmpWDeg guarantees only the sign of result
  if (GradingDim(P)>0)
  {
    cout << "LF(f3)   gives  " << LF(f3) << endl;
    cout << "CmpWDegPartial(f1,f2) gives  "
         << sign(CmpWDegPartial(f1,f2,1)) << endl;  // CmpWDeg guarantees only the sign of result
    cout << "IsHomogPartial(f3,1)  gives  "
                   << IsHomogPartial(f3,1) << endl;
  }
  if (GradingDim(P)<2)
  {
    cout << "  -- homogenizing with h = " << h << endl;
    cout << "homog(f1, h)   gives  " << homog(f1,h) << endl;
    cout << "homog(f2, h)   gives  " << homog(f2,h) << endl;
    cout << "homog(f3, h)   gives  " << homog(f3,h) << endl;
  }

  if ( IsField(CoeffRing(P)) || IsTrueGCDDomain(CoeffRing(P)))
  {
    TEST_ASSERT( IsInvertible(gcd(f1*f2, f1*f3)/f1) );
    TEST_ASSERT( IsInvertible(gcd(f1*f2, f2*f3)/f2) );
    TEST_ASSERT( IsInvertible(gcd(f1*f3, f2*f3)/f3) );
  }
  if ( IsFractionFieldOfGCDDomain(CoeffRing(P)) )
  {
    TEST_ASSERT( CommonDenom(x[0]/3+x[1]/2) ==  6);
    TEST_ASSERT( ClearDenom(x[0]/3+x[1]/2) ==  2*x[0]+3*x[1]);
  }
  
  // Check some computations (without printing any results).
  cout << "one(P)  " << one(P) << endl;
  TEST_ASSERT(IsOne(one(P)));
  TEST_ASSERT(!IsOne(f1));
  TEST_ASSERT(!IsOne(f2));

  RingElem f1f2 = f1*f2;
  RingElem f2f1 = f2*f1;
  TEST_ASSERT(f1f2 == f2f1);
  TEST_ASSERT(f1 != f2);
  TEST_ASSERT(f1f2 != f1);
  TEST_ASSERT(f1f2 != f2);
  TEST_ASSERT(IsDivisible(f1f2, f1));
  TEST_ASSERT(IsDivisible(f1f2, f2));
  TEST_ASSERT(!IsDivisible(f1, f2));
  TEST_ASSERT(!IsDivisible(f2, f1));
  TEST_ASSERT(f1f2/f1 == f2);
  TEST_ASSERT(f1f2/f2 == f1);
  TEST_ASSERT(power(f1,2) == f1*f1);
  TEST_ASSERT(power(f2,2) == f2*f2);
  TEST_ASSERT((power(f1f2,2))/(f1f2*f1f2) == 1);
  TEST_ASSERT(deriv(f1f2, x[0]) == deriv(f1, x[0])*f2 + f1*deriv(f2, x[0]));
  TEST_ASSERT(deriv(x[1]+2*x[0], x[0]) == 2);

  //  TEST_ASSERT(gcd(f1*f1f2, f2*f1f2) == f1f2);

  TEST_ASSERT(CmpWDeg(f1, one(P)) > 0);
  TEST_ASSERT(CmpWDeg(f2, one(P)) > 0);
  TEST_ASSERT(CmpWDeg(one(P), f1) < 0);
  TEST_ASSERT(CmpWDeg(one(P), f2) < 0);
  TEST_ASSERT(CmpWDeg(f1f2, f1) > 0);
  TEST_ASSERT(CmpWDeg(f1f2, f2) > 0);
  TEST_ASSERT(CmpWDeg(f1, f1f2) < 0);
  TEST_ASSERT(CmpWDeg(f2, f1f2) < 0);

  TEST_ASSERT(CmpWDeg(power(f1f2, 3), f1f2) > 0);

  f2f1 = 1;
  TEST_ASSERT(IsOne(f2f1));

  TEST_ASSERT(wdeg(f1) + wdeg(one(P)) == wdeg(f1));
  TEST_ASSERT(wdeg(f2) + wdeg(one(P)) == wdeg(f2));
  TEST_ASSERT(wdeg(f1) + wdeg(f2) == wdeg(f1f2));

  vector<long> exps0(NumIndets(P));  exps0[0] = 1;
  vector<long> exps1(NumIndets(P));  exps1[1] = 1;
  RingElem g0(x[0]), g1(x[1]);
  PushBack(g0, -one(CoeffRing(P)), exps1);
  PushFront(g1, one(CoeffRing(P)), exps0);
  TEST_ASSERT(g0 == x[0]-x[1]);
  TEST_ASSERT(g1 == x[0]+x[1]);

  cout << "------------------------------------------------" << endl << endl;
}


void program()
{
  GlobalManager CoCoAFoundations(UseSymmResidues);
  const vector<RingElem> v(3, one(RingZZ()));
  ConstMatrixView I = IdentityMat(RingZZ(), 3);
  PPOrdering Deg2(NewMatrixOrdering(3, 2, NewMatCompleteOrd(ConcatVer(RowMat(v), I))));

  const QuotientRing Fp = NewZZmod(101);
  const SparsePolyRing ZZxy = NewPolyRing(RingZZ(), 2);
  const SparsePolyRing Fpy = NewPolyRing_DMPI(Fp, 3);
  const SparsePolyRing FpyII = NewPolyRing_DMPII(Fp, 4);
  const SparsePolyRing Qx = NewPolyRing_DMPI(RingQQ(), 4);
  const SparsePolyRing Q2x = NewPolyRing_DMPI(RingQQ(), symbols("x","y","z"), Deg2);

  TestSparsePolyRing(ZZxy);
  TestSparsePolyRing(Fpy);
  TestSparsePolyRing(FpyII);
  TestSparsePolyRing(Qx);
  TestSparsePolyRing(Q2x);
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
