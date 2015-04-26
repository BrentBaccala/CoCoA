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
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/PPMonoid.H"
#include "CoCoA/PPMonoidEv.H"
#include "CoCoA/PPMonoidEvOv.H"
#include "CoCoA/PPMonoidEvZZ.H"
#include "CoCoA/PPMonoidOv.H"
#include "CoCoA/PPMonoidSparse.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/IntOperations.H"
#include "CoCoA/degree.H"
#include "CoCoA/error.H"
#include "CoCoA/matrix.H"
#include "CoCoA/symbol.H"

#include<iostream>
using std::cout;
using std::cerr;
using std::endl;
#include<vector>
using std::vector;

//----------------------------------------------------------------------
// Test for PPMonoidOv functions on PPMonoidElem
// functions: *, wdeg, CmpWDeg, >, IsOne, ==, IsDivisible, /, IsCoprime, colon, 
//            gcd, lcm, power, cmp, log, IndetPower, wdeg
// environments: PPMonoidEv, PPMonoidOv, PPMonoidEvOv, PPMonoidBigEv;
//               lex, DegLex, DegRevLex, MatOrd
//----------------------------------------------------------------------

using namespace CoCoA;

#define TEST_ASSERT(cond)   do { if (!(cond)) CoCoA::AssertionFailed(#cond, __FILE__, __LINE__); } while(0)


// convention: a function containing a "new" should be named "New.."
matrix NewMatrixFromC(ring Z, int cmat[4][4])
{
  matrix M(NewDenseMat(Z,4,4));
  
  for (int i=0; i < 4; ++i)
    for (int j=0; j < 4; ++j)
      SetEntry(M, i, j, cmat[i][j]);
  return M;
}


//-- TestPPMonoid ------------------------------------------------------
// behaviour of different orderings and gradings on PPMonoid and PolyRing

void TestPPMonoid(PPMonoid PPM)
{
  cout << "TESTING: " << PPM << endl << endl;
  cout << std::boolalpha; // prints true/false for bool

  vector<PPMonoidElem> x;
  for (long i=0; i < NumIndets(PPM); ++i)
    x.push_back(indet(PPM, i));

  
  const PPMonoidElem t1 = x[0]*x[314];
  const PPMonoidElem t2 = x[200];
  cout << "Given t1 = " << t1 << endl;
  cout << "  and t2 = " << t2 << endl << endl;

  cout << "t1*t2*t2*t1   gives  " << t1*t2*t2*t1 << endl;
  cout << "t1*t2*t2*t1*x[385]   gives  " << t1*t2*t2*t1*x[385] << endl;
//   cout << "gcd(t1,t2)   gives  " << gcd(t1,t2) << endl;
//   cout << "lcm(t1,t2)   gives  " << lcm(t1,t2) << endl;
//   cout << "wdeg(t1)   gives  " << wdeg(t1) << endl;
//   cout << "wdeg(t2)   gives  " << wdeg(t2) << endl;
//   cout << "CmpWDeg(t1,t2)   gives  " << CmpWDeg(t1,t2) << endl;
//   if ( GradingDim(PPM) > 0 )
//     cout << "CmpWDegPartial(t1,t2,1)   gives  " << CmpWDegPartial(t1,t2,1) << endl;
//   cout << "t1 > t2   gives  " << (t1 > t2) << endl;

  // Check some computations (without printing any results).
  PPMonoidElem one(PPM);
  TEST_ASSERT(IsOne(one));
  TEST_ASSERT(!IsOne(t1));
  TEST_ASSERT(!IsOne(t2));

  TEST_ASSERT(IsIndet(x[2]));
  TEST_ASSERT(IsIndetPosPower(x[2]));
  TEST_ASSERT(!IsIndet(t1));
  TEST_ASSERT(!IsIndetPosPower(t1));
  PPMonoidElem t1t2 = t1*t2;
  PPMonoidElem t2t1 = t2*t1;
  TEST_ASSERT(t1t2 == t2t1);
  TEST_ASSERT(t1 != t2);
  TEST_ASSERT(t1t2 != t1);
  TEST_ASSERT(t1t2 != t2);
  TEST_ASSERT(IsDivisible(t1t2, t1));
  TEST_ASSERT(IsDivisible(t1t2, t2));
  TEST_ASSERT(!IsDivisible(t1, t2));
  TEST_ASSERT(!IsDivisible(t2, t1));
//   TEST_ASSERT(t1t2/t1 == t2);
//   TEST_ASSERT(t1t2/t2 == t1);
  TEST_ASSERT(IsCoprime(t1,t2));
  TEST_ASSERT(!IsCoprime(t1t2,t1));
  TEST_ASSERT(!IsCoprime(t1t2, t2));
//   TEST_ASSERT(colon(t1, t2) == t1);
//   TEST_ASSERT(colon(t2, t1) == t2);
//   TEST_ASSERT(colon(t1t2, t1) == t2);
//   TEST_ASSERT(colon(t1t2, t2) == t1);
//   TEST_ASSERT(gcd(t1*t1t2, t2*t1t2) == t1t2);
//   TEST_ASSERT(lcm(t1*t1t2, t2*t1t2) == power(t1t2,2));

//   TEST_ASSERT(cmp(t1, one) > 0);
//   TEST_ASSERT(cmp(t2, one) > 0);
//   TEST_ASSERT(cmp(one, t1) < 0);
//   TEST_ASSERT(cmp(one, t2) < 0);
//   TEST_ASSERT(cmp(t1t2, t1) > 0);
//   TEST_ASSERT(cmp(t1t2, t2) > 0);
//   TEST_ASSERT(cmp(t1, t1t2) < 0);
//   TEST_ASSERT(cmp(t2, t1t2) < 0);
//   TEST_ASSERT((t1 > t2)^(t1 <= t2));
//   TEST_ASSERT((t2 > t1)^(t2 <= t1));
//   TEST_ASSERT((t1 >= t2)^(t1 < t2));
//   TEST_ASSERT((t2 >= t1)^(t2 < t1));

//   TEST_ASSERT(power(t1t2, 3) > t1t2);
//   TEST_ASSERT(exponent(power(t1,3),0) == 3*exponent(t1,0));
//   TEST_ASSERT(exponent(power(t1,3),1) == 3*exponent(t1,1));
//   TEST_ASSERT(exponent(power(t1,3),2) == 3*exponent(t1,2));
//   TEST_ASSERT(exponent(power(t1,3),3) == 3*exponent(t1,3));
//   TEST_ASSERT(t1 == IndetPower(PPM, 0, exponent(t1,0))*
//               IndetPower(PPM, 314, exponent(t1,314)));
  
  AssignOne(t2t1);
  TEST_ASSERT(IsOne(t2t1));

//   TEST_ASSERT(wdeg(t1) + wdeg(one) == wdeg(t1));
//   TEST_ASSERT(wdeg(t2) + wdeg(one) == wdeg(t2));
//   TEST_ASSERT(wdeg(t1) + wdeg(t2) == wdeg(t1t2));

  cout << "------------------------------------------------" << endl << endl;
}


//-- program --------------------------------------------------------------
// we run TestPolyRing on predefined and user-defined orderings

void program()
{
  GlobalManager CoCoAFoundations;

  // each ordering is degree-compatible with grading over Z^GradingDim
  // i.e. the grading is given by the first GradingDim rows
  // of the ordering matrix
 
  const int n = 400;

  // predefined orderings and gradings
  PPOrdering lex =       NewLexOrdering(n);        // GradingDim = 0
  PPOrdering DegLex =    NewStdDegLexOrdering(n);  // GradingDim = 1
  PPOrdering DegRevLex = NewStdDegRevLexOrdering(n); // GradingDim = 1

  // NB: next version we should have a better way of constructing indet names
  const vector<symbol> X = SymbolRange("x", 0, n-1);

  TestPPMonoid(NewPPMonoidEv(X, lex));

  TestPPMonoid(NewPPMonoidSparse(X, lex));
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
