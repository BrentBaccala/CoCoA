//   Copyright (c)  2014  John Abbott

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
#include "CoCoA/BigRat.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingTwinFloat.H"
#include "CoCoA/IntOperations.H"
//#include "CoCoA/assert.H"  // for using TEST_ASSERT
#include "CoCoA/error.H"
#include "CoCoA/utils.H"


#include <iostream>
using std::cerr;
using std::endl;
#include <vector>
using std::vector;

using namespace CoCoA;

// Handy macro for making assertions.
#define TEST_ASSERT(cond)   do { if (!(cond)) CoCoA::AssertionFailed(#cond, __FILE__, __LINE__); } while(0)


void TestCmpAbs(ring R)
{
  vector<RingElem> v;
  for (int i=0; i < 10; ++i)
    v.push_back(RingElem(R,i));

  const int n = len(v);
  for (int i=0; i < n; ++i)
    for (int j=0; j < n; ++j)
    {
      const int sgn = cmp(v[i],v[j]);
      TEST_ASSERT(CmpAbs(v[i],v[j]) == sgn);
      TEST_ASSERT(CmpAbs(-v[i],v[j]) == sgn);
      TEST_ASSERT(CmpAbs(v[i],-v[j]) == sgn);
      TEST_ASSERT(CmpAbs(-v[i],-v[j]) == sgn);
    }
}


void TestCmpAbs2(ring R)
{
  vector<RingElem> v;
  for (int i=1; i <= 100; ++i)
    v.push_back(RingElem(R, BigRat(fibonacci(i+1),fibonacci(i))));

  const int n = len(v);
  for (int i=0; i < n; ++i)
    for (int j=0; j < n; ++j)
    {
      const int sgn = cmp(v[i],v[j]);
      TEST_ASSERT(CmpAbs(v[i],v[j]) == sgn);
      TEST_ASSERT(CmpAbs(-v[i],v[j]) == sgn);
      TEST_ASSERT(CmpAbs(v[i],-v[j]) == sgn);
      TEST_ASSERT(CmpAbs(-v[i],-v[j]) == sgn);
    }
}


void program()
{
  GlobalManager CoCoAFoundations;
  ring ZZ = RingZZ();
  ring QQ = RingQQ();

  TestCmpAbs(ZZ);
  TestCmpAbs(QQ);
  for (int prec=20; prec < 1000; prec += 20)
  {
    ring RR = NewRingTwinFloat(prec);
    TestCmpAbs(RR);
  }

  TestCmpAbs2(QQ);
  for (int prec=20; prec <= 500; prec += 20)
  {
    ring RR = NewRingTwinFloat(prec);
    try { TestCmpAbs2(RR); }
    catch (const RingTwinFloat::InsufficientPrecision&) {}
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
