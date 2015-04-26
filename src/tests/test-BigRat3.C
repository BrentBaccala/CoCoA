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
//#include "CoCoA/assert.H"  // for using TEST_ASSERT
#include "CoCoA/error.H"
#include "CoCoA/BigRat.H"
#include "CoCoA/NumTheory.H"
#include "CoCoA/utils.H"


#include <iostream>
using std::cerr;
using std::endl;
#include <vector>
using std::vector;

using namespace CoCoA;

// Handy macro for making assertions.
#define TEST_ASSERT(cond)   do { if (!(cond)) CoCoA::AssertionFailed(#cond, __FILE__, __LINE__); } while(0)


// This test checks that CmpAbs(q1,q2) gives the same answer as cmp(abs(q1),abs(q2))
void TestCmpAbs()
{
  const int MAX = 40; // magic number -- takes a reasonable time (~2s) on my computer.
  vector<BigRat> q;
  q.push_back(BigRat(0,1));
  for (int n=1; n <= MAX; ++n)
    for (int d=1; d <= MAX; ++d)
      if (gcd(n,d)==1) q.push_back(BigRat(n,d));

  const int N = len(q);
  for (int i=0; i < N; ++i)
    for (int j=0; j < N; ++j)
    {
      const int ans = cmp(q[i],q[j]); // we know q[i] & q[j] are non-neg
      TEST_ASSERT(ans == CmpAbs(q[i],q[j]));
      TEST_ASSERT(ans == CmpAbs(-q[i],q[j]));
      TEST_ASSERT(ans == CmpAbs(q[i],-q[j]));
      TEST_ASSERT(ans == CmpAbs(-q[i],-q[j]));
    }
}


void TestILogBase()
{
  const int MaxBase = 100;
  const int MaxExp = 100;
  const BigRat SlightlyBigger(101,100);
  const BigRat SlightlySmaller(99,100);
  for (int base=2; base <= MaxBase; ++base)
    for (int exp=-MaxExp; exp <= MaxExp; ++exp)
    {
      BigRat b(base,1);
      TEST_ASSERT(ILogBase(power(b,exp),base) == exp);
      TEST_ASSERT(ILogBase(SlightlyBigger*power(b,exp),base) == exp);
      TEST_ASSERT(ILogBase(SlightlySmaller*power(b,exp),base) == exp-1);
    }
}

void program()
{
  GlobalManager CoCoAFoundations;
  TestCmpAbs();
  TestILogBase();
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
