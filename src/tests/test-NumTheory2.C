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
#include "CoCoA/IntOperations.H"
#include "CoCoA/NumTheory.H"
#include "CoCoA/error.H"
#include "CoCoA/time.H"

#include <iostream>
using std::cerr;
using std::endl;

using namespace CoCoA;

// Handy macro for making assertions.
#define TEST_ASSERT(cond)   do { if (!(cond)) CoCoA::AssertionFailed(#cond, __FILE__, __LINE__); } while(0)

void TestPrime(long p)
{
  TEST_ASSERT(IsPrime(p));

  const long g = PrimitiveRoot(p);
  const long order = p-1;
  for (int i=1; i < p; ++i)
  {
    const long x = PowerMod(g,i,p);
    TEST_ASSERT(MultiplicativeOrderMod(x,p) == order/long(gcd(order,i)));
  }
}

void TestPrime(const BigInt& P)
{
  TEST_ASSERT(IsProbPrime(P));

  const double StartTime = CpuTime();
  const long g = PrimitiveRoot(P);
  const BigInt order(P-1);
  for (int i=1; i <= 9; ++i)
  {
    if (CpuTime() > StartTime+5) break;
    const BigInt x = PowerMod(g,i,P);
    TEST_ASSERT(MultiplicativeOrderMod(x,P) == order/gcd(order,i));
  }
}

void program()
{
  // This test does virtually nothing, but is a handy template if you want
  // to create your own test code: just copy this file and add your code
  // after the line below -- remember that you must not use CoCoALib functions
  // without first creating a GlobalManager.
  GlobalManager CoCoAFoundations;

  int p = 2;
  while (p < 2000)
  {
    TestPrime(p);
    p = NextPrime(p);
  }


  BigInt P(2);
  while (P < 2000)
  {
    TestPrime(P);
    P = NextProbPrime(P);
  }

  TestPrime(NextProbPrime(power(2,32)));
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
