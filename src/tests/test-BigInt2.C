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
//#include "CoCoA/assert.H"  // for using TEST_ASSERT
#include "CoCoA/error.H"

#include <iostream>
using std::cerr;
using std::endl;

using namespace CoCoA;

// Handy macro for making assertions.
#define TEST_ASSERT(cond)   do { if (!(cond)) CoCoA::AssertionFailed(#cond, __FILE__, __LINE__); } while(0)

void program()
{
  GlobalManager CoCoAFoundations(UseGMPAllocator); // for speed specify GMPAllocator (in case the default changes)

  const int Nmax = 200;

  // factorial function
  TEST_ASSERT(factorial(0) == 1);
  for (int i=1; i <= Nmax; ++i)
  {
    TEST_ASSERT(factorial(i) == factorial(BigInt(i)));
    TEST_ASSERT(i*factorial(i-1) == factorial(i));
  }

  // binomial function
  for (int i=-Nmax; i <= Nmax; ++i)
    for (int j=0; j <= Nmax; ++j)
    {
      if (j == 0)
      {
        TEST_ASSERT(binomial(i,j) == 1);
        TEST_ASSERT(binomial(BigInt(i),j) == 1);
        TEST_ASSERT(binomial(i,BigInt(j)) == 1);
        TEST_ASSERT(binomial(BigInt(i),BigInt(j)) == 1);
        continue;
      }
      if (i >= 0 && i < j)
      {
        TEST_ASSERT(binomial(i,j) == 0);
        TEST_ASSERT(binomial(BigInt(i),j) == 0);
        TEST_ASSERT(binomial(i,BigInt(j)) == 0);
        TEST_ASSERT(binomial(BigInt(i),BigInt(j)) == 0);
        continue;
      }

      const BigInt b = binomial(i,j);
      TEST_ASSERT(b == binomial(BigInt(i),j));
      TEST_ASSERT(b == binomial(i,BigInt(j)));
      TEST_ASSERT(b == binomial(BigInt(i),BigInt(j)));
      TEST_ASSERT(b == binomial(i-1,j-1) + binomial(i-1,j));
    }

  // fibonacci
  TEST_ASSERT(fibonacci(0) == 0);
  TEST_ASSERT(fibonacci(1) == 1);
  for (int i=0; i <= Nmax; ++i)
  {
    const BigInt a = fibonacci(i);
    const BigInt b = fibonacci(i+1);
    const BigInt c = fibonacci(i+2);
    TEST_ASSERT(a+b == c);
    TEST_ASSERT(a == fibonacci(BigInt(i)));
  }

  // RoundDiv
  for (int i = -Nmax; i <= Nmax; ++i)
    for (int j = -Nmax; j <= Nmax; ++j)
    {
      if (j == 0) continue;
      const int q = RoundDiv(i,j);
      TEST_ASSERT(q == RoundDiv(BigInt(i), j));
      TEST_ASSERT(q == RoundDiv(i, BigInt(j)));
      TEST_ASSERT(q == RoundDiv(BigInt(i), BigInt(j)));
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
