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

  BigInt N1;
  TEST_ASSERT(IsZero(N1));

  BigInt N2(0);
  TEST_ASSERT(N1 == N2);
  TEST_ASSERT(N2 == N1);

  BigInt N3(1);
  TEST_ASSERT(IsOne(N3));

  const int LIMIT = 100;

  // Check some equalities/inequalities
  for (int i = -LIMIT; i <= LIMIT; ++i)
    for (int j = -LIMIT; j <= LIMIT; ++j)
    {
      if (i == j)
      {
        TEST_ASSERT(i == BigInt(j));
        TEST_ASSERT(BigInt(i) == j);
        TEST_ASSERT(BigInt(i) == BigInt(j));
      }
      if (i != j)
      {
        TEST_ASSERT(i != BigInt(j));
        TEST_ASSERT(BigInt(i) != j);
        TEST_ASSERT(BigInt(i) != BigInt(j));
      }
      if (i < j)
      {
        TEST_ASSERT(i < BigInt(j));
        TEST_ASSERT(BigInt(i) < j);
        TEST_ASSERT(BigInt(i) < BigInt(j));
      }
      if (i <= j)
      {
        TEST_ASSERT(i <= BigInt(j));
        TEST_ASSERT(BigInt(i) <= j);
        TEST_ASSERT(BigInt(i) <= BigInt(j));
      }
      if (i > j)
      {
        TEST_ASSERT(i > BigInt(j));
        TEST_ASSERT(BigInt(i) > j);
        TEST_ASSERT(BigInt(i) > BigInt(j));
      }
      if (i >= j)
      {
        TEST_ASSERT(i >= BigInt(j));
        TEST_ASSERT(BigInt(i) >= j);
        TEST_ASSERT(BigInt(i) >= BigInt(j));
      }
    }

  // Check some cmp function calls.
  for (int i = -LIMIT; i <= LIMIT; ++i)
    for (int j = -LIMIT; j <= LIMIT; ++j)
    {
      if (cmp(i,j) > 0)
      {
        TEST_ASSERT(cmp(i, BigInt(j)) > 0);
        TEST_ASSERT(cmp(BigInt(i), j) > 0);
        TEST_ASSERT(cmp(BigInt(i), BigInt(j)) > 0);
      }
      if (cmp(i,j) == 0)
      {
        TEST_ASSERT(cmp(i, BigInt(j)) == 0);
        TEST_ASSERT(cmp(BigInt(i), j) == 0);
        TEST_ASSERT(cmp(BigInt(i), BigInt(j)) == 0);
      }
      if (cmp(i,j) < 0)
      {
        TEST_ASSERT(cmp(i, BigInt(j)) < 0);
        TEST_ASSERT(cmp(BigInt(i), j) < 0);
        TEST_ASSERT(cmp(BigInt(i), BigInt(j)) < 0);
      }
    }


  // Check some basic arithmetic operations.
  for (int i = -LIMIT; i <= LIMIT; ++i)
    for (int j = -LIMIT; j <= LIMIT; ++j)
    {
      TEST_ASSERT(i+j == i+BigInt(j));
      TEST_ASSERT(i+j == BigInt(i)+j);
      TEST_ASSERT(i+j == BigInt(i)+BigInt(j));

      TEST_ASSERT(i-j == i-BigInt(j));
      TEST_ASSERT(i-j == BigInt(i)-j);
      TEST_ASSERT(i-j == BigInt(i)-BigInt(j));

      TEST_ASSERT(i*j == i*BigInt(j));
      TEST_ASSERT(i*j == BigInt(i)*j);
      TEST_ASSERT(i*j == BigInt(i)*BigInt(j));

      // Restrict to non-negative for division as C++ is not well-defined with negative numbers.
      if (i >= 0 && j > 0)
      {
        TEST_ASSERT(i/j == i/BigInt(j));
        TEST_ASSERT(i/j == BigInt(i)/j);
        TEST_ASSERT(i/j == BigInt(i)/BigInt(j));
      }
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
