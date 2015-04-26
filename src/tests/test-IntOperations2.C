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
#include "CoCoA/IntOperations.H"
#include "CoCoA/error.H"


#include <iostream>
using std::cerr;
using std::endl;
#include <limits>
using std::numeric_limits;

using namespace CoCoA;

// Handy macro for making assertions.
#define TEST_ASSERT(cond)   do { if (!(cond)) CoCoA::AssertionFailed(#cond, __FILE__, __LINE__); } while(0)

void program()
{
  // Sundry tests for edge cases for several integer functions.
  GlobalManager CoCoAFoundations;

  const int zero = 0;
  const BigInt ZERO;
  TEST_ASSERT(IsOne(power(zero,zero)));
  TEST_ASSERT(IsOne(power(ZERO,zero)));
  TEST_ASSERT(IsOne(power(zero,ZERO)));
  TEST_ASSERT(IsOne(power(ZERO,ZERO)));

  TEST_ASSERT(SmallPower(zero,zero) == 1);

  const int IntMax = numeric_limits<int>::max();
  const int IntMin = numeric_limits<int>::min();
  const long LongMax = numeric_limits<long>::max();
  const long LongMin = numeric_limits<long>::min();
  const unsigned long ULongMax = numeric_limits<unsigned long>::max();
  const unsigned long ULongMin = numeric_limits<unsigned long>::min();

  TEST_ASSERT(IsPowerOf2(1));
  TEST_ASSERT(IsPowerOf2(2));
  TEST_ASSERT(!IsPowerOf2(3));
  TEST_ASSERT(IsPowerOf2(4));

  TEST_ASSERT(!IsPowerOf2(-1));
  TEST_ASSERT(!IsPowerOf2(-2));
  TEST_ASSERT(!IsPowerOf2(3));
  TEST_ASSERT(!IsPowerOf2(-4));

  TEST_ASSERT(!IsPowerOf2(IntMax));
  TEST_ASSERT(!IsPowerOf2(IntMin));
  TEST_ASSERT(!IsPowerOf2(LongMax));
  TEST_ASSERT(!IsPowerOf2(LongMin));
  TEST_ASSERT(!IsPowerOf2(ULongMax));
  TEST_ASSERT(!IsPowerOf2(ULongMin));

  TEST_ASSERT(!IsPowerOf2(BigInt(0)));
  for (int n=0; n <= 1024; ++n)
  {
    TEST_ASSERT(IsPowerOf2(power(2,n)));
    TEST_ASSERT(!IsPowerOf2(-power(2,n)));
    TEST_ASSERT(n==0 || !IsPowerOf2(1+power(2,n)));
    TEST_ASSERT(n==1 || !IsPowerOf2(-1+power(2,n)));
    TEST_ASSERT(!IsPowerOf2(3*power(2,n)));
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
