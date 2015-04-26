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
#include "CoCoA/utils.H"


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
  // Test the fn LongRange.
  // LongRange(lo,hi) should produce a vector<long> filled with
  // lo, lo+1, ..., hi-1, hi
  // If lo > hi then it produces an empty vector.
  GlobalManager CoCoAFoundations;

  const long MinLong = numeric_limits<long>::min();
  const long MaxLong = numeric_limits<long>::max();

  // Test when lo <= hi
  TEST_ASSERT(len(LongRange(MaxLong,MaxLong)) == 1);
  TEST_ASSERT(len(LongRange(MinLong,MinLong)) == 1);
  TEST_ASSERT(len(LongRange(0,0)) == 1);
  TEST_ASSERT(len(LongRange(1,1)) == 1);
  TEST_ASSERT(len(LongRange(-1,-1)) == 1);
  TEST_ASSERT(len(LongRange(0,1)) == 2);
  TEST_ASSERT(len(LongRange(-1,0)) == 2);
  TEST_ASSERT(len(LongRange(-1,1)) == 3);

  // Test when lo > hi
  TEST_ASSERT(len(LongRange(0,-1)) == 0);
  TEST_ASSERT(len(LongRange(1,0)) == 0);
  TEST_ASSERT(len(LongRange(MinLong+1,MinLong)) == 0);
  TEST_ASSERT(len(LongRange(MaxLong,MaxLong-1)) == 0);
  TEST_ASSERT(len(LongRange(MaxLong,MinLong)) == 0);

  // Test special "impossible" case: LongRange(MinLong,MaxLong)
  try
  {
    std::vector<long> VeryBig = LongRange(MinLong,MaxLong);
    TEST_ASSERT(!"NEVER GET HERE!");
  }
  catch (const ErrorInfo& err)
  {
    TEST_ASSERT(err == ERR::ArgTooBig);
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
