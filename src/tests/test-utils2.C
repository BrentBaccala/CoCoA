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
#include "CoCoA/IntOperations.H"
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


// This is the test:
// let X be the limit value
// confirm X^2 fits
// confirm (X+1)^2 doesn't
// confirm that X^2 is not zero
template<typename T>
void check()
{
  const T MaxVal = numeric_limits<T>::max();
  const T MaxSqrt = MaxSquarableInteger<T>();
  TEST_ASSERT(MaxVal/MaxSqrt >= MaxSqrt);
  TEST_ASSERT(MaxVal/(MaxSqrt+1) < MaxSqrt+1);
  const T MaxSqrt2 = MaxSqrt*MaxSqrt;
  TEST_ASSERT(MaxSqrt2 != 0);
}


void program()
{
  // Test the fn NumericLimits<T>::MaxSquarableInteger()

  GlobalManager CoCoAFoundations;

  check<signed char>();
  check<unsigned char>();

  check<short>();
  check<unsigned short>();

  check<int>();
  check<unsigned int>();

  check<long>();
  check<unsigned long>();

  // long long?
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
