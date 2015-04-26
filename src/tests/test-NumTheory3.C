//   Copyright (c)  2012  John Abbott

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
#include "CoCoA/convert.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/NumTheory.H"
#include "CoCoA/error.H"


#include <iostream>
using std::cerr;
using std::cout;
using std::endl;

using namespace CoCoA;

// Handy macro for making assertions.
#define TEST_ASSERT(cond)   do { if (!(cond)) CoCoA::AssertionFailed(#cond, __FILE__, __LINE__); } while(0)

void program()
{
  // Test for functions related to continued fractions.
  GlobalManager CoCoAFoundations;

  // Generator/iterator for continued fraction "quotients"
  BigRat x(-101,100);
  while (x < 0)
  {
    ContFracIter cfi(x);
    cout << "The continued fraction quotients of " << x << " are:";
    while (!IsEnded(cfi))
      cout << "  " << *cfi++;
    cout << endl;
    x += BigRat(1,10);
  }

  // Make a rational from its cont frac quotients...
  ContFracApproximant approx;
  for (ContFracIter it(x); !IsEnded(it); ++it)
    approx.myAppendQuot(quot(it));
  TEST_ASSERT(approx.myRational() == x);

  // Compute the approximants for a rational...
  x = BigRat(-123,100);
  for (CFApproximantsIter it(x); !IsEnded(it); ++it)
    cout << "Approximant: " << *it << endl;

  // CFApprox...
  TEST_ASSERT(CFApprox(ConvertTo<BigRat>(1.4142135), BigRat(1,1)) == BigRat(1,1));
  TEST_ASSERT(CFApprox(ConvertTo<BigRat>(-1.4142135), BigRat(1,1)) == BigRat(-2,1));
  TEST_ASSERT(CFApprox(ConvertTo<BigRat>(1.4142135), ConvertTo<BigRat>(0.01)) == BigRat(17,12));
  TEST_ASSERT(CFApprox(ConvertTo<BigRat>(-1.4142135), ConvertTo<BigRat>(0.01)) == BigRat(-17,12));

  TEST_ASSERT(SimplestBigRatBetween(BigRat(1,1), BigRat(2,1)) == BigRat(1,1));
  TEST_ASSERT(SimplestBigRatBetween(BigRat(-2,1), BigRat(-1,1)) == BigRat(-1,1));
  TEST_ASSERT(SimplestBigRatBetween(BigRat(-1,1), BigRat(1,1)) == BigRat(0,1));
  TEST_ASSERT(SimplestBigRatBetween(BigRat(0,1), BigRat(1,1)) == BigRat(0,1));
  TEST_ASSERT(SimplestBigRatBetween(BigRat(-1,1), BigRat(0,1)) == BigRat(0,1));
  TEST_ASSERT(SimplestBigRatBetween(BigRat(1,1), BigRat(4,3)) == BigRat(1,1));
  TEST_ASSERT(SimplestBigRatBetween(BigRat(4,3), BigRat(2,1)) == BigRat(2,1));
  TEST_ASSERT(SimplestBigRatBetween(BigRat(-4,3), BigRat(-1,1)) == BigRat(-1,1));
  TEST_ASSERT(SimplestBigRatBetween(BigRat(-2,1), BigRat(-4,3)) == BigRat(-2,1));
  TEST_ASSERT(SimplestBigRatBetween(BigRat(4,3), BigRat(10,3)) == BigRat(2,1));
  TEST_ASSERT(SimplestBigRatBetween(BigRat(-10,3), BigRat(-4,3)) == BigRat(-2,1));
  TEST_ASSERT(SimplestBigRatBetween(BigRat(123,100), BigRat(135,100)) == BigRat(4,3));
  TEST_ASSERT(SimplestBigRatBetween(BigRat(-135,100), BigRat(-123,100)) == BigRat(-4,3));
  TEST_ASSERT(SimplestBigRatBetween(BigRat(5,2), BigRat(299,100)) == BigRat(5,2));
  TEST_ASSERT(SimplestBigRatBetween(BigRat(-299,100), BigRat(-5,2)) == BigRat(-5,2));
  TEST_ASSERT(SimplestBigRatBetween(BigRat(4181, 2584), BigRat(6765,4181)) == BigRat(4181,2584));
  TEST_ASSERT(SimplestBigRatBetween(BigRat(-6765,4181), BigRat(-4181, 2584)) == BigRat(-4181,2584));
  TEST_ASSERT(SimplestBigRatBetween(BigRat(225,157), BigRat(268,187)) == BigRat(225,157));
  TEST_ASSERT(SimplestBigRatBetween(BigRat(-268,187), BigRat(-225,157)) == BigRat(-225,157));
  TEST_ASSERT(SimplestBigRatBetween(BigRat(1393,972), BigRat(268,187)) == BigRat(268,187));
  TEST_ASSERT(SimplestBigRatBetween(BigRat(-268,187), BigRat(-1393,972)) == BigRat(-268,187));
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
