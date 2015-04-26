//   Copyright (c)  2013  Anna M. Bigatti

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
#include "CoCoA/PolyRing.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/error.H"
#include "CoCoA/symbol.H"
#include "CoCoA/VectorOperations.H"
#include "CoCoA/utils.H"
//#include "CoCoA/assert.H"  // for using TEST_ASSERT



#include <iostream>
using std::cerr;
using std::endl;

using namespace CoCoA;

// Handy macro for making assertions.
#define TEST_ASSERT(cond)   do { if (!(cond)) CoCoA::AssertionFailed(#cond, __FILE__, __LINE__); } while(0)


void test_Qi()
{
  //  std::cout << "------- Qi --------" << std::endl;
  PolyRing Qi = NewPolyRing(RingQQ(), symbols("i"));
  QuotientRing QR = NewQuotientRing(Qi, ideal(power(indet(Qi,0),2)+1));
  RingElem i(QR, symbol("i"));
  RingElem z1 = 2+3*i;
  //  std::cout << "z1 = " << z1 << std::endl;
  TEST_ASSERT(IsInvertible(z1));
  TEST_ASSERT(!IsZeroDivisor(z1));
  //  std::cout << "inverse(z1) = " << inverse(z1) << std::endl;
  TEST_ASSERT(IsIntegralDomain(QR));
  TEST_ASSERT(IsField(QR));
  ideal I(2+i, 3*i);
  //  std::cout << "I = " << I << std::endl;
  TEST_ASSERT(IsOne(I));
  TEST_ASSERT(len(TidyGens(I)) == 1);
  // std::cout << "IsPrime(I) = " << IsPrime(I) << std::endl;

  //  PolyRing QR_x = NewPolyRing(Q(i), symbols("x","y"));
}


void test_xy()
{
  //  std::cout << "------- Qxy --------" << std::endl;
  PolyRing Qxy = NewPolyRing(RingQQ(), symbols("x","y"));
  QuotientRing QR = NewQuotientRing(Qxy, ideal(product(indets(Qxy))));
  RingElem x(QR, symbol("x"));
  RingElem y(QR, symbol("y"));
  RingElem z1 = x;
  //  std::cout << "z1 = " << z1 << std::endl;
  TEST_ASSERT(!IsInvertible(z1));
  TEST_ASSERT(IsZeroDivisor(z1));
  z1 = x+y+1;
  TEST_ASSERT(!IsInvertible(z1));
  TEST_ASSERT(!IsZeroDivisor(z1));
  //  std::cout << "inverse(z1) = " << inverse(z1) << std::endl;
  //  std::cout << "IsIntegralDomain(QR) = " << IsIntegralDomain(QR) << std::endl;
  //  std::cout << "IsField(QR) = " << IsField(QR) << std::endl;
  ideal I(2+x);
  //  std::cout << "I = " << I << std::endl;
  TEST_ASSERT(!IsOne(z1));
  //  std::cout << "TidyGens(I) = " << TidyGens(I) << std::endl;
  TEST_ASSERT(len(TidyGens(I)) == 2);
  // std::cout << "IsPrime(I) = " << IsPrime(I) << std::endl;

  //  PolyRing QR_x = NewPolyRing(Q(i), symbols("x","y"));
}


void program()
{
  GlobalManager CoCoAFoundations;
  std::cout << std::boolalpha;

  test_Qi();
  test_xy();
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
