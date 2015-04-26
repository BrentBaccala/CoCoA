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
#include "CoCoA/GlobalManager.H"
#include "CoCoA/PolyRing.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/error.H"
#include "CoCoA/ring.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"


#include <iostream>
using std::cerr;
using std::endl;

using namespace CoCoA;

// Handy macro for making assertions.
#define TEST_ASSERT(cond)   do { if (!(cond)) CoCoA::AssertionFailed(#cond, __FILE__, __LINE__); } while(0)

void TestGCD(ring R)
{
  if (!IsTrueGCDDomain(R)) return;
  TEST_ASSERT(IsZero(gcd(zero(R),zero(R))));
  TEST_ASSERT(IsOne(gcd(zero(R),one(R))));
  TEST_ASSERT(IsOne(gcd(0,one(R))));
  TEST_ASSERT(IsOne(gcd(zero(R),1)));
  TEST_ASSERT(IsOne(gcd(one(R),zero(R))));
  TEST_ASSERT(IsOne(gcd(1,zero(R))));
  TEST_ASSERT(IsOne(gcd(one(R),0)));
  TEST_ASSERT(IsOne(gcd(one(R),one(R))));
  TEST_ASSERT(IsOne(gcd(1,one(R))));
  TEST_ASSERT(IsOne(gcd(one(R),1)));

  RingElem two(R,2);
  if (!IsZero(two))
  {
    TEST_ASSERT(IsInvertible(gcd(zero(R),two)/two));
    TEST_ASSERT(IsInvertible(gcd(0,two)/two));
    TEST_ASSERT(IsInvertible(gcd(zero(R),2)/two));
    TEST_ASSERT(IsInvertible(gcd(two,two)/two));
    TEST_ASSERT(IsInvertible(gcd(2,two)/two));
    TEST_ASSERT(IsInvertible(gcd(two,2)/two));
    TEST_ASSERT(IsInvertible(gcd(two*two,two)/two));
    TEST_ASSERT(IsInvertible(gcd(4,two)/two));
    TEST_ASSERT(IsInvertible(gcd(two*two,2)/two));
  }

  if (!IsPolyRing(R)) return;
  const RingElem x = indet(R,0);
  TEST_ASSERT(IsInvertible(gcd(x,x)/x));
  TEST_ASSERT(IsInvertible(gcd(0,x)/x));
  TEST_ASSERT(IsInvertible(gcd(x,0)/x));
}

void TestLCM(ring R)
{
  if (!IsTrueGCDDomain(R)) return;
  TEST_ASSERT(IsZero(lcm(zero(R),zero(R))));
  TEST_ASSERT(IsZero(lcm(zero(R),one(R))));
  TEST_ASSERT(IsZero(lcm(0,one(R))));
  TEST_ASSERT(IsZero(lcm(zero(R),1)));
  TEST_ASSERT(IsZero(lcm(one(R),zero(R))));
  TEST_ASSERT(IsZero(lcm(1,zero(R))));
  TEST_ASSERT(IsZero(lcm(one(R),0)));
  TEST_ASSERT(IsOne(lcm(one(R),one(R))));
  TEST_ASSERT(IsOne(lcm(1,one(R))));
  TEST_ASSERT(IsOne(lcm(one(R),1)));

  RingElem two(R,2);
  if (!IsZero(two))
  {
    TEST_ASSERT(IsZero(lcm(zero(R),two)));
    TEST_ASSERT(IsZero(lcm(0,two)));
    TEST_ASSERT(IsZero(lcm(zero(R),2)));
    TEST_ASSERT(IsInvertible(lcm(two,two)/two));
    TEST_ASSERT(IsInvertible(lcm(2,two)/two));
    TEST_ASSERT(IsInvertible(lcm(two,2)/two));
    TEST_ASSERT(IsInvertible(lcm(two*two,two)/4));
    TEST_ASSERT(IsInvertible(lcm(4,two)/4));
    TEST_ASSERT(IsInvertible(lcm(two*two,2)/4));
  }

  if (!IsPolyRing(R)) return;
  const RingElem x = indet(R,0);
  TEST_ASSERT(IsInvertible(lcm(x,x)/x));
  TEST_ASSERT(IsZero(lcm(0,x)));
  TEST_ASSERT(IsZero(lcm(x,0)));
}


void program()
{
  // Test GCD and LCM in various rings.
  GlobalManager CoCoAFoundations;

  ring R1 = RingZZ();
  TestGCD(R1);   TestLCM(R1);

  ring R2 = RingQQ();
  TestGCD(R2);   TestLCM(R2);

  ring R3 = NewPolyRing(RingQQ(),1);
  TestGCD(R3);   TestLCM(R3);

  ring R4 = NewPolyRing(RingQQ(),2);
  TestGCD(R4);   TestLCM(R4);

  ring R5 = NewPolyRing(RingZZ(),1);
  TestGCD(R5);   TestLCM(R5);

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
