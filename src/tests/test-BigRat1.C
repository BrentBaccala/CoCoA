//   Copyright (c)  2009  John Abbott

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


#include "CoCoA/BigRat.H"
#include "CoCoA/BuildInfo.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/IntOperations.H"
#include "CoCoA/error.H"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

using namespace CoCoA;

#define TEST_ASSERT(cond)   do { if (!(cond)) CoCoA::AssertionFailed(#cond, __FILE__, __LINE__); } while(0)

void program()
{
  // This test does virtually nothing, but is a handy template if you want
  // to create your own test code: just copy this file and add your code
  // after the line below -- remember that you must not use CoCoALib functions
  // without first creating a GlobalManager.
  GlobalManager CoCoAFoundations;

  const int a=13;
  const int b=21;

  cout << 101 << " = " << BigRat(101,1) << endl;
  cout << power(101, 9) << " = " << BigRat(power(101,9),1) << endl;

  // Test the ctors.
  BigRat frac(a,b);
  TEST_ASSERT(frac == BigRat(a,BigInt(b)));
  TEST_ASSERT(frac == BigRat(BigInt(a),b));
  TEST_ASSERT(frac == BigRat(BigInt(a),BigInt(b)));

  // Generators are negative
  TEST_ASSERT(frac == BigRat(-a,-b));
  TEST_ASSERT(frac == BigRat(-a,BigInt(-b)));
  TEST_ASSERT(frac == BigRat(BigInt(-a),-b));
  TEST_ASSERT(frac == BigRat(BigInt(-a),BigInt(-b)));

  // Generators are not coprime, but positive.
  TEST_ASSERT(frac == BigRat(a*b, b*b));
  TEST_ASSERT(frac == BigRat(a*b, BigInt(b*b)));
  TEST_ASSERT(frac == BigRat(BigInt(a*b), b*b));
  TEST_ASSERT(frac == BigRat(BigInt(a*b), BigInt(b*b)));

  // Generators are not coprime, and negative.
  TEST_ASSERT(frac == BigRat(-a*b, -b*b));
  TEST_ASSERT(frac == BigRat(-a*b, BigInt(-b*b)));
  TEST_ASSERT(frac == BigRat(BigInt(-a*b), -b*b));
  TEST_ASSERT(frac == BigRat(BigInt(-a*b), BigInt(-b*b)));

  // Now check that a/b - a/b gives zero.
  TEST_ASSERT(IsZero(frac-frac));
  TEST_ASSERT(IsZero(frac + BigRat(-a,b)));
  TEST_ASSERT(IsZero(frac + BigRat(a,-b)));
  TEST_ASSERT(IsZero(frac + BigRat(-a,BigInt(b))));
  TEST_ASSERT(IsZero(frac + BigRat(a,BigInt(-b))));
  TEST_ASSERT(IsZero(frac + BigRat(BigInt(-a),b)));
  TEST_ASSERT(IsZero(frac + BigRat(BigInt(a),-b)));
  TEST_ASSERT(IsZero(frac + BigRat(BigInt(-a),BigInt(b))));
  TEST_ASSERT(IsZero(frac + BigRat(BigInt(a),BigInt(-b))));

  // Non trivial GCD between the ctor args.
  TEST_ASSERT(IsZero(frac + BigRat(-a*b,b*b)));
  TEST_ASSERT(IsZero(frac + BigRat(a*b,-b*b)));
  TEST_ASSERT(IsZero(frac + BigRat(-a*b,BigInt(b*b))));
  TEST_ASSERT(IsZero(frac + BigRat(a*b,BigInt(-b*b))));
  TEST_ASSERT(IsZero(frac + BigRat(BigInt(-a*b),b*b)));
  TEST_ASSERT(IsZero(frac + BigRat(BigInt(a*b),-b*b)));
  TEST_ASSERT(IsZero(frac + BigRat(BigInt(-a*b),BigInt(b*b))));
  TEST_ASSERT(IsZero(frac + BigRat(BigInt(a*b),BigInt(-b*b))));


  BigRat frac2;
  TEST_ASSERT(IsZero(frac2));
  frac2 += frac;
  TEST_ASSERT(frac2 == frac);
  frac2 -= frac;
  TEST_ASSERT(IsZero(frac2));
  frac2 = frac;
  frac2 *= frac;
  TEST_ASSERT(frac2 == BigRat(a*a, b*b));
  frac2 /= frac2;
  TEST_ASSERT(IsOne(frac2));

  frac2 += BigInt(99);
  TEST_ASSERT(frac2 == 100);
  frac2 -= BigInt(-100);
  TEST_ASSERT(frac2 == 200);
  frac2 *= BigInt(-5);
  TEST_ASSERT(frac2 == -1000);
  frac2 /= BigInt(1000);
  TEST_ASSERT(IsMinusOne(frac2));

  frac2 += 11;
  TEST_ASSERT(frac2 == 10);
  frac2 -= -6;
  TEST_ASSERT(frac2 == 16);
  frac2 *= -5;
  TEST_ASSERT(frac2 == -80);
  frac2 /= 9;
  TEST_ASSERT(frac2 == BigRat(80,-9));

  TEST_ASSERT(BigRat(-71,9) == ++frac2);
  TEST_ASSERT(frac2++ == BigRat(71,-9));
  TEST_ASSERT(frac2 == BigRat(62,-9));
  TEST_ASSERT(--frac2 == BigRat(-71,9));
  TEST_ASSERT(BigRat(71,-9) == frac2--);
  TEST_ASSERT(frac2 == BigRat(-80,9));

  swap(frac, frac2);
  TEST_ASSERT(frac == BigRat(-80,9));
  TEST_ASSERT(frac2 == BigRat(a,b));

  TEST_ASSERT(num(frac) == -80);
  TEST_ASSERT(den(frac) == 9);
  TEST_ASSERT(num(abs(frac)) == 80);
  TEST_ASSERT(den(abs(frac)) == 9);
  TEST_ASSERT(num(-frac) == 80);
  TEST_ASSERT(den(-frac) == 9);

  TEST_ASSERT(frac+frac2 == BigRat(9*a-80*b,9*b));
  TEST_ASSERT(frac-frac2 == BigRat(-9*a-80*b,9*b));
  TEST_ASSERT(frac*frac2 == BigRat(-80*a,9*b));
  TEST_ASSERT(frac/frac2 == BigRat(-80*b,9*a));

  TEST_ASSERT(frac+BigInt(10) == BigRat(10,9));
  TEST_ASSERT(frac-BigInt(10) == BigRat(-170,9));
  TEST_ASSERT(frac*BigInt(10) == BigRat(-800,9));
  TEST_ASSERT(frac/BigInt(10) == BigRat(-8,9));

  TEST_ASSERT(BigInt(-10)+frac == BigRat(-170,9));
  TEST_ASSERT(BigInt(-10)-frac == BigRat(-10,9));
  TEST_ASSERT(BigInt(-10)*frac == BigRat(800,9));
  TEST_ASSERT(BigInt(-10)/frac == BigRat(9,8));

  TEST_ASSERT(frac+10 == BigRat(10,9));
  TEST_ASSERT(frac-10 == BigRat(-170,9));
  TEST_ASSERT(frac*10 == BigRat(-800,9));
  TEST_ASSERT(frac/10 == BigRat(-8,9));

  TEST_ASSERT(10+frac == BigRat(10,9));
  TEST_ASSERT(10-frac == BigRat(170,9));
  TEST_ASSERT(10*frac == BigRat(-800,9));
  TEST_ASSERT(10/frac == BigRat(-9,8));

  BigRat frac3 = power(frac2, 17);
  TEST_ASSERT(frac3 == BigRat(power(a,17), power(b,17)));

  frac3 = power(frac2, BigInt(255));
  TEST_ASSERT(frac3 == BigRat(power(a,255), power(b,255)));

  TEST_ASSERT(!IsZero(frac));
  TEST_ASSERT(!IsOne(frac));
  TEST_ASSERT(IsOne(frac/frac));
  TEST_ASSERT(!IsMinusOne(frac/frac));
  TEST_ASSERT(IsMinusOne(frac*((-1)/frac)));
  TEST_ASSERT(sign(frac) == -1 && sign(frac2) == 1);

  // Comparisons
  TEST_ASSERT(cmp(frac, frac) == 0);
  TEST_ASSERT(cmp(frac, -frac) < 0);
  TEST_ASSERT(cmp(-frac, frac) > 0);

  TEST_ASSERT(cmp(frac, -BigInt(9)) > 0);
  TEST_ASSERT(cmp(frac, -BigInt(8)) < 0);
  TEST_ASSERT(cmp(BigInt(9), -frac) > 0);
  TEST_ASSERT(cmp(BigInt(8), -frac) < 0);

  TEST_ASSERT(cmp(frac, -9) > 0);
  TEST_ASSERT(cmp(frac, -8) < 0);
  TEST_ASSERT(cmp(9, -frac) > 0);
  TEST_ASSERT(cmp(8, -frac) < 0);

  TEST_ASSERT(frac == frac);
  TEST_ASSERT(!(frac != frac));
  TEST_ASSERT(!(frac == frac2));
  TEST_ASSERT(frac != frac2);
  TEST_ASSERT(!(frac < frac));
  TEST_ASSERT(frac < frac2);
  TEST_ASSERT(frac <= frac);
  TEST_ASSERT(frac <= frac2);
  TEST_ASSERT(!(frac > frac));
  TEST_ASSERT(frac2 > frac);
  TEST_ASSERT(frac >= frac);
  TEST_ASSERT(frac2 >= frac);

  TEST_ASSERT(!(frac == BigInt(-9)));
  TEST_ASSERT(frac != BigInt(-9));
  TEST_ASSERT(frac > BigInt(-9));
  TEST_ASSERT(frac >= BigInt(-9));
  TEST_ASSERT(frac < BigInt(-8));
  TEST_ASSERT(frac <= BigInt(-8));

  TEST_ASSERT(!(BigInt(-9) == frac));
  TEST_ASSERT(BigInt(-9) != frac);
  TEST_ASSERT(BigInt(-8) > frac);
  TEST_ASSERT(BigInt(-8) >= frac);
  TEST_ASSERT(-BigInt(9) < frac);
  TEST_ASSERT(-BigInt(9) <= frac);

  TEST_ASSERT(!(frac == -9));
  TEST_ASSERT(frac != -9);
  TEST_ASSERT(frac > -9);
  TEST_ASSERT(frac >= -9);
  TEST_ASSERT(frac < -8);
  TEST_ASSERT(frac <= -8);

  TEST_ASSERT(!(-9 == frac));
  TEST_ASSERT(-9 != frac);
  TEST_ASSERT(-8 > frac);
  TEST_ASSERT(-8 >= frac);
  TEST_ASSERT(-9 < frac);
  TEST_ASSERT(-9 <= frac);

  TEST_ASSERT(round(frac) == -9);
  TEST_ASSERT(round(BigRat(1,2)) == RoundDiv(1,2));
  TEST_ASSERT(round(BigRat(-1,2)) == RoundDiv(-1,2));
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
