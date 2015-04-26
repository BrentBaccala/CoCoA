//   Copyright (c)  2013  John Abbott

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
#include "CoCoA/CanonicalHom.H"
#include "CoCoA/FractionField.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/PolyRing.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/SparsePolyRing.H"
//#include "CoCoA/assert.H"  // for using TEST_ASSERT
#include "CoCoA/error.H"
#include "CoCoA/symbol.H"


#include <iostream>
using std::cerr;
using std::endl;

using namespace CoCoA;

// Handy macro for making assertions.
#define TEST_ASSERT(cond)   do { if (!(cond)) CoCoA::AssertionFailed(#cond, __FILE__, __LINE__); } while(0)

void program()
{
  // This test does virtually nothing, but is a handy template if you want
  // to create your own test code: just copy this file and add your code
  // after the line below -- remember that you must not use CoCoALib functions
  // without first creating a GlobalManager.
  GlobalManager CoCoAFoundations;

  ring ZZ = RingZZ();
  FractionField QQ = RingQQ();
  ring Q = NewFractionField(ZZ);
  TEST_ASSERT(Q == QQ);
  RingElem q(QQ,2);
  q /= 3;
  TEST_ASSERT(num(q) == 2 && den(q) == 3);

  ring QQx = NewPolyRing(QQ, symbols("x"));
  FractionField FrFQQx = NewFractionField(QQx);
  RingElem x(FrFQQx, symbol("x"));
  RingElem f = (x-1)/(x+1);
  RingHom embed = EmbeddingHom(FrFQQx);
  TEST_ASSERT(embed(num(f)) == x-1);
  TEST_ASSERT(embed(den(f)) == x+1);

  RingElem f1 = deriv(f, x);
  TEST_ASSERT(f1 == 2/power(x+1,2));

  RingHom ZZtoQQx = CanonicalHom(ZZ, QQx);
  RingHom QQtoQQx = InducedHom(QQ, ZZtoQQx);
  RingElem image = QQtoQQx(q);
  TEST_ASSERT(image == BigRat(2,3));
  RingHom QQtoFrF = InducedHom(QQ, CanonicalHom(ZZ, FrFQQx));
  TEST_ASSERT(QQtoFrF(q) == BigRat(2,3));
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
