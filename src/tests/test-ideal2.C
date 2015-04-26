//   Copyright (c)  2007  John Abbott

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
#include "CoCoA/FractionField.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingDistrMPolyClean.H"
#include "CoCoA/RingDistrMPolyInlFpPP.H"
#include "CoCoA/RingDistrMPolyInlPP.H"
#include "CoCoA/RingFp.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/ideal.H"

using namespace CoCoA;

#include <iostream>
using std::cerr;
using std::endl;
#include <vector>
using std::vector;

//----------------------------------------------------------------------
// Test for ideals in particular cases
// functions: IsElem, IsPrime
// environments: RingZZ,
//----------------------------------------------------------------------

#define TEST_ASSERT(cond)   do { if (!(cond)) CoCoA::AssertionFailed(#cond, __FILE__, __LINE__); } while(0)


void TestZZ()
{
  ring ZZ = RingZZ();
  
  ideal Z = ideal(zero(ZZ));
  ideal J3 = ideal(RingElem(ZZ,3));
  ideal J35 = ideal(RingElem(ZZ,3), RingElem(ZZ,5));
  ideal J46 = ideal(RingElem(ZZ,4), RingElem(ZZ,6));
  ideal J6 = ideal(RingElem(ZZ,6));

  TEST_ASSERT( !IsElem(RingElem(ZZ, 7), Z) );
  TEST_ASSERT( IsPrime(Z) );
  TEST_ASSERT( IsPrime(J3) );
  TEST_ASSERT( IsPrime(J46) );
  TEST_ASSERT( IsPrime(intersect(J35, J46)) );
  TEST_ASSERT( !IsPrime(intersect(J3, J46)) );
  TEST_ASSERT( IsPrime(intersect(J46, Z)) );
  TEST_ASSERT( intersect(J3, J6) == J6 );
}


void program()
{
  GlobalManager CoCoAFoundations;

  TestZZ();
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
