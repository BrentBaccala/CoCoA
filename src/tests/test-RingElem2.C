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


#include "CoCoA/BigRat.H"
#include "CoCoA/BuildInfo.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/error.H"
#include "CoCoA/IntOperations.H"
#include "CoCoA/ring.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"


#include <iostream>
using std::cerr;
using std::endl;

using namespace CoCoA;

// Handy macro for making assertions.
#define TEST_ASSERT(cond)   do { if (!(cond)) CoCoA::AssertionFailed(#cond, __FILE__, __LINE__); } while(0)


void program()
{
  GlobalManager CoCoAFoundations;

  // Test default ring, and automatic mapping into ring when assigning.
  ring R;  // default is ZZ
  TEST_ASSERT(R == RingZZ());
  RingElem x; // default is 0 in ZZ
  TEST_ASSERT(owner(x) == RingZZ());
  x = 1;  // will not change ring
  TEST_ASSERT(owner(x) == RingZZ());
  x = power(10,100);  // will not change ring
  TEST_ASSERT(owner(x) == RingZZ());
  x = BigRat(2,1);    // will not change ring
  TEST_ASSERT(owner(x) == RingZZ());

  // We expect error if rational is not an integer
  try { x = BigRat(2,3); TEST_ASSERT(!"NEVER GET HERE!"); }
  catch (const ErrorInfo& err) { TEST_ASSERT(err == ERR::EmbedBigRatFailed); }

  R = NewFractionField(R);  // now R is QQ
  TEST_ASSERT(R == RingQQ());
  TEST_ASSERT(owner(x) != R);

  RingElem y(R); // y belongs to QQ
  TEST_ASSERT(owner(y) == RingQQ());
  y = BigRat(2,3);   // works because y is in QQ
  TEST_ASSERT(owner(y) == RingQQ());

  TEST_ASSERT(owner(x) != owner(y));
  swap(x, y);       // now x is in QQ, and y in ZZ
  TEST_ASSERT(owner(x) == RingQQ());
  TEST_ASSERT(owner(y) == RingZZ());

  x = y;             // x is now back in ZZ
  TEST_ASSERT(owner(x) == owner(y));
  TEST_ASSERT(owner(x) == RingZZ());
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
