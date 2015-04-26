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
#include "CoCoA/RingTwinFloat.H"
#include "CoCoA/RingFp.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/BigInt.H"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <vector>
using std::vector;

using namespace CoCoA;

#define TEST_ASSERT(cond)   do { if (!(cond)) CoCoA::AssertionFailed(#cond, __FILE__, __LINE__); } while(0)


void TestInteger(const RingElem& x)
{
  BigInt n;
  try
  {
  if (IsInteger(n, x))
    cout << "The ring element " << x << " is the image of the integer " << n << endl;
  else
    cout << "The ring element " << x << " is not the image of an integer." << endl;
  }
  catch (const RingTwinFloat::InsufficientPrecision&)
  {
    cout << "Failed to decide whether the ring element " << x << " is or is not the image of an integer." << endl;
  }

}


void trial(const ring& R)
{

  cout << "START OF TEST for IsInteger over " << R << endl;

  const RingElem x(R, 29);
  TestInteger(x);
  if (IsDivisible(x, 2)) TestInteger(x/2);
  if (IsDivisible(x, 3)) TestInteger(x/3);
  if (IsDivisible(x, x+1)) TestInteger(x/(x+1));
  TestInteger(power(x, 10));

  cout << "END OF TEST for IsInteger over " << R << endl << endl << endl;
}

void program()
{
  GlobalManager CoCoAFoundations(UseNonNegResidues);

  trial(RingZZ());
  trial(RingQQ());
  trial(NewZZmod(2));
  trial(NewZZmod(3));
  trial(NewZZmod(32003));   // large prime
  trial(NewZZmod(1000003)); // larger prime
  trial(NewZZmod(6*29)); // ring has zero divisors
  trial(NewZZmod(1048576)); // ring has zero divisors
  trial(NewRingTwinFloat(32));
  trial(NewRingTwinFloat(16)); // this will cause a "Failed..." message to be printed for 29^10
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
