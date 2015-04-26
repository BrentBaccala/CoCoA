//   Copyright (c)  2010  John Abbott

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
#include "CoCoA/RingTwinFloat.H"
#include "CoCoA/VectorOperations.H"

using namespace CoCoA;

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <vector>
using std::vector;

#define TEST_ASSERT(cond)   do { if (!(cond)) CoCoA::AssertionFailed(#cond, __FILE__, __LINE__); } while(0)

// This test illustrates the "soft transitions" from false to uncertain and from
// uncertain to true in the equality test between twin floats.

void trial(long BitPrec)
{
  const int IterMax = 1000;
  const ring RR = NewRingTwinFloat(BitPrec);

  cout << "RING is " << RR << endl;
  BigRat eps = BigRat(1,1);
  for (long i = 1; i < 99; ++i)
  {
    eps /= 2;

    // Subprogram A
    vector<int> SubprogA(3);
    for (int j=0; j < IterMax; ++j)
    {
      try
      {
        const RingElem X(RR, 1);
        if (X+eps == X)
          ++SubprogA[0];
        else
          ++SubprogA[2];
      }
      catch (const RingTwinFloat::InsufficientPrecision&)
      {
        ++SubprogA[1];
      }
    }

    // Subprogram B
    vector<int> SubprogB(3);
    for (int j=0; j < IterMax; ++j)
    {
      try
      {
        const RingElem X(RR, 1);
        const RingElem Y(RR, 1+eps);
        if (X == Y)
          ++SubprogB[0];
        else
          ++SubprogB[2];
      }
      catch (const RingTwinFloat::InsufficientPrecision&)
      {
        ++SubprogB[1];
      }
    }

    // Subprogram C
    vector<int> SubprogC(3);
    for (int j=0; j < IterMax; ++j)
    {
      try
      {
        const RingElem X(RR, 1);
        const RingElem Y(RR, 1+eps);
        if (IsZero(X-Y))
          ++SubprogC[0];
        else
          ++SubprogC[2];
      }
      catch (const RingTwinFloat::InsufficientPrecision&)
      {
        ++SubprogC[1];
      }
    }


    cout << "precision=" << BitPrec
                   << "  iter=" << i
                   << "  SubprogA=" << SubprogA
                   << "  SubprogB=" << SubprogB
                   << "  SubprogC=" << SubprogC
                   << endl;
  }
  cout << "---------------------------------" << endl;
}


void program()
{
  GlobalManager CoCoAFoundations;

  for (long prec=16; prec <= 32; prec+=16)
    trial(prec);

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
