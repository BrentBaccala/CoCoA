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
#include "CoCoA/assert.H"
#include "CoCoA/error.H"
#include "CoCoA/IntOperations.H"


#include <iostream>
using std::cerr;
using std::endl;

using namespace CoCoA;

// Handy macro for making assertions.
#define TEST_ASSERT(cond)   do { if (!(cond)) CoCoA::AssertionFailed(#cond, __FILE__, __LINE__); } while(0)

void program()
{
  // Test for integer division and remainder.
  GlobalManager CoCoAFoundations;

  const int n1 = 99;
  const BigInt N1(99);
  const int n2 = 101;
  const BigInt N2(101);
  const int d = 100;
  const BigInt D(100);
  const int half = d/2;
  const BigInt HALF(half);

  ///////////////////////////////////////////////////////
  // INTEGER DIVISION

  // MachineInt & MachineInt
  TEST_ASSERT(n1/d == 0);
  TEST_ASSERT((-n1)/d == 0);
  TEST_ASSERT(n1/(-d) == 0);
  TEST_ASSERT((-n1)/(-d) == 0);
  TEST_ASSERT(n2/d == 1);
  TEST_ASSERT((-n2)/d == -1);
  TEST_ASSERT(n2/(-d) == -1);
  TEST_ASSERT((-n2)/(-d) == 1);

  // BigInt & MachineInt
  TEST_ASSERT(N1/d == 0);
  TEST_ASSERT((-N1)/d == 0);
  TEST_ASSERT(N1/(-d) == 0);
  TEST_ASSERT((-N1)/(-d) == 0);
  TEST_ASSERT(N2/d == 1);
  TEST_ASSERT((-N2)/d == -1);
  TEST_ASSERT(N2/(-d) == -1);
  TEST_ASSERT((-N2)/(-d) == 1);

  // MachineInt & BigInt
  TEST_ASSERT(n1/D == 0);
  TEST_ASSERT((-n1)/D == 0);
  TEST_ASSERT(n1/(-D) == 0);
  TEST_ASSERT((-n1)/(-D) == 0);
  TEST_ASSERT(n2/D == 1);
  TEST_ASSERT((-n2)/D == -1);
  TEST_ASSERT(n2/(-D) == -1);
  TEST_ASSERT((-n2)/(-D) == 1);

  // BigInt & BigInt
  TEST_ASSERT(N1/D == 0);
  TEST_ASSERT((-N1)/D == 0);
  TEST_ASSERT(N1/(-D) == 0);
  TEST_ASSERT((-N1)/(-D) == 0);
  TEST_ASSERT(N2/D == 1);
  TEST_ASSERT((-N2)/D == -1);
  TEST_ASSERT(N2/(-D) == -1);
  TEST_ASSERT((-N2)/(-D) == 1);

  ///////////////////////////////////////////////////////
  // REMAINDER

  // MachineInt & MachineInt
  TEST_ASSERT(n1%d == 99);
//   TEST_ASSERT((-n1)%d == 1); // depends on compiler & operating system
//   TEST_ASSERT(n1%(-d) == 99); // depends on compiler & operating system
//   TEST_ASSERT((-n1)%(-d) == 1); // depends on compiler & operating system
  TEST_ASSERT(n2%d == 1);
//   TEST_ASSERT((-n2)%d == 99); // depends on compiler & operating system
//   TEST_ASSERT(n2%(-d) == 1); // depends on compiler & operating system
//   TEST_ASSERT((-n2)%(-d) == 99); // depends on compiler & operating system

  // BigInt & MachineInt
  TEST_ASSERT(N1%d == 99);
  TEST_ASSERT((-N1)%d == -99);
  TEST_ASSERT(N1%(-d) == 99);
  TEST_ASSERT((-N1)%(-d) == -99);
  TEST_ASSERT(N2%d == 1);
  TEST_ASSERT((-N2)%d == -1);
  TEST_ASSERT(N2%(-d) == 1);
  TEST_ASSERT((-N2)%(-d) == -1);

  // MachineInt & BigInt
  TEST_ASSERT(n1%D == 99);
  TEST_ASSERT((-n1)%D == -99);
  TEST_ASSERT(n1%(-D) == 99);
  TEST_ASSERT((-n1)%(-D) == -99);
  TEST_ASSERT(n2%D == 1);
  TEST_ASSERT((-n2)%D == -1);
  TEST_ASSERT(n2%(-D) == 1);
  TEST_ASSERT((-n2)%(-D) == -1);

  // BigInt & BigInt
  TEST_ASSERT(N1%D == 99);
  TEST_ASSERT((-N1)%D == -99);
  TEST_ASSERT(N1%(-D) == 99);
  TEST_ASSERT((-N1)%(-D) == -99);
  TEST_ASSERT(N2%D == 1);
  TEST_ASSERT((-N2)%D == -1);
  TEST_ASSERT(N2%(-D) == 1);
  TEST_ASSERT((-N2)%(-D) == -1);

  ///////////////////////////////////////////////////////
  // LeastNNegRemainder

  TEST_ASSERT(LeastNNegRemainder(d,d) == 0);
  TEST_ASSERT(LeastNNegRemainder(-d,d) == 0);
  TEST_ASSERT(LeastNNegRemainder(d,-d) == 0);
  TEST_ASSERT(LeastNNegRemainder(-d,-d) == 0);

  TEST_ASSERT(LeastNNegRemainder(D,d) == 0);
  TEST_ASSERT(LeastNNegRemainder(-D,d) == 0);
  TEST_ASSERT(LeastNNegRemainder(D,-d) == 0);
  TEST_ASSERT(LeastNNegRemainder(-D,-d) == 0);

  TEST_ASSERT(LeastNNegRemainder(d,D) == 0);
  TEST_ASSERT(LeastNNegRemainder(-d,D) == 0);
  TEST_ASSERT(LeastNNegRemainder(d,-D) == 0);
  TEST_ASSERT(LeastNNegRemainder(-d,-D) == 0);

  TEST_ASSERT(LeastNNegRemainder(D,D) == 0);
  TEST_ASSERT(LeastNNegRemainder(-D,D) == 0);
  TEST_ASSERT(LeastNNegRemainder(D,-D) == 0);
  TEST_ASSERT(LeastNNegRemainder(-D,-D) == 0);


  // MachineInt & MachineInt
  TEST_ASSERT(LeastNNegRemainder(n1,d) == 99);
  TEST_ASSERT(LeastNNegRemainder(-n1,d) == 1);
  TEST_ASSERT(LeastNNegRemainder(n1,-d) == 99);
  TEST_ASSERT(LeastNNegRemainder(-n1,-d) == 1);
  TEST_ASSERT(LeastNNegRemainder(n2,d) == 1);
  TEST_ASSERT(LeastNNegRemainder(-n2,d) == 99);
  TEST_ASSERT(LeastNNegRemainder(n2,-d) == 1);
  TEST_ASSERT(LeastNNegRemainder(-n2,-d) == 99);
  TEST_ASSERT(LeastNNegRemainder(half,d) == half);
  TEST_ASSERT(LeastNNegRemainder(-half,d) == half);
  TEST_ASSERT(LeastNNegRemainder(half,-d) == half);
  TEST_ASSERT(LeastNNegRemainder(-half,-d) == half);
  TEST_ASSERT(LeastNNegRemainder(half+1,d) == half+1);
  TEST_ASSERT(LeastNNegRemainder(-half-1,d) == half-1);
  TEST_ASSERT(LeastNNegRemainder(half+1,-d) == half+1);
  TEST_ASSERT(LeastNNegRemainder(-half-1,-d) == half-1);
  TEST_ASSERT(LeastNNegRemainder(half-1,d) == half-1);
  TEST_ASSERT(LeastNNegRemainder(-half+1,d) == half+1);
  TEST_ASSERT(LeastNNegRemainder(half-1,-d) == half-1);
  TEST_ASSERT(LeastNNegRemainder(-half+1,-d) == half+1);

  // BigInt & MachineInt
  TEST_ASSERT(LeastNNegRemainder(N1,d) == 99);
  TEST_ASSERT(LeastNNegRemainder(-N1,d) == 1);
  TEST_ASSERT(LeastNNegRemainder(N1,-d) == 99);
  TEST_ASSERT(LeastNNegRemainder(-N1,-d) == 1);
  TEST_ASSERT(LeastNNegRemainder(N2,d) == 1);
  TEST_ASSERT(LeastNNegRemainder(-N2,d) == 99);
  TEST_ASSERT(LeastNNegRemainder(N2,-d) == 1);
  TEST_ASSERT(LeastNNegRemainder(-N2,-d) == 99);
  TEST_ASSERT(LeastNNegRemainder(HALF,d) == HALF);
  TEST_ASSERT(LeastNNegRemainder(-HALF,d) == HALF);
  TEST_ASSERT(LeastNNegRemainder(HALF,-d) == HALF);
  TEST_ASSERT(LeastNNegRemainder(-HALF,-d) == HALF);
  TEST_ASSERT(LeastNNegRemainder(HALF+1,d) == HALF+1);
  TEST_ASSERT(LeastNNegRemainder(-HALF-1,d) == HALF-1);
  TEST_ASSERT(LeastNNegRemainder(HALF+1,-d) == HALF+1);
  TEST_ASSERT(LeastNNegRemainder(-HALF-1,-d) == HALF-1);
  TEST_ASSERT(LeastNNegRemainder(HALF-1,d) == HALF-1);
  TEST_ASSERT(LeastNNegRemainder(-HALF+1,d) == HALF+1);
  TEST_ASSERT(LeastNNegRemainder(HALF-1,-d) == HALF-1);
  TEST_ASSERT(LeastNNegRemainder(-HALF+1,-d) == HALF+1);

  // MachineInt & BigInt
  TEST_ASSERT(LeastNNegRemainder(n1,D) == 99);
  TEST_ASSERT(LeastNNegRemainder(-n1,D) == 1);
  TEST_ASSERT(LeastNNegRemainder(n1,-D) == 99);
  TEST_ASSERT(LeastNNegRemainder(-n1,-D) == 1);
  TEST_ASSERT(LeastNNegRemainder(n2,D) == 1);
  TEST_ASSERT(LeastNNegRemainder(-n2,D) == 99);
  TEST_ASSERT(LeastNNegRemainder(n2,-D) == 1);
  TEST_ASSERT(LeastNNegRemainder(-n2,-D) == 99);
  TEST_ASSERT(LeastNNegRemainder(half,D) == half);
  TEST_ASSERT(LeastNNegRemainder(-half,D) == half);
  TEST_ASSERT(LeastNNegRemainder(half,-D) == half);
  TEST_ASSERT(LeastNNegRemainder(-half,-D) == half);
  TEST_ASSERT(LeastNNegRemainder(half+1,D) == half+1);
  TEST_ASSERT(LeastNNegRemainder(-half-1,D) == half-1);
  TEST_ASSERT(LeastNNegRemainder(half+1,-D) == half+1);
  TEST_ASSERT(LeastNNegRemainder(-half-1,-D) == half-1);
  TEST_ASSERT(LeastNNegRemainder(half-1,D) == half-1);
  TEST_ASSERT(LeastNNegRemainder(-half+1,D) == half+1);
  TEST_ASSERT(LeastNNegRemainder(half-1,-D) == half-1);
  TEST_ASSERT(LeastNNegRemainder(-half+1,-D) == half+1);

  // BigInt & BigInt
  TEST_ASSERT(LeastNNegRemainder(N1,D) == 99);
  TEST_ASSERT(LeastNNegRemainder(-N1,D) == 1);
  TEST_ASSERT(LeastNNegRemainder(N1,-D) == 99);
  TEST_ASSERT(LeastNNegRemainder(-N1,-D) == 1);
  TEST_ASSERT(LeastNNegRemainder(N2,D) == 1);
  TEST_ASSERT(LeastNNegRemainder(-N2,D) == 99);
  TEST_ASSERT(LeastNNegRemainder(N2,-D) == 1);
  TEST_ASSERT(LeastNNegRemainder(-N2,-D) == 99);
  TEST_ASSERT(LeastNNegRemainder(HALF,D) == HALF);
  TEST_ASSERT(LeastNNegRemainder(-HALF,D) == HALF);
  TEST_ASSERT(LeastNNegRemainder(HALF,-D) == HALF);
  TEST_ASSERT(LeastNNegRemainder(-HALF,-D) == HALF);
  TEST_ASSERT(LeastNNegRemainder(HALF+1,D) == HALF+1);
  TEST_ASSERT(LeastNNegRemainder(-HALF-1,D) == HALF-1);
  TEST_ASSERT(LeastNNegRemainder(HALF+1,-D) == HALF+1);
  TEST_ASSERT(LeastNNegRemainder(-HALF-1,-D) == HALF-1);
  TEST_ASSERT(LeastNNegRemainder(HALF-1,D) == HALF-1);
  TEST_ASSERT(LeastNNegRemainder(-HALF+1,D) == HALF+1);
  TEST_ASSERT(LeastNNegRemainder(HALF-1,-D) == HALF-1);
  TEST_ASSERT(LeastNNegRemainder(-HALF+1,-D) == HALF+1);

  ///////////////////////////////////////////////////////
  // SymmRemainder

  TEST_ASSERT(SymmRemainder(d,d) == 0);
  TEST_ASSERT(SymmRemainder(-d,d) == 0);
  TEST_ASSERT(SymmRemainder(d,-d) == 0);
  TEST_ASSERT(SymmRemainder(-d,-d) == 0);

  TEST_ASSERT(SymmRemainder(D,d) == 0);
  TEST_ASSERT(SymmRemainder(-D,d) == 0);
  TEST_ASSERT(SymmRemainder(D,-d) == 0);
  TEST_ASSERT(SymmRemainder(-D,-d) == 0);

  TEST_ASSERT(SymmRemainder(d,D) == 0);
  TEST_ASSERT(SymmRemainder(-d,D) == 0);
  TEST_ASSERT(SymmRemainder(d,-D) == 0);
  TEST_ASSERT(SymmRemainder(-d,-D) == 0);

  TEST_ASSERT(SymmRemainder(D,D) == 0);
  TEST_ASSERT(SymmRemainder(-D,D) == 0);
  TEST_ASSERT(SymmRemainder(D,-D) == 0);
  TEST_ASSERT(SymmRemainder(-D,-D) == 0);


  // MachineInt & MachineInt
  TEST_ASSERT(SymmRemainder(n1,d) == -1);
  TEST_ASSERT(SymmRemainder(-n1,d) == 1);
  TEST_ASSERT(SymmRemainder(n1,-d) == -1);
  TEST_ASSERT(SymmRemainder(-n1,-d) == 1);
  TEST_ASSERT(SymmRemainder(n2,d) == 1);
  TEST_ASSERT(SymmRemainder(-n2,d) == -1);
  TEST_ASSERT(SymmRemainder(n2,-d) == 1);
  TEST_ASSERT(SymmRemainder(-n2,-d) == -1);
  TEST_ASSERT(SymmRemainder(half,d) == half);
  TEST_ASSERT(SymmRemainder(-half,d) == half);
  TEST_ASSERT(SymmRemainder(half,-d) == half);
  TEST_ASSERT(SymmRemainder(-half,-d) == half);
  TEST_ASSERT(SymmRemainder(half+1,d) == -half+1);
  TEST_ASSERT(SymmRemainder(-half-1,d) == half-1);
  TEST_ASSERT(SymmRemainder(half+1,-d) == -half+1);
  TEST_ASSERT(SymmRemainder(-half-1,-d) == half-1);
  TEST_ASSERT(SymmRemainder(half-1,d) == half-1);
  TEST_ASSERT(SymmRemainder(-half+1,d) == -half+1);
  TEST_ASSERT(SymmRemainder(half-1,-d) == half-1);
  TEST_ASSERT(SymmRemainder(-half+1,-d) == -half+1);

  // BigInt & MachineInt
  TEST_ASSERT(SymmRemainder(N1,d) == -1);
  TEST_ASSERT(SymmRemainder(-N1,d) == 1);
  TEST_ASSERT(SymmRemainder(N1,-d) == -1);
  TEST_ASSERT(SymmRemainder(-N1,-d) == 1);
  TEST_ASSERT(SymmRemainder(N2,d) == 1);
  TEST_ASSERT(SymmRemainder(-N2,d) == -1);
  TEST_ASSERT(SymmRemainder(N2,-d) == 1);
  TEST_ASSERT(SymmRemainder(-N2,-d) == -1);
  TEST_ASSERT(SymmRemainder(HALF,d) == HALF);
  TEST_ASSERT(SymmRemainder(-HALF,d) == HALF);
  TEST_ASSERT(SymmRemainder(HALF,-d) == HALF);
  TEST_ASSERT(SymmRemainder(-HALF,-d) == HALF);
  TEST_ASSERT(SymmRemainder(HALF+1,d) == -HALF+1);
  TEST_ASSERT(SymmRemainder(-HALF-1,d) == HALF-1);
  TEST_ASSERT(SymmRemainder(HALF+1,-d) == -HALF+1);
  TEST_ASSERT(SymmRemainder(-HALF-1,-d) == HALF-1);
  TEST_ASSERT(SymmRemainder(HALF-1,d) == HALF-1);
  TEST_ASSERT(SymmRemainder(-HALF+1,d) == -HALF+1);
  TEST_ASSERT(SymmRemainder(HALF-1,-d) == HALF-1);
  TEST_ASSERT(SymmRemainder(-HALF+1,-d) == -HALF+1);

  // MachineInt & BigInt
  TEST_ASSERT(SymmRemainder(n1,D) == -1);
  TEST_ASSERT(SymmRemainder(-n1,D) == 1);
  TEST_ASSERT(SymmRemainder(n1,-D) == -1);
  TEST_ASSERT(SymmRemainder(-n1,-D) == 1);
  TEST_ASSERT(SymmRemainder(n2,D) == 1);
  TEST_ASSERT(SymmRemainder(-n2,D) == -1);
  TEST_ASSERT(SymmRemainder(n2,-D) == 1);
  TEST_ASSERT(SymmRemainder(-n2,-D) == -1);
  TEST_ASSERT(SymmRemainder(half,D) == half);
  TEST_ASSERT(SymmRemainder(-half,D) == half);
  TEST_ASSERT(SymmRemainder(half,-D) == half);
  TEST_ASSERT(SymmRemainder(-half,-D) == half);
  TEST_ASSERT(SymmRemainder(half+1,D) == -half+1);
  TEST_ASSERT(SymmRemainder(-half-1,D) == half-1);
  TEST_ASSERT(SymmRemainder(half+1,-D) == -half+1);
  TEST_ASSERT(SymmRemainder(-half-1,-D) == half-1);
  TEST_ASSERT(SymmRemainder(half-1,D) == half-1);
  TEST_ASSERT(SymmRemainder(-half+1,D) == -half+1);
  TEST_ASSERT(SymmRemainder(half-1,-D) == half-1);
  TEST_ASSERT(SymmRemainder(-half+1,-D) == -half+1);

  // BigInt & BigInt
  TEST_ASSERT(SymmRemainder(N1,D) == -1);
  TEST_ASSERT(SymmRemainder(-N1,D) == 1);
  TEST_ASSERT(SymmRemainder(N1,-D) == -1);
  TEST_ASSERT(SymmRemainder(-N1,-D) == 1);
  TEST_ASSERT(SymmRemainder(N2,D) == 1);
  TEST_ASSERT(SymmRemainder(-N2,D) == -1);
  TEST_ASSERT(SymmRemainder(N2,-D) == 1);
  TEST_ASSERT(SymmRemainder(-N2,-D) == -1);
  TEST_ASSERT(SymmRemainder(HALF,D) == HALF);
  TEST_ASSERT(SymmRemainder(-HALF,D) == HALF);
  TEST_ASSERT(SymmRemainder(HALF,-D) == HALF);
  TEST_ASSERT(SymmRemainder(-HALF,-D) == HALF);
  TEST_ASSERT(SymmRemainder(HALF+1,D) == -HALF+1);
  TEST_ASSERT(SymmRemainder(-HALF-1,D) == HALF-1);
  TEST_ASSERT(SymmRemainder(HALF+1,-D) == -HALF+1);
  TEST_ASSERT(SymmRemainder(-HALF-1,-D) == HALF-1);
  TEST_ASSERT(SymmRemainder(HALF-1,D) == HALF-1);
  TEST_ASSERT(SymmRemainder(-HALF+1,D) == -HALF+1);
  TEST_ASSERT(SymmRemainder(HALF-1,-D) == HALF-1);
  TEST_ASSERT(SymmRemainder(-HALF+1,-D) == -HALF+1);

  // RoundDiv; incl. check that halves round correctly.
  // Implicitly assume that rounding is symm about zero!

  const int RoundHalf = 1; // halves round AWAY FROM ZERO

  TEST_ASSERT(RoundDiv(1,d) == 0);
  TEST_ASSERT(RoundDiv(-1,d) == 0);
  TEST_ASSERT(RoundDiv(half,d) == RoundHalf);
  TEST_ASSERT(RoundDiv(half-1,d) == 0);
  TEST_ASSERT(RoundDiv(half+1,d) == 1);
  TEST_ASSERT(RoundDiv(-half,d) == -RoundHalf);
  TEST_ASSERT(RoundDiv(-half+1,d) == 0);
  TEST_ASSERT(RoundDiv(-half-1,d) == -1);
  TEST_ASSERT(RoundDiv(d+half,d) == 1+RoundHalf);
  TEST_ASSERT(RoundDiv(d+half-1,d) == 1);
  TEST_ASSERT(RoundDiv(d+half+1,d) == 2);
  TEST_ASSERT(RoundDiv(-d-half,d) == -1-RoundHalf);
  TEST_ASSERT(RoundDiv(-d-half+1,d) == -1);
  TEST_ASSERT(RoundDiv(-d-half-1,d) == -2);

  TEST_ASSERT(RoundDiv(1,D) == 0);
  TEST_ASSERT(RoundDiv(-1,D) == 0);
  TEST_ASSERT(RoundDiv(half,D) == RoundHalf);
  TEST_ASSERT(RoundDiv(half-1,D) == 0);
  TEST_ASSERT(RoundDiv(half+1,D) == 1);
  TEST_ASSERT(RoundDiv(-half,D) == -RoundHalf);
  TEST_ASSERT(RoundDiv(-half+1,D) == 0);
  TEST_ASSERT(RoundDiv(-half-1,D) == -1);
  TEST_ASSERT(RoundDiv(d+half,D) == 1+RoundHalf);
  TEST_ASSERT(RoundDiv(d+half-1,D) == 1);
  TEST_ASSERT(RoundDiv(d+half+1,D) == 2);
  TEST_ASSERT(RoundDiv(-d-half,D) == -1-RoundHalf);
  TEST_ASSERT(RoundDiv(-d-half+1,D) == -1);
  TEST_ASSERT(RoundDiv(-d-half-1,D) == -2);

  TEST_ASSERT(RoundDiv(BigInt(1),d) == 0);
  TEST_ASSERT(RoundDiv(BigInt(-1),d) == 0);
  TEST_ASSERT(RoundDiv(HALF,d) == RoundHalf);
  TEST_ASSERT(RoundDiv(HALF-1,d) == 0);
  TEST_ASSERT(RoundDiv(HALF+1,d) == 1);
  TEST_ASSERT(RoundDiv(-HALF,d) == -RoundHalf);
  TEST_ASSERT(RoundDiv(-HALF+1,d) == 0);
  TEST_ASSERT(RoundDiv(-HALF-1,d) == -1);
  TEST_ASSERT(RoundDiv(d+HALF,d) == 1+RoundHalf);
  TEST_ASSERT(RoundDiv(d+HALF-1,d) == 1);
  TEST_ASSERT(RoundDiv(d+HALF+1,d) == 2);
  TEST_ASSERT(RoundDiv(-d-HALF,d) == -1-RoundHalf);
  TEST_ASSERT(RoundDiv(-d-HALF+1,d) == -1);
  TEST_ASSERT(RoundDiv(-d-HALF-1,d) == -2);

  TEST_ASSERT(RoundDiv(BigInt(1),D) == 0);
  TEST_ASSERT(RoundDiv(BigInt(-1),D) == 0);
  TEST_ASSERT(RoundDiv(HALF,D) == RoundHalf);
  TEST_ASSERT(RoundDiv(HALF-1,D) == 0);
  TEST_ASSERT(RoundDiv(HALF+1,D) == 1);
  TEST_ASSERT(RoundDiv(-HALF,D) == -RoundHalf);
  TEST_ASSERT(RoundDiv(-HALF+1,D) == 0);
  TEST_ASSERT(RoundDiv(-HALF-1,D) == -1);
  TEST_ASSERT(RoundDiv(D+HALF,D) == 1+RoundHalf);
  TEST_ASSERT(RoundDiv(D+HALF-1,D) == 1);
  TEST_ASSERT(RoundDiv(D+HALF+1,D) == 2);
  TEST_ASSERT(RoundDiv(-D-HALF,D) == -1-RoundHalf);
  TEST_ASSERT(RoundDiv(-D-HALF+1,D) == -1);
  TEST_ASSERT(RoundDiv(-D-HALF-1,D) == -2);

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
