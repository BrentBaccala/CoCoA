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


#include "CoCoA/BuildInfo.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/random.H"
#include "CoCoA/assert.H"  // for using TEST_ASSERT
#include "CoCoA/error.H"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <vector>
using std::vector;

using namespace CoCoA;

// Handy macro for making assertions.
#define TEST_ASSERT(cond)   do { if (!(cond)) CoCoA::AssertionFailed(#cond, __FILE__, __LINE__); } while(0)

void program()
{
  GlobalManager CoCoAFoundations;

  const int N = 10;
  RandomSeqLong RndLong(-N,N);
  cout << "RndLong=" << RndLong << endl;

  const int ExpectedFreq = 256; // should be larger than 100ish
  vector<int> histo(2*N+1);
  for (int i=0; i < ExpectedFreq*(2*N+1); ++i)
  {
    const long rnd = NextValue(RndLong);
    TEST_ASSERT(rnd >= -N && rnd <= N);
    ++histo[N+rnd];
  }

  // Check histogram: each entry is > 0, and greatest and least freqs are within a factor of 2
  int MostFreq = histo[0];
  int LeastFreq = histo[0];
  for (int i=0; i < 2*N+1; ++i)
  {
    TEST_ASSERT(histo[i] != 0);
    if (histo[i] > MostFreq) { MostFreq = histo[i]; continue; }
    if (histo[i] < LeastFreq) { LeastFreq = histo[i]; }
  }
  const double ratio = double(MostFreq)/LeastFreq;
  TEST_ASSERT(ratio < 2.0);

  cout << "RndLong=" << RndLong << endl;
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
