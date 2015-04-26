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
#include "CoCoA/error.H"
#include "CoCoA/random.H"

#include <iostream>
using std::cerr;
using std::endl;
#include <vector>
using std::vector;


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

  RandomSource RndSrc(0);

  const int NumBits = 11;
  const int N = 1<<NumBits;
  const int freq = 2000; // if you change this, then change 0.8 in TEST_ASSERT after loop
  vector<int> hist(N);

  // Make a histogram of the NumBits blocks of random bits.
  const int NumTrials = N*freq;
  for (int i=0; i < NumTrials; ++i)
  {
    int sample = 0;
    for (int j=0; j < NumBits; ++j)
    {
      sample <<= 1;
      if (RandomBool(RndSrc))
        ++sample;
    }
    ++hist[sample];
  }

  // Now find highest and lowest frequencies.
  int max = 0;
  int min = freq;
  for (int i=0; i < N; ++i)
  {
    if (hist[i] > max)  max = hist[i];
    else if (hist[i] < min)  min = hist[i];
  }

  // Check that highest and lowest frequencies are not too different.
  TEST_ASSERT(min > 0.8*max);
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
