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
#include "CoCoA/GlobalManager.H"
#include "CoCoA/PPMonoidEv.H"
#include "CoCoA/RingDistrMPolyClean.H"
#include "CoCoA/RingTwinFloat.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/ring.H"
#include "CoCoA/symbol.H"

using namespace CoCoA;

#include <iostream>
using std::cerr;
using std::endl;
#include <vector>
using std::vector;

#define TEST_ASSERT(cond)   do { if (!(cond)) CoCoA::AssertionFailed(#cond, __FILE__, __LINE__); } while(0)


void solve(long BuffBits)
{
  const ring RR = NewRingTwinFloat(20, BuffBits, 32);

  const PolyRing P = NewPolyRing(RR, SymbolRange("x", 0, 2), NewStdDegLexOrdering(3));

  const RingElem x = indet(P, 0);
  const RingElem y = indet(P, 1);
  const RingElem z = indet(P, 2);
  const RingElem f = power(x,41)-power(z,40)*(x-z);
//  RingElem g = power(x,5)*y-2*power(x,3)-3*x+1;
  const RingElem g = power(x,6)- 9*power(x,5)*z+x*power(z,5)+11*power(z,6)-power(y,6);

  vector<RingElem> gens;
  gens.push_back(f);
  gens.push_back(g);
  const ideal I = ideal(gens);

//  vector<RingElem> GI = CoCoA::gens(I);
  const vector<RingElem> GBasis = TidyGens(I); // might throw an exception
  TEST_ASSERT(BuffBits >= 143);
  TEST_ASSERT(GBasis.size() == 64);
  TEST_ASSERT(GBasis[0] == g);
}


void program()
{
  GlobalManager CoCoAFoundations;

  for (long BuffBits = 32; BuffBits < 160; ++BuffBits)
  {
    try
    {
      solve(BuffBits);
    }
    catch (const RingTwinFloat::InsufficientPrecision&)
    {
      TEST_ASSERT(BuffBits < 143);
    }
  }
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
