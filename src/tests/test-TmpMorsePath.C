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
#include "CoCoA/TmpMorsePaths.H"

#include "CoCoA/error.H"

using namespace CoCoA;
using namespace std;

// Handy macro for making assertions.
#define TEST_ASSERT(cond)   do { if (!(cond)) CoCoA::AssertionFailed(#cond, __FILE__, __LINE__); } while(0)

void test_myAddPath()
{
  ring Q = RingQQ();
  SparsePolyRing PolyRing = NewPolyRing(Q, SymbolRange("x",0,4), StdDegLex);
  const vector<RingElem> x = indets(PolyRing);
  vector<RingElem> PolInput;
  PolInput.push_back(x[0]);
  PolInput.push_back(x[1]);
  PolInput.push_back(x[2]);
  PolInput.push_back(x[3]);
  PolInput.push_back(x[4]);
  JBMill mill = ExtendedJanetBasis(PolInput);
  std::vector<std::pair<RingElem, DynamicBitset> > TestBasis;
  std::map<MorseElement, MorsePaths> TestResolution;
  std::vector<RingElem> JBBasis(JBReturnJB(mill));
  for (std::vector<RingElem>::iterator i = JBBasis.begin(); i != JBBasis.end(); ++i)
  {
    TestBasis.push_back(std::pair<RingElem, DynamicBitset>(*i, DynamicBitset(5)));
  }
  for (std::vector<std::pair<RingElem, DynamicBitset> >::iterator i = TestBasis.begin(); i != TestBasis.end(); ++i)
  {
    TestResolution.insert(std::pair<MorseElement, MorsePaths>(MorseElement(DynamicBitset(5), i), MorsePaths())); 
  }
  MorsePaths test;
  ConstResIter iter(TestResolution.begin());
  test.myAddPath(iter, x[0]);
  TEST_ASSERT(x[0] == test.myGetPath(iter));
  ++iter;
  test.myAddPath(iter, x[1]);

  TEST_ASSERT(x[1] == test.myGetPath(iter));

  --iter;
  test.myAddPath(iter, -x[0]);
  PathMap paths(test.myGetPaths());
  TEST_ASSERT(paths.size() == 1);
  ++iter;
  TEST_ASSERT(paths[iter] == x[1]);
}

void program()
{
  GlobalManager CoCoAFoundations;
  test_myAddPath();

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
