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


#include "CoCoA/BuildInfo.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/SparsePolyRing.H"
//#include "CoCoA/assert.H"  // for using TEST_ASSERT
#include "CoCoA/error.H"
#include "CoCoA/VectorOperations.H"
#include "CoCoA/symbol.H"


#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
#include <vector>
using std::vector;

using namespace CoCoA;

// Handy macro for making assertions.
#define TEST_ASSERT(cond)   do { if (!(cond)) CoCoA::AssertionFailed(#cond, __FILE__, __LINE__); } while(0)

void program()
{
  GlobalManager CoCoAFoundations;

  // Ripped off from ex-PolyRing3 :-)
  ring ZZ = RingZZ();
  SparsePolyRing ZZxy = NewPolyRing(ZZ, symbols("x","y"));
  const RingElem x = indet(ZZxy,0);
  const RingElem y = indet(ZZxy,1);

  RingElem f = 2*x*x*y - 4*y*y + 6*x*x + 36;

  cout << "In the following we have f = " << f << endl;

  // Accessing coeffs via SparsePolyIter:
  cout << "Using a SparsePolyIter we decompose f as follows:" << endl;
  for (SparsePolyIter it = BeginIter(f); !IsEnded(it); ++it)
  {
    cout << "coeff = " << coeff(it) << "  in ring " << owner(coeff(it)) << endl;
    cout << "PP    = " << PP(it)    << "  in " << owner(PP(it)) << endl;
    cout << endl;
  }

  cout << endl;
  cout << "-------------------------------------------------------" << endl;
  cout << endl;

  // Regard f as a poly in just "x" or just "y" we obtain:
  cout << "Coefficients with respect to certain indeterminates:" << endl << endl;
  cout << "Case (1) considering f as a *univariate* polynomial in..." << endl;
  cout << "...the indeterminate x the coeffs are " << CoeffVecWRT(f, x) << endl;
  cout << "...the indeterminate y the coeffs are " << CoeffVecWRT(f, y) << endl;
  cout << endl;

  cout << "Case (2) considering f as a sparse multivariate polynomial in..." << endl;
  cout << "...the indet x its structure is " << CoefficientsWRT(f, x) << endl;
  cout << "...the indet y its structure is " << CoefficientsWRT(f, y) << endl;
  vector<long> XandY; XandY.push_back(0); XandY.push_back(1);
  cout << "...the indets x & y its structure is " << CoefficientsWRT(f, XandY) << endl;
  cout << endl;
  cout << "-------------------------------------------------------" << endl;
  cout << "The content of a polynomial" << endl << endl;

  // Content of f
  RingElem ContF = content(f);
  cout << "The \"numerical\" content of f is " << ContF << "  -- an element of ring " << owner(ContF) << endl;
  cout << endl;
  RingElem ContWRTx = ContentWRT(f, x);
  cout << "Content WRT x is " << ContWRTx << "  -- element of " << owner(ContWRTx) << endl;

  RingElem ContWRTy = ContentWRT(f, y);
  cout << "Content WRT y is " << ContWRTy << "  -- element of " << owner(ContWRTy) << endl;

  RingElem ContWRTxy = ContentWRT(f, XandY);
  cout << "Content WRT x & y is " << ContWRTxy << "  -- element of " << owner(ContWRTxy) << endl;

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
