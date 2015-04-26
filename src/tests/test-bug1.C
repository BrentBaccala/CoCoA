//   Copyright (c)  2005  John Abbott

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


// This bug originally reported by Dag Arneson 2005/04/01
// It used to provoke SEGV.

#include "CoCoA/BuildInfo.H"
#include "CoCoA/FractionField.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/PPMonoidEv.H"
#include "CoCoA/PPOrdering.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingDistrMPolyClean.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/symbol.H"

using namespace CoCoA;

// Includes from the standard C++ library
#include <iostream>
using std::cerr;
using std::endl;


void program()
{
  GlobalManager CoCoAFoundations;

  //To make the polynomial ring we have to make a PP monoid
  PolyRing Qp = NewPolyRing(RingQQ(), SymbolRange("y", 0, 2)); // Q[y1..yn]

  RingElem foo(Qp);
  foo = 5*indet(Qp,0)*indet(Qp,1)+2*indet(Qp,0)*indet(Qp,2);
  ideal I(foo);
  QuotientRing R = NewQuotientRing(Qp,I);
}


//----------------------------------------------------------------------
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
