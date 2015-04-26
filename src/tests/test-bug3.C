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


// Don't remember who originally reported this.
// The GBasis used to "ignore" the denominator in f.

#include "CoCoA/BuildInfo.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/IntOperations.H"
#include "CoCoA/NumTheory.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingDistrMPolyClean.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/BigInt.H"
#include "CoCoA/symbol.H"

using namespace CoCoA;

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <vector>
using std::vector;

void program()
{
  GlobalManager CoCoAFoundations;

  const int p = 13; // a small prime number
  PolyRing P1 = NewPolyRing(NewZZmod(p), symbols("x"));
  const RingElem x = indet(P1,0);
  cout << x-1 << endl;

  const int n = 20; // a small composite number
  PolyRing P2 = NewPolyRing(NewZZmod(n), symbols("y"));
  const RingElem y = indet(P2,0);
  cout << y-1 << endl;

  const BigInt P = NextProbPrime(power(10,20)); // a large prime number
  PolyRing P3 = NewPolyRing(NewZZmod(P), symbols("z"));
  const RingElem z = indet(P3, 0);
  cout << z-1 << endl;

  const BigInt N = power(10,20); // a large composite number
  PolyRing P4 = NewPolyRing(NewZZmod(N), symbols("t"));
  const RingElem t = indet(P4, 0);
  cout << t-1 << endl;
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
    cerr << "***ERROR***  UNCAUGHT CoCoA error";
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
  return 1;
}

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/tests/test-bug3.C,v 1.6 2012/05/28 09:18:20 abbott Exp $
// $Log: test-bug3.C,v $
// Revision 1.6  2012/05/28 09:18:20  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.5  2012/02/08 17:36:48  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.4  2011/08/23 06:40:31  bigatti
// -- fixed with new name "BigInt" for old "ZZ"
//
// Revision 1.3  2010/12/17 16:06:25  abbott
// Ensured that all i/o is on standard C++ streams (instead of GlobalInput, etc)
//
// Revision 1.2  2010/10/08 08:19:51  bigatti
// -- RingDistrMPoly.H --> RingDistrMPolyClean.H
//
// Revision 1.1  2009/12/29 22:45:17  abbott
// Added two new tests.
//
// Revision 1.2  2007/10/30 17:14:05  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:12  abbott
// Imported files
//
// Revision 1.6  2007/03/03 14:13:21  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.5  2007/03/02 17:46:40  bigatti
// -- unique RingZ and RingQ
// -- requires foundations.H ;  foundations blah;  (thik of a better name)
//
// Revision 1.4  2007/02/26 17:11:58  bigatti
// -- getting ready for unique ring Z: using NewZmod(N), NewRingQ()
//
// Revision 1.3  2007/02/10 18:44:03  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.2  2006/12/07 17:25:32  cocoa
// -- minimal set of #include's instead of library.H
//
// Revision 1.1  2006/10/06 14:04:57  cocoa
// The new assert header and implementation files.
// A new test.
//
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.1  2006/03/12 21:28:34  cocoa
// Major check in after many changes
//
