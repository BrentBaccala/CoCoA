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
#include "CoCoA/FractionField.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/RingDistrMPolyClean.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/ideal.H"

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

  ring Q = RingQQ();
  SparsePolyRing Qx = NewPolyRing(Q, 1);
  const vector<RingElem>& x = indets(Qx);

  RingElem half(Q, 1);
  half /= 2;
  RingHom phi = CoeffEmbeddingHom(Qx);
  RingElem f(Qx);
  f = x[0] + phi(half);

  ideal I = ideal(f);
  cout << "f=" << f << endl;
  vector<RingElem> GBasis = TidyGens(I);
  cout << GBasis.front() << "\n";
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

  BuildInfo::PrintAll(cerr);
  return 1;
}

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/tests/test-bug2.C,v 1.8 2012/02/10 11:57:12 bigatti Exp $
// $Log: test-bug2.C,v $
// Revision 1.8  2012/02/10 11:57:12  bigatti
// -- changed RingZ.H, RingQ.H --> RingZZ.H, RingQQ.H
//
// Revision 1.7  2012/02/08 17:35:45  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.6  2010/12/17 16:06:25  abbott
// Ensured that all i/o is on standard C++ streams (instead of GlobalInput, etc)
//
// Revision 1.5  2010/10/22 12:42:03  abbott
// Added missing call to BuildInfo::PrintAll.
//
// Revision 1.4  2010/10/08 08:19:51  bigatti
// -- RingDistrMPoly.H --> RingDistrMPolyClean.H
//
// Revision 1.3  2009/07/30 15:41:25  bigatti
// -- now using the new nice constructors for ideals
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
