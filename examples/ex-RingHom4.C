// Copyright (c) 2005 Anna Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"
using namespace CoCoA;
using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This program shows how we define a ring homomorphism       \n"
  "to evaluate polynomials.                                   \n";
//----------------------------------------------------------------------

void program()
{
  GlobalManager CoCoAFoundations;

  cout << ShortDescription << endl;

  // Create some coefficient rings
  ring Fp = NewZZmod(32003);

  PolyRing Fpx = NewPolyRing(Fp, 4); // Fpx is Z/(32003)[x[0..3]]
  cout << "  --- PolyRing is " << Fpx << endl;

  const vector<RingElem>& x = indets(Fpx);
  RingElem f = 4*x[1]*x[1] + x[0]*x[2] - 2*x[0];
  cout << "  --- f = " << f << endl;

  RingElem TwoThirds(Fp);
  TwoThirds = RingElem(Fp,2)/3; // C++ reads 2/3 as an integer division! 
  vector<RingElem> images;
  images.push_back(one(Fp));          // x[0] |-> 1
  images.push_back(zero(Fp));         // x[1] |-> 0
  images.push_back(TwoThirds);        // x[2] |-> 2/3
  images.push_back(RingElem(Fp,100)); // x[3] |-> 100

  RingHom eval = PolyAlgebraHom(Fpx, Fp, images);
  // same as PolyRingHom(Fpx, Fp, IdentityHom(Fp), images);
  cout << x[0] << " |-> " << eval(x[0]) << endl;
  cout << x[1] << " |-> " << eval(x[1]) << endl;
  cout << x[2] << " |-> " << eval(x[2]) << endl;
  cout << x[3] << " |-> " << eval(x[3]) << endl;
  cout << "f    |-> " << eval(f) << endl;
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
    cerr << "***ERROR***  UNCAUGHT CoCoAError";
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-RingHom4.C,v 1.4 2012/07/04 12:18:02 abbott Exp $
// $Log: ex-RingHom4.C,v $
// Revision 1.4  2012/07/04 12:18:02  abbott
// Improved printing: actually use the indet names instead of hardwiring the names in a string.
//
// Revision 1.3  2012/02/08 17:45:51  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.2  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.9  2007/03/08 14:38:07  cocoa
// Added new range function in symbol.H, and tidied many calls to PolyRing
// pseudo ctors (as a consequence).
//
// Revision 1.8  2007/03/07 14:52:56  bigatti
// -- minor: text and comments
//
// Revision 1.7  2007/03/07 14:32:53  bigatti
// -- now uses PolyAlgebraHom
// -- minor cleaning
//
// Revision 1.6  2007/03/05 21:06:07  cocoa
// New names for homomorphism pseudo-ctors: removed the "New" prefix.
//
// Revision 1.5  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.4  2007/02/26 17:39:46  bigatti
// -- getting ready for unique ring Z: using NewZmod(N), NewRingQ()
//
// Revision 1.3  2007/02/12 16:15:37  bigatti
// -- added strings ShortDescription and LongDescription for indexing
//
// Revision 1.2  2007/02/10 18:44:04  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.2  2006/04/27 16:12:07  cocoa
// -- changed NewIdentityRingHom --> NewIdentityHom
//
// Revision 1.1.1.1  2005/10/17 10:46:53  cocoa
// Imported files
//
// Revision 1.2  2005/10/14 15:25:07  cocoa
// Major tidying and cleaning to small prime finite fields.
// Several consequential changes.  Improved their documentation.
//
// Added Makefile and script to include/CoCoA/ directory to
// keep library.H up to date.
//
// Revision 1.1  2005/09/30 13:01:34  cocoa
// -- first import
//
// Revision 1.3  2005/09/22 18:04:17  cocoa
// It compiles; the tests run OK.  The examples compile.
// No documentation -- the mindless eurocrats have rendered
// me mindless too.
//
// Revision 1.2  2005/06/22 14:47:56  cocoa
// PPMonoids and PPMonoidElems updated to mirror the structure
// used for rings and RingElems.  Many consequential changes.
//
// Revision 1.1  2005/05/04 16:35:06  cocoa
// -- new examples
//
