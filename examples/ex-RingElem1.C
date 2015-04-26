// Copyright (c) 2005  John Abbott
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"
using namespace CoCoA;
using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Example showing operations on RingElem for a ring or a PolyRing.\n";

const string LongDescription =
  "This is a long list of function calls from different rings.  \n";
//----------------------------------------------------------------------

// Includes from the standard C++ library
// #include <iostream> // using std::endl;


//-- TestRing ----------------------------------------------------------
// testing some ring operations over elements of a given ring

void TestRing(ring R, RingElem a, RingElem b)
{
  RingElem one(R, 1);
  cout << "  initial ring elements are:" << endl;
  cout << "    a = " << a << endl;
  cout << "    b = " << b << endl;
  cout << endl;

  cout << "a == b gives " << (a == b) << endl;
  cout << "a != b gives " << (a != b) << endl;
  cout << "a * b  gives " << a*b << endl;
  cout << "-a     gives " << -a << endl;
  cout << "a += b gives a = " << (a += b) << endl;
  cout << endl;

  cout << "IsZero(one)     gives " << IsZero(one) << endl;
  cout << "IsOne(one)      gives " << IsOne(one) << endl;
  cout << "IsMinusOne(one) gives " << IsMinusOne(one) << endl;
  cout << endl;
  cout << "power(a, 10)  gives "  << power(a, 10) << endl;
  cout << endl;

  cout << "IsDivisible(a, b) gives " << IsDivisible(a, b) << endl;
  if (!IsDivisible(a, b)) 
    cout << "  so we CANNOT compute b / a" << endl;
  else 
    cout << "  b / a gives " << b/a << endl;
  cout << endl;

  cout << "IsFractionField(R) gives " << IsFractionField(R) << endl;
  if (!IsFractionField(R)) 
    cout << "  so we CANNOT compute num(a), den(a)" << endl;
  else 
  {
    cout << "  num(a)   gives " << num(a) << endl;
    cout << "  den(a)   gives " << den(a) << endl;
  }
  cout << endl;

  cout << "IsTrueGCDDomain(R) gives " << IsTrueGCDDomain(R) << endl;
  if (!IsTrueGCDDomain(R)) 
    cout << "  so we CANNOT compute gcd(a,b)" << endl;
  else 
    cout << "  gcd(a, b)   gives " << gcd(a, b) << endl;
  cout << endl;

  cout << "IsPolyRing(R) gives " << IsPolyRing(R) << endl;
  if (!IsPolyRing(R)) 
    cout << "  so we CANNOT compute deg(a) or StdDeg(a)" << endl;
  else 
  {
    cout << "  NB deg(a) and StdDeg(a) are synonyms;" << endl
         << "  StdDeg is a more precise name but is also more cumbersome" << endl;
    cout << "  deg(a)    gives " << deg(a) << endl;
    cout << "  StdDeg(a) gives " << StdDeg(a) << endl;
  }
  cout << endl;

  cout << "IsSparsePolyRing(R) gives " << IsSparsePolyRing(R) << endl;
  if (!IsSparsePolyRing(R)) 
    cout << "  so we CANNOT compute LPP(a), wdeg(a)" << endl;
  else 
  {
    cout << "  LPP(a)  gives " << LPP(a) << endl;
    cout << "  wdeg(a) gives " << wdeg(a) << endl;
    cout << "  NB wdeg(a) agrees with deg(a) only if the ring is standard graded." << endl;
    if (GradingDim(R)>0)
      cout << "  LF(a) gives " << LF(a) << endl;
  }
  cout << endl;

}


//-- main --------------------------------------------------------------
// we run TestRing on some rings

void program()
{
  GlobalManager CoCoAFoundations;

  cout << ShortDescription << endl;

  //------------------------------------------------------------
  // Z
  ring ZZ = RingZZ();
  {
    RingElem a(ZZ), b(ZZ);
    
    a = 5;
    b = 7;
    cout << "----------------------------------------" << endl;
    cout << "  ring is ZZ" << endl;
    TestRing(ZZ, a, b);
  }
  
  //------------------------------------------------------------
  // Q
  ring QQ = RingQQ();  
  {
    RingElem a(QQ), b(QQ);
    
    a = 5; a /= 2;
    b = 7; b /= 3;
    cout << "----------------------------------------" << endl;
    cout << "  ring is QQ" << endl;
    TestRing(QQ, a, b);
  }

  //------------------------------------------------------------
  // QQ[x,y]

  PolyRing P = NewPolyRing(QQ, symbols("x","y"));
  RingElem f = ReadExpr(P, "15*y + y^3");
  RingElem g = ReadExpr(P, "4*y - 3*y^7 + x^7");
  cout << "----------------------------------------" << endl;
  cout << "  ring is QQ[x,y]" << endl;
  TestRing(P, f, g);
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-RingElem1.C,v 1.8 2014/07/08 12:47:26 abbott Exp $
// $Log: ex-RingElem1.C,v $
// Revision 1.8  2014/07/08 12:47:26  abbott
// Summary: Removed AsPolyRing, AsSparsePolyRing, AsQuotientRing
// Author: JAA
//
// Revision 1.7  2014/03/21 18:16:00  bigatti
// -- now using ReadExpr
//
// Revision 1.6  2012/10/05 10:21:39  bigatti
// -- added LF (leading form)
//
// Revision 1.5  2012/05/22 10:02:38  abbott
// Removed IsGCDDomain; substituted by IsTrueGCDDomain.
// Added IsFractionFieldOfGCDDomain.
//
// Revision 1.4  2012/05/20 09:43:54  abbott
// Corrected layout.
//
// Revision 1.3  2012/02/08 17:43:14  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.2  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.6  2007/03/08 16:55:06  cocoa
// Changed name of "range" function to "SymbolRange".
//
// Revision 1.5  2007/03/08 14:38:07  cocoa
// Added new range function in symbol.H, and tidied many calls to PolyRing
// pseudo ctors (as a consequence).
//
// Revision 1.4  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.3  2007/02/12 15:45:15  bigatti
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
// Revision 1.3  2006/03/01 14:26:41  cocoa
// -- removed some useless "include"
//
// Revision 1.2  2006/02/13 12:08:04  cocoa
// Fixed a problem with some missing assignment ops for certain PPMonoidElems.
// Fixed a bug in RingDistrMPoly::myIndetPower.
//
// Revision 1.1.1.1  2005/10/17 10:46:53  cocoa
// Imported files
//
// Revision 1.3  2005/09/22 18:04:17  cocoa
// It compiles; the tests run OK.  The examples compile.
// No documentation -- the mindless eurocrats have rendered
// me mindless too.
//
// Revision 1.2  2005/07/19 15:30:20  cocoa
// A first attempt at iterators over sparse polynomials.
// Main additions are to SparsePolyRing, DistrMPoly*.
// Some consequential changes to PPMonoid*.
//
// Revision 1.1.1.1  2005/05/03 15:47:30  cocoa
// Imported files
//
// Revision 1.3  2005/04/27 16:14:56  cocoa
// Cleaned up example programs -- added "free use" permit.
// Changed a couple of ErrorInfo object names, and added
// ERR::NotTrueGCDDomain.
//
// Revision 1.2  2005/04/21 15:12:19  cocoa
// Revised NewPolyRing as Dag Arneson suggested (perhaps just an interim
// measure).
// Brought example programs up to date (new name for CoCoA error
// information objects).
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.2  2004/12/09 15:08:42  cocoa
// -- added log info
//
