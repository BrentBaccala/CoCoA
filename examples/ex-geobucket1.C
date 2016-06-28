// Copyright (c) 2006  John Abbott
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"
using namespace CoCoA;
using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This program shows how to use geobuckets for long sums.          \n"
  "It compares timings between normal sum \"+=\" and geobucket sum. \n";

const string LongDescription =
  "We simulate a long sum:\n"
  "we add the summands of a long polynomial one at a time,               \n"
  "and do it twice to consider also the arithmetics on the coefficients. \n"
  "Geobuckets are to be used as a temporary value;                       \n"
  "the final result is then copied into a RingElem.";

//----------------------------------------------------------------------


void CompareSums(SparsePolyRing P)
{
  RingElem f = power(sum(indets(P)) + 1, 10);
  PPMonoidElem x(indets(PPM(P))[0]);
  RingElem f_copy(P);
  RingElem f_gbk(P); // just the copy of gbk final value
  geobucket gbk(P);
  cout << "----------------------------------------------------" << endl;
  cout << "---- P is          " << P << endl;
  cout << "---- NumTerms(f) = " << NumTerms(f) << endl;
  cout << "----------------------------------------------------" << endl;

  double t0 = CpuTime();
  for (long i=0; i<2; ++i)
    for (SparsePolyIter it=BeginIter(f); !IsEnded(it); ++it)
      f_copy += monomial(P, coeff(it), PP(it)*x);
  cout << "sum using \"+=\"      time = " << CpuTime()-t0 << endl;  

  t0 = CpuTime();
  for (long i=0; i<2; ++i)
    for (SparsePolyIter it=BeginIter(f); !IsEnded(it); ++it)
      gbk.myAddMul(monomial(P, coeff(it), PP(it)), one(P), 1);
  AddClear(f_gbk, gbk); // copies gbk to f_gbk and clears gbk
  cout << "sum using geobucket time = " << CpuTime()-t0 << endl;

  cout <<  "check:  2*f == f_copy is " << (2*f == f_copy);
  cout <<  ";  2*f == f_gbk is " << (2*f == f_gbk) << endl;
}


void program()
{
  GlobalManager CoCoAFoundations;

  cout << ShortDescription << endl;
  cout << LongDescription << endl;
  cout << boolalpha; // prints true/false for bool

  CompareSums(NewPolyRing(RingQQ(), 3));
  CompareSums(NewPolyRing(RingQQ(), 4));
  CompareSums(NewPolyRing(RingQQ(), 6));
  CompareSums(NewPolyRing_DMPI(RingQQ(), 6));
  CompareSums(NewPolyRing_DMPII(NewRingFp(101), 6));
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-geobucket1.C,v 1.4 2014/07/08 12:21:44 abbott Exp $
// $Log: ex-geobucket1.C,v $
// Revision 1.4  2014/07/08 12:21:44  abbott
// Summary: Removed AsSparsePolyRing
// Author: JAA
//
// Revision 1.3  2014/01/28 13:08:10  bigatti
// -- minor improvements
//
// Revision 1.2  2012/10/05 15:32:52  bigatti
// -- minor changes after improving +=
//
// Revision 1.1  2012/10/05 06:45:29  bigatti
// -- first import
//
// Revision 1.5  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.4  2008/10/07 12:12:54  abbott
// Removed useless commented out #include.
//
// Revision 1.3  2007/05/31 16:06:16  bigatti
// -- removed previous unwanted checked-in version
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.9  2007/03/07 11:51:40  bigatti
// -- improved test alignment
//
// Revision 1.8  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.7  2007/03/02 17:46:40  bigatti
// -- unique RingZ and RingQ
// -- requires foundations.H ;  foundations blah;  (thik of a better name)
//
// Revision 1.6  2007/03/02 10:47:53  cocoa
// First stage of RingZ modifications -- tests do not compile currently, Anna will fix this.
//
// Revision 1.5  2007/03/01 13:52:59  bigatti
// -- minor: fixed typo
//
// Revision 1.4  2007/02/28 15:15:56  bigatti
// -- minor: removed quotes in description
//
// Revision 1.3  2007/02/12 16:27:43  bigatti
// -- added strings ShortDescription and LongDescription for indexing
//
// Revision 1.2  2007/02/10 18:44:03  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.1  2006/03/12 21:28:34  cocoa
// Major check in after many changes
//
