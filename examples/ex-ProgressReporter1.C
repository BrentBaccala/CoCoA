// Copyright (c) 2014  John Abbott
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"
using namespace CoCoA;
using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This is a simple example showing how to use a ProgressReporter. \n";

const string LongDescription =
  "An example of how to use a ProgressReporter to print occasional \n"
  "updates during a long iterative computation.                    \n";
//----------------------------------------------------------------------

bool IsSquare(long n)
{
  BigInt junk;
  return IsExactIroot(junk, n, 2);
}

void program()
{
  GlobalManager CoCoAFoundations;

  cout << ShortDescription << endl;
  cout << boolalpha; // so that bools print out as true/false

  cout << "-------------------------------------------------------" << endl;
  cout << "First loop:" << endl;

  ProgressReporter ProgRep(2);
  for (long n=1; n < 100000; ++n)
  {
    ProgRep(); // report simply says says how many times ProgRep() was called
    if (n%4 == 0 || n%4 == 3) continue;
    long NumReprs = 0;
    for (long j=1; 2*j*j <= n; ++j)
      if (IsSquare(n-j*j)) ++NumReprs;
    if (NumReprs > 8)
      cout << n << " has " << NumReprs << " representations of the form A^2+B^2" << endl;
  }

  cout << endl;
  cout << "-------------------------------------------------------" << endl;
  cout << "Second loop:" << endl;

  ProgressReporter report(2); // print reports every 5 seconds
  for (long p=101; ; p = NextPrime(p))
  {
    report(p); // report will say which prime we have reached
    if (PrimitiveRoot(p) <= 50) continue;
    cout << "Prime " << p << " has least primitive root " << PrimitiveRoot(p) << endl;
    break;
  }
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-ProgressReporter1.C,v 1.3 2014/07/28 14:43:55 abbott Exp $
// $Log: ex-ProgressReporter1.C,v $
// Revision 1.3  2014/07/28 14:43:55  abbott
// Summary: Reduced a parameter to make the example faster (was too slow on my netbook)
// Author: JAA
//
// Revision 1.2  2014/04/22 14:52:23  abbott
// Summary: Made slightly clearer
// Author: JAA
//
// Revision 1.1  2014/04/22 14:09:12  abbott
// Summary: Example for ProgressReporter
// Author: JAA
//
// Revision 1.7  2013/05/28 07:07:04  bigatti
// -- added "cout << boolalpha": useful for testing
//
// Revision 1.6  2012/11/30 14:04:55  abbott
// Increased visibility of comment saying "put your code here".
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
