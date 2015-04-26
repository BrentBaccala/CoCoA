// Copyright (c) 2008  Anna Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"
using namespace CoCoA;
using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Example of use of power products in different PPMonoids.   \n"
  "Program exhibiting timings of the different implementations.\n";

const string LongDescription =
  "The implementations of PPMonoids are optimized for different uses:  \n"
  "PPMonoidEv:   stores the Exponent vector                            \n"
  "              good for accessing the exponents, slow for ordering   \n"
  "PPMonoidOv:   stores the Order vector                               \n"
  "              good for ordering, slow for accessing the exponents   \n"
  "PPMonoidEvOv: stores the Exponent vector and the Order vector       \n"
  "              good for accessing the exponents and for ordering     \n"
  "              uses more memory and takes more time to assign        \n"
  "PPMonoidBigEv: stores the Exponent vector as BigInt's               \n"
  "              necessary if you use big exponents (>2^10)            \n"
  "              obviously slow   \n";

//----------------------------------------------------------------------

// Includes from the standard C++ library
// #include <iostream> // using std::endl;


void TryPPMonoid(const PPMonoid PPM)
{
  //  cout << PPM << endl;
  const vector<PPMonoidElem>& x = indets(PPM);

  // at least 20 indets!
  PPMonoidElem t1 = product(x);
  PPMonoidElem t2 = x[0] * x[1] * x[20];

  double t0;  // for CpuTime

  const long NumIters = 70000;

  // if your optimizer is too clever it might realise these operations
  // are "useless" (return value not used and no side effect)
  // so you will get constant timings
  t0 = CpuTime();
  for (long i=1; i < NumIters; ++i)    t2 < t1;
  cout << " <    -->          " << CpuTime()-t0 << " s" << endl;
  t0 = CpuTime();
  for (long i=1; i < NumIters; ++i)    t1 * t2;
  cout << " *    -->          " << CpuTime()-t0 << " s" << endl;
  t0 = CpuTime();
  for (long i=1; i < NumIters; ++i)    gcd(t1, t2);
  cout << " gcd  -->          " << CpuTime()-t0 << " s" << endl;  
  for (long i=1; i < NumIters; ++i)    IsDivisible(t1, t2);
  cout << " IsDivisible  -->  " << CpuTime()-t0 << " s" << endl;  
}


void TryOrd(const PPOrdering ord)
{
  const vector<symbol> x = SymbolRange("x", 0, NumIndets(ord)-1);
  cout << "-- =========================================== --" << endl;
  cout << "   " << ord << endl;
  cout << "-- =========================================== --" << endl;
  cout << "---------   PPMonoidBigEv  ----------" << endl;
  TryPPMonoid(NewPPMonoidEv(x, ord, BigExps));
  cout << "----------   PPMonoidEv  -----------" << endl;
  TryPPMonoid(NewPPMonoidEv(x, ord));
  cout << "----------   PPMonoidOv  -----------" << endl;
  TryPPMonoid(NewPPMonoidOv(x, ord));
  cout << "---------   PPMonoidEvOv  ----------" << endl;
  TryPPMonoid(NewPPMonoidEvOv(x, ord));
}


void program()
{
  GlobalManager CoCoAFoundations;

  cout << ShortDescription << endl;
  const int n = 40;
  PPOrdering lex = NewLexOrdering(n);
  PPOrdering DegRevLex = NewStdDegRevLexOrdering(n);
  vector<long> y(1,1);  // y = [1]
  PPOrdering MatOrd = NewMatrixOrdering(n, 1, ElimMat(n, y));

  TryOrd(lex);
  TryOrd(DegRevLex);
  TryOrd(MatOrd);
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-PPMonoidElem2.C,v 1.10 2013/03/15 11:04:56 abbott Exp $
// $Log: ex-PPMonoidElem2.C,v $
// Revision 1.10  2013/03/15 11:04:56  abbott
// Updated to new interface for NewPPMonoidEv pseudo-ctor.
//
// Revision 1.9  2013/02/14 17:49:14  bigatti
// -- renamed NewMatElim --> ElimMat
//
// Revision 1.8  2012/05/04 20:01:19  abbott
// Minor modifications to reduce execution time (usu. by reducing number of iterations).
//
// Revision 1.7  2011/08/23 12:04:04  bigatti
// -- updated after renaming ZZ --> BigInt
//
// Revision 1.6  2011/03/08 18:04:16  bigatti
// -- changed size_t into long
//
// Revision 1.5  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.4  2009/09/22 15:27:07  bigatti
// -- NewMatrixElim --> NewMatElim
//
// Revision 1.3  2009/05/14 09:11:39  abbott
// Added a non-essential const.
//
// Revision 1.2  2008/11/19 11:12:07  bigatti
// -- added PPMonoidEvZZ and more explanations
//
// Revision 1.1  2008/11/19 09:20:50  bigatti
// -- first import
//
// Revision 1.2  2007/12/04 14:27:07  bigatti
// -- changed "log(pp, i)" into "exponent(pp, i)"
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.8  2007/03/08 22:26:27  cocoa
// Removed try..catch constructs from some "simple" examples.
//
// Revision 1.7  2007/03/08 17:43:11  cocoa
// Swapped order of args to the NewPPMonoid pseudo ctors.
//
// Revision 1.6  2007/03/08 16:55:06  cocoa
// Changed name of "range" function to "SymbolRange".
//
// Revision 1.5  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.4  2007/02/28 14:00:13  bigatti
// -- minor: just a comment
//
// Revision 1.3  2007/02/12 15:31:57  bigatti
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
// Revision 1.1.1.1  2005/10/17 10:46:53  cocoa
// Imported files
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
