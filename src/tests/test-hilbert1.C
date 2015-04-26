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


#include "CoCoA/BuildInfo.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/TmpHilbert.H"
#include "CoCoA/time.H"

using namespace CoCoA;

#include <algorithm>
using std::min;
#include <iostream>
using std::cerr;
using std::endl;
#include <vector>
using std::vector;

//----------------------------------------------------------------------
// First test for Hilbert.
//----------------------------------------------------------------------

#define TEST_ASSERT(cond)   do { if (!(cond)) CoCoA::AssertionFailed(#cond, __FILE__, __LINE__); } while(0)

// ---- chess tests -----------------------------------------------------

ConstRefRingElem CsbSquareIndet(SparsePolyRing P, long l, long sq1, long sq2)
{
  TEST_ASSERT( l*l <= NumIndets(P) );
  TEST_ASSERT( sq1 <= l && sq2 <= l );
  return indet(P, (sq1-1)*l + (sq2-1));
}


ideal NewQueenMovesFrom(SparsePolyRing P, long Csb, long sq1, long sq2)
{
  ConstRefRingElem x = CsbSquareIndet(P, Csb, sq1, sq2);
  vector<RingElem> g;
  for ( long i=sq2+1 ; i<=Csb ; ++i )
    g.push_back(x * CsbSquareIndet(P, Csb, sq1, i));
  for ( long i=sq1+1 ; i<=Csb ; ++i )
    g.push_back(x * CsbSquareIndet(P, Csb, i, sq2));
  for ( long i=min(Csb-sq1,Csb-sq2) ; i>0 ; --i )
    g.push_back(x * CsbSquareIndet(P, Csb, sq1+i, sq2+i));
  for ( long i=min(Csb-sq1, sq2-1) ; i>0 ; --i )
    g.push_back(x * CsbSquareIndet(P, Csb, sq1+i, sq2-i));
  return ideal(P, g); // ideal(P, g) because g might be empty
}


ideal NewQueenIdeal(SparsePolyRing P, long Csb)
{
  ideal I = ideal(zero(P));
  for ( long sq1=1 ; sq1<=Csb ; ++sq1 )
    for ( long sq2=1 ; sq2<=Csb ; ++sq2 )
      I += NewQueenMovesFrom(P, Csb, sq1, sq2);
  return I;
}


// ----------------------------------------------------------------------

void program()
{
  GlobalManager CoCoAFoundations;

///  double t0 = CpuTime();

  SparsePolyRing P = NewPolyRing(RingQQ(), 4);
  const vector<RingElem>& x = indets(P);

  ideal I = ideal(zero(P));
  RingElem HS = HilbertNumQuot(I);
  TEST_ASSERT(HS == 1);
  
  I = ideal(x[1], x[2], x[3]);
  HS = HilbertNumQuot(I);

  const PolyRing HPRing = owner(HS);
  const RingElem lambda = indet(HPRing, 0);
  //  SparsePolyRing HPRingS = NewPolyRing(RingZZ(), symbols("lambda"));
  PolyRing QQt = RingQQt(1);
  RingElem t = indet(QQt, 0);
  RingHom UniToPoly=PolyAlgebraHom(HPRing, QQt, indets(QQt));

  TEST_ASSERT(HS == -power(lambda,3) +3*power(lambda,2) -3*lambda +1);
  TEST_ASSERT(HilbertNumQuot_C(I) == -power(lambda,3) +3*power(lambda,2) -3*lambda +1);
  TEST_ASSERT(MGHilbertNumQuot(I) == -power(t,3) +3*power(t,2) -3*t +1);
  
  SparsePolyRing CsbRing = NewPolyRing(RingQQ(), 9*9);

  ideal Q3 = NewQueenIdeal(CsbRing, 3);
  RingElem HN3 = HilbertNumQuot_C(Q3);
  RingElem HN3_CPP = HilbertNumQuot(Q3);
  RingElem MGHN3_CPP = MGHilbertNumQuot(Q3);
  TEST_ASSERT(HN3 == HN3_CPP);
  TEST_ASSERT(UniToPoly(HN3) == MGHN3_CPP);
  
  ideal Q4 = NewQueenIdeal(CsbRing, 4);
  RingElem HN4 = HilbertNumQuot_C(Q4);
  RingElem HN4_CPP = HilbertNumQuot(Q4);
  RingElem MGHN4_CPP = MGHilbertNumQuot(Q4);
  TEST_ASSERT(HN4 == HN4_CPP);
  TEST_ASSERT(UniToPoly(HN4) == MGHN4_CPP);

  ideal Q6 = NewQueenIdeal(CsbRing, 6);
///  t0 = CpuTime();
  RingElem HN = HilbertNumQuot_C(Q6);
  TEST_ASSERT(LC(HN) == 19);
  TEST_ASSERT(deg(HN) == 36);
///  t0 = CpuTime();
  RingElem HN_CPP = HilbertNumQuot(Q6);
  TEST_ASSERT(LC(HN_CPP) == 19);
  TEST_ASSERT(deg(HN_CPP) == 36);
  RingElem MGHN_CPP = MGHilbertNumQuot(Q6);
  TEST_ASSERT(UniToPoly(HN_CPP) == MGHN_CPP);

  //  std::cout << HilbertSeriesQuot(I) << std::endl;
  
  /*
  ideal Q9 = NewQueenIdeal(CsbRing, 9);
  t0 = CpuTime();
  RingElem HN9 = HilbertNumQuot_C(Q9);
  std::clog << CpuTime()-t0 << std::endl;
  t0 = CpuTime();
  RingElem HN9_CPP = HilbertNumQuot(Q9);
  std::clog << CpuTime()-t0 << std::endl;
  TEST_ASSERT(HN9 == HN9_CPP);
  */

  // this line is just to avoid mempool complaints with -DCoCoA_MEMPOOL_DEBUG
  EndPoincare_C();
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/tests/test-hilbert1.C,v 1.24 2014/09/05 16:13:58 abbott Exp $
// $Log: test-hilbert1.C,v $
// Revision 1.24  2014/09/05 16:13:58  abbott
// Summary: Commented out unused variable t0
// Author: JAA
//
// Revision 1.23  2014/07/07 13:40:28  abbott
// Summary: Removed AsPolyRing
// Author: JAA
//
// Revision 1.22  2013/07/30 17:40:26  bigatti
// -- commented out temporary (printing) test
//
// Revision 1.21  2013/07/30 16:53:26  bigatti
// -- simplified input for hilbert
//
// Revision 1.20  2012/04/02 16:45:09  abbott
// Changed CoCoA_ASSERT into TEST_ASSERT.
//
// Revision 1.19  2012/02/10 11:57:12  bigatti
// -- changed RingZ.H, RingQ.H --> RingZZ.H, RingQQ.H
//
// Revision 1.18  2012/02/08 17:36:47  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.17  2011/05/26 16:34:18  bigatti
// -- added test for ideal(0)
//
// Revision 1.16  2011/05/24 14:56:58  abbott
// Consequential changes from removal of several ctors for principal ideals.
//
// Revision 1.15  2011/04/27 09:48:37  bigatti
// -- updated syntax
//
// Revision 1.14  2011/04/26 10:32:21  bigatti
// -- added tests for multigraded case
//
// Revision 1.13  2011/04/08 14:07:37  bigatti
// -- renamed HilbertNumeratorMod into HilbertNumQuot
//
// Revision 1.12  2011/03/10 17:58:33  bigatti
// -- using long instead of size_t
//
// Revision 1.11  2010/12/17 16:06:25  abbott
// Ensured that all i/o is on standard C++ streams (instead of GlobalInput, etc)
//
// Revision 1.10  2010/10/29 09:43:09  bigatti
// -- manually freeing global memory for C implementation
//    to removed complaints when compiled with -DCoCoA_MEMPOOL_DEBUG
//
// Revision 1.9  2010/10/08 08:19:51  bigatti
// -- RingDistrMPoly.H --> RingDistrMPolyClean.H
//
// Revision 1.8  2009/07/30 15:41:25  bigatti
// -- now using the new nice constructors for ideals
//
// Revision 1.7  2007/10/30 17:14:05  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.6  2007/10/19 10:04:23  bigatti
// -- RingDenseUPolyClean now allow to specify the MinCapacity for all
//    coeff vectors (to avoid too many reallocations)
//
// Revision 1.5  2007/10/18 11:37:33  bigatti
// -- added Q9 for timings tests
//
// Revision 1.4  2007/10/15 12:59:58  bigatti
// -- added include
//
// Revision 1.3  2007/10/15 12:45:59  bigatti
// -- HP computed in Z[lambda] instead of Q[lambda]
//
// Revision 1.2  2007/10/10 14:40:48  bigatti
// -- added: test for comparing old C code with development version using CoCoAlib
//
// Revision 1.1.1.1  2007/03/09 15:16:12  abbott
// Imported files
//
// Revision 1.9  2007/03/08 17:41:36  bigatti
// -- improved NewPolyRing
//
// Revision 1.8  2007/03/08 14:38:07  cocoa
// Added new range function in symbol.H, and tidied many calls to PolyRing
// pseudo ctors (as a consequence).
//
// Revision 1.7  2007/03/03 14:13:21  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.6  2007/03/02 17:46:40  bigatti
// -- unique RingZ and RingQ
// -- requires foundations.H ;  foundations blah;  (thik of a better name)
//
// Revision 1.5  2007/02/26 17:11:58  bigatti
// -- getting ready for unique ring Z: using NewZmod(N), NewRingQ()
//
// Revision 1.4  2007/02/10 18:44:02  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.3  2007/01/17 17:38:11  bigatti
// -- moved all cocoa-4 code for hilbert into src/TmpHilbertDir
//
// Revision 1.2  2006/12/07 17:25:32  cocoa
// -- minimal set of #include's instead of library.H
//
// Revision 1.1  2006/11/16 18:16:10  cocoa
// -- added test for HilbertNumQuot
//
