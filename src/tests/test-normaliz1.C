//   Copyright (c)  2009 Anna Bigatti
//   Copyright (c)  2012-2014 Christof Soeger

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
#include "CoCoA/IntOperations.H"
#include "CoCoA/BigInt.H"
#include "CoCoA/BigRat.H"
#include "CoCoA/error.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/RingElemInput.H"
#include "CoCoA/TmpHilbert.H"
#include "CoCoA/QuasiPoly.H"

#include "CoCoA/utils.H"
#include "CoCoA/ExternalLibs-Normaliz.H"


using namespace CoCoA;

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <vector>
using std::vector;

//----------------------------------------------------------------------
// First test for Normaliz library.
// Simple computation to test proper integration.
// This tests is for the functions working directly on cone.
//----------------------------------------------------------------------

#define TEST_ASSERT(cond)   do { if (!(cond)) CoCoA::AssertionFailed(#cond, __FILE__, __LINE__); } while(0)

vector<BigInt> BigIntVec(const int* CVector, int len)
{
  vector<BigInt> v(len);
  for (int i=0; i<len; ++i)  v[i] = (BigInt(CVector[i]));
  return v;
}

void program()
{
  GlobalManager CoCoAFoundations;

#ifdef CoCoA_WITH_NORMALIZ
  using namespace CoCoA::Normaliz;

// this is the polytope example from the Normaliz examples
  const int M[4][4] = {{0, 0, 0, 1},
                       {2, 0, 0, 1},
                       {0, 3, 0, 1},
                       {0, 0, 5, 1}};
  const int g[4] = {0, 0, 0, 1};
  vector<vector<BigInt> > l;
  for (int i=0; i<4; ++i)
    l.push_back(BigIntVec(M[i], 4));

//  cout << "l -> " << len(l) << endl;

  std::map< libnormaliz::InputType, std::vector<std::vector<BigInt> > > m;
  m[libnormaliz::Type::integral_closure] = l;
  m[libnormaliz::Type::grading] = vector<vector<BigInt> >(1,BigIntVec(g, 4));

  vector<vector<BigInt> > res_m;
  vector<BigInt>  res_v;

  cone C(m);
  vector<vector<BigInt> > hb = HilbertBasis(C);

  TEST_ASSERT(len(hb) == 19);
  for (long i=0; i<len(hb); ++i)
  {
    TEST_ASSERT(len(hb[i]) == 4);
    //check if the points are inside the cone by evaluating the support hyperplanes
    TEST_ASSERT(hb[i][0] >= 0);
    TEST_ASSERT(hb[i][1] >= 0);
    TEST_ASSERT(hb[i][2] >= 0);
    TEST_ASSERT(30*hb[i][3]-15*hb[i][0]-10*hb[i][1]-6*hb[i][2] >= 0);
  }

  const vector< vector<BigInt> > deg1 = Deg1Elements(C);
  TEST_ASSERT(len(deg1) == 18);

  const vector< vector<BigInt> > sh = SupportHyperplanes(C);
  TEST_ASSERT(len(sh) == 4);

  const vector< vector<BigInt> > rays = ExtremeRays(C);
  TEST_ASSERT(len(rays) == 4);

  const vector< vector<BigInt> > equ = Equations(C);
  TEST_ASSERT(len(equ) == 0);

  const vector< vector<BigInt> > cong = Congruences(C);
  TEST_ASSERT(len(cong) == 0);

  TEST_ASSERT(IsPointed(C));
  TEST_ASSERT(!IsInhomogeneous(C));
  TEST_ASSERT(!IsIntegrallyClosed(C));
  TEST_ASSERT(!IsDeg1HilbertBasis(C));

  TEST_ASSERT(multiplicity(C) == 30);

  HPSeries hs = HilbertSeries(C);
  // compare it with:  (1 + 14*t + 15*t^2) / (1-t)^4
  PolyRing QQt = RingQQt(1);
  const RingElem t = indet(QQt,0);
  RingElem ref_num(ReadExpr(QQt, "1 + 14*t + 15*t^2"));
  TEST_ASSERT(num(hs) == ref_num);
  factorization<RingElem> den = DenFactors(hs);
  TEST_ASSERT(den.myFactors() == vector<RingElem>(1,ReadExpr(QQt, "1-t")));
  TEST_ASSERT(den.myMultiplicities() == std::vector<long>(1,4));
  TEST_ASSERT(den.myRemainingFactor() == one(QQt));

  RingElem hp = HilbertPoly(C);
  // compare with Hilbert polynomial:   1 +4*t +8*t^2 +5*t^3
  RingElem ref_hp = ReadExpr(QQt, "1 + 4*t + 8*t^2 +5*t^3");
  TEST_ASSERT(hp == ref_hp);

  // test QuasiPoly of period 1
  const QuasiPoly qp = HilbertQuasiPoly(C);
  vector<RingElem> qpv = constituents(qp);
  TEST_ASSERT(len(qpv) == 1);
  TEST_ASSERT(qpv[0] == ref_hp);
  TEST_ASSERT(qp(BigInt(0)) == 1);
  TEST_ASSERT(qp(BigInt(1)) == 18); // deg1
  TEST_ASSERT(qp(BigInt(2)) == 81);

  // new very simple example not generated in deg1


  std::map< libnormaliz::InputType, std::vector<std::vector<BigInt> > > m2;
  vector<vector<BigInt> > l2 = vector<vector<BigInt> >(2,vector<BigInt>(2));
  l2[0][0] = 1; l2[0][1] = 2;
  l2[1][0] = 2; l2[1][1] = 1;
  m2[libnormaliz::Type::integral_closure] = l2;
  m2[libnormaliz::Type::grading] = vector< vector<BigInt> >(1,vector<BigInt>(2,BigInt(1)));


  cone C2(m2);

  const QuasiPoly qp2 = HilbertQuasiPoly(C2);
  //cout << qp2;
  const vector<RingElem>& qpv2 = constituents(qp2);
  TEST_ASSERT(len(qpv2) == 3);
  TEST_ASSERT(qp2(BigInt(0)) == 1);
  TEST_ASSERT(qp2(BigInt(1)) == 0); // deg1
  TEST_ASSERT(qp2(BigInt(2)) == 1);
  TEST_ASSERT(qp2(BigInt(3)) == 2);
  TEST_ASSERT(qp2(BigInt(4)) == 1);
  TEST_ASSERT(qp2(BigInt(5)) == 2);
  TEST_ASSERT(qp2(BigInt(6)) == 3);

#endif // #ifdef CoCoA_WITH_NORMALIZ
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/tests/test-normaliz1.C,v 1.21 2014/07/14 13:17:13 abbott Exp $
// $Log: test-normaliz1.C,v $
// Revision 1.21  2014/07/14 13:17:13  abbott
// Summary: Added two consts
// Author: JAA
//
// Revision 1.20  2014/07/14 12:24:04  abbott
// Summary: Added a const
// Author: JAA
//
// Revision 1.19  2014/07/14 10:08:23  abbott
// Summary: Christof added test for quasipolys
// Author: JAA
//
// Revision 1.18  2014/06/17 10:22:07  abbott
// Summary: Removed pointless call to AsPolyRing (on RingQQt(...))
// Author: JAA
//
// Revision 1.17  2014/05/12 10:20:42  abbott
// Summary: Added some consts
// Author: JAA
//
// Revision 1.16  2014/05/09 14:56:58  bigatti
// -- new fn by Christof Soeger
//
// Revision 1.15  2012/10/08 13:53:44  bigatti
// -- more cleaning and updates by Christof Soeger
//
// Revision 1.14  2012/08/03 16:31:22  bigatti
// -- changed: procedural --> functional (by C.Soeger)
//
// Revision 1.13  2012/07/19 17:15:28  abbott
// NEEDS TO BE REWRIITEN.
// Replaced calls to NewConeLong and NewConeBigInt by calls to unified pseudo-ctor NewCone.
//
// Revision 1.12  2012/06/19 16:11:00  bigatti
// -- changed Ht1 --> Deg1 (by C.Soeger)
//
// Revision 1.11  2012/05/28 09:18:20  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.10  2012/03/12 11:36:42  abbott
// Added two "const"s and improved the indentation.
//
// Revision 1.9  2011/11/07 11:29:59  bigatti
// -- changed syntax for HilbertBasis
//
// Revision 1.8  2011/10/12 15:49:21  abbott
// Simplified use of CPP macro CoCoA_WITH_NORMALIZ, but compilation may be a little slower.
//
// Revision 1.7  2011/09/30 15:57:54  bigatti
// -- moved namespace CoCoA::Normaliz inside #ifdef
//
// Revision 1.6  2011/09/30 12:55:44  bigatti
// -- introduced namespace "Normaliz" and removed Normaliz from function names
// -- input of Normaliz functions in CoCoA-5 is now a matrix instead of
//    vector<vector<BigInt>>
//
// Revision 1.5  2011/08/23 06:40:30  bigatti
// -- fixed with new name "BigInt" for old "ZZ"
//
// Revision 1.4  2011/07/20 13:49:37  bigatti
// -- added "Normaliz" postfix to Normaliz function calls
//
// Revision 1.3  2011/07/20 12:45:12  bigatti
// -- new normaliz interface (not yet public)
//
// Revision 1.2  2011/07/20 10:12:17  bigatti
// -- fixed compilation without normaliz
//
// Revision 1.1  2011/07/19 16:24:16  bigatti
// -- first import
//
// Revision 1.8  2011/07/05 15:15:33  bigatti
// -- changed AlexanderDual --> AlexanderDualFrobby
// -- changed PrimaryDecomposition --> PrimaryDecompositionFrobby
//
// Revision 1.7  2010/12/17 16:06:25  abbott
// Ensured that all i/o is on standard C++ streams (instead of GlobalInput, etc)
//
// Revision 1.6  2010/11/26 15:32:16  bigatti
// -- renamed Dimension into dimension (coding conventions)
//
// Revision 1.5  2010/11/22 17:39:15  bigatti
// -- updated "TmpFrobby" --> "ExternalLibs-Frobby"
//
// Revision 1.4  2010/02/04 09:38:28  bigatti
// -- new syntax for frobby (more CoCoALib-like)
//
// Revision 1.3  2010/02/03 18:00:22  bigatti
// -- more functions from frobby
//
// Revision 1.2  2009/07/30 15:41:25  bigatti
// -- now using the new nice constructors for ideals
//
// Revision 1.1  2009/02/11 15:08:26  bigatti
// -- first import
//

