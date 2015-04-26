//   Copyright (c)  2009  Anna Bigatti

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
#include "CoCoA/FractionField.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/PPWithMask.H"
#include "CoCoA/PPMonoidEv.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingDistrMPolyClean.H"
#include "CoCoA/RingDistrMPolyInlFpPP.H"
#include "CoCoA/RingDistrMPolyInlPP.H"
#include "CoCoA/RingFp.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/TmpPPVector.H"
#include "CoCoA/assert.H"  // for using TEST_ASSERT
#include "CoCoA/error.H"
#include "CoCoA/symbol.H"

using namespace CoCoA;

#include <algorithm>
using std::min;
#include <iostream>
using std::cerr;
using std::endl;
#include <vector>
using std::vector;

//----------------------------------------------------------------------
// Test for PPVector operations
// see also test-MonomialIdeal1
// functions: 
// environments: 
//----------------------------------------------------------------------

// Handy macro for making assertions.
#define TEST_ASSERT(cond)   do { if (!(cond)) CoCoA::AssertionFailed(#cond, __FILE__, __LINE__); } while(0)

ConstRefPPMonoidElem CsbSquareIndet(PPMonoid PPM, long l, long sq1, long sq2)
{
  TEST_ASSERT( l*l <= NumIndets(PPM) );
  TEST_ASSERT( sq1 <= l && sq2 <= l );
  return indet(PPM, (sq1-1)*l + (sq2-1));
}


void AppendQueenMovesFrom(PPVector& PPV, long Csb, long sq1, long sq2)
{
  ConstRefPPMonoidElem x = CsbSquareIndet(PPM(PPV), Csb, sq1, sq2);
  for ( long i=sq2+1 ; i<=Csb ; ++i )
    PushBack(PPV, x * CsbSquareIndet(PPM(PPV), Csb, sq1, i));
  for ( long i=sq1+1 ; i<=Csb ; ++i )
    PushBack(PPV, x * CsbSquareIndet(PPM(PPV), Csb, i, sq2));
  for ( long i=min(Csb-sq1,Csb-sq2) ; i>0 ; --i )
    PushBack(PPV, x * CsbSquareIndet(PPM(PPV), Csb, sq1+i, sq2+i));
  for ( long i=min(Csb-sq1, sq2-1) ; i>0 ; --i )
    PushBack(PPV, x * CsbSquareIndet(PPM(PPV), Csb, sq1+i, sq2-i));
}


void QueenPPVector(PPVector& PPV, long Csb)
{
  for ( long sq1=1 ; sq1<=Csb ; ++sq1 )
    for ( long sq2=1 ; sq2<=Csb ; ++sq2 )
      AppendQueenMovesFrom(PPV, Csb, sq1, sq2);
}


// ----------------------------------------------------------------------

void program()
{
  GlobalManager CoCoAFoundations;

  int N = 7;
  PPMonoid PPM1 = NewPPMonoidEv(SymbolRange("x", 0, N*N-1), lex);
  DivMaskRule DMR1 = NewDivMaskEvenPowers();

  PPVector PPV(PPM1, DMR1);
  QueenPPVector(PPV, N);
  interreduce(PPV);


  const std::vector<PPMonoidElem> x = indets(PPM1);
  PPVector f(PPM1, DMR1);
  PPVector g(PPM1, DMR1);
  PPVector h(PPM1, DMR1);
  PushBack(f, x[1]*x[2]);  PushBack(f, x[3]*x[2]);
  PushBack(g, x[1]*x[5]);  PushBack(g, x[6]*x[7]);
  
  lcms(h, f, g);

  ring P = NewPolyRing(RingQQ(), NumIndets(PPM1));
  std::vector<RingElem> v;
  convert(v, P, h);
  PPVector l(PPM1, DMR1);
  convert(l, v);

  TEST_ASSERT(len(l) == len(h));
  TEST_ASSERT(PPM(l) == PPM(h));
  TEST_ASSERT(owner(PP(l[0])) == owner(PP(h[0])));
  for (long i=0; i<len(l) ; ++i)    TEST_ASSERT(l[i] == h[i]);
}


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
