//   Copyright (c)  2013  John Abbott

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
#include "CoCoA/IdealOfPoints.H"
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/ideal.H"
#include "CoCoA/matrix.H"
#include "CoCoA/PolyRing.H"
#include "CoCoA/ring.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/error.H"
#include "CoCoA/utils.H"


#include <iostream>
using std::cerr;
using std::endl;
#include <vector>
using std::vector;

using namespace CoCoA;

// Handy macro for making assertions.
#define TEST_ASSERT(cond)   do { if (!(cond)) CoCoA::AssertionFailed(#cond, __FILE__, __LINE__); } while(0)

void CheckEvalToZero(const vector<RingElem>& g, const matrix& pts)
{
  const PolyRing P = owner(g[0]);
  const ring k = CoeffRing(P);
  const int NumVars = NumIndets(P);
  const int NumPts = NumRows(pts);

  vector<RingHom> EvalMap;
  for (int i=0; i < NumPts; ++i)
  {
    vector<RingElem> images(NumVars, zero(k));
    for (int j=0; j < NumCols(pts); ++j)
      images[j] = pts(i,j);
    EvalMap.push_back(PolyRingHom(P, k, IdentityHom(k), images));
  }

  for (int i=0; i < len(g); ++i)
    for(int j=0; j < len(EvalMap); ++j)
      TEST_ASSERT(IsZero(EvalMap[j](g[i])));
}

void test1()
{
  ring k = RingQQ();
  SparsePolyRing P = NewPolyRing(k, 3);
  matrix pts = NewDenseMat(k,2,2);
  SetEntry(pts,0,0,1);
  SetEntry(pts,0,1,2);
  SetEntry(pts,1,0,3);
  SetEntry(pts,1,1,4);

  ideal I = IdealOfPoints(P, pts);
  CheckEvalToZero(gens(I), pts);
}

void program()
{
  // This test does virtually nothing, but is a handy template if you want
  // to create your own test code: just copy this file and add your code
  // after the line below -- remember that you must not use CoCoALib functions
  // without first creating a GlobalManager.
  GlobalManager CoCoAFoundations;


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
