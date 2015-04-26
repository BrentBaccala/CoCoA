//   Copyright (c)  2012  Anna M Bigatti

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
#include "CoCoA/BigRat.H"
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/MatrixForOrdering.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/error.H"
#include "CoCoA/matrix.H"

using namespace CoCoA;
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

//----------------------------------------------------------------------
// Test for matrix fro orderings
// functions: NewPositiveMat, NewIntegerOrdMat
// environments: ZZ, QQ
//----------------------------------------------------------------------

// Handy macro for making assertions.
#define TEST_ASSERT(cond)   do { if (!(cond)) CoCoA::AssertionFailed(#cond, __FILE__, __LINE__); } while(0)


// convention: a function containing a "new" should be named "New.."
matrix NewMatrixFromC(ring K, int* cmat, long NumRows, long NumCols)
{
  matrix M(NewDenseMat(K,NumRows,NumCols));

  for (long i=0; i < NumRows; ++i)
    for (long j=0; j < NumCols; ++j)
      SetEntry(M, i, j, cmat[i*NumCols+j]);
  return M;
}


void program()
{
  GlobalManager CoCoAFoundations;

  const ring ZZ = RingZZ();
  const ring QQ = RingQQ();

  int M1[4*4] = {1, 0, 0, 4,
                 0, 3, 2, 0,
                 0,-3, 2, 1,
                 1, 0, 0, 0};
  int M1Pos[4*4] =  {1, 0, 0, 4,
                     0, 3, 2, 0,
                     0, 0, 4, 1,
                     1, 0, 0, 0};
  
  //-- ZZ ----
  matrix M = NewMatrixFromC(ZZ, M1, 4,4);
  matrix MPos = NewMatrixFromC(ZZ, M1Pos, 4,4);
  TEST_ASSERT(NewPositiveMat(M) == MPos);
  TEST_ASSERT(NewIntegerOrdMat(M) == M);

  //-- QQ ----
  M = NewMatrixFromC(QQ, M1, 4,4);
  MPos = NewMatrixFromC(QQ, M1Pos, 4,4);
  SetEntry(M,2,3, BigRat(1,3));
  SetEntry(MPos,2,3, BigRat(1,3));
  TEST_ASSERT(NewPositiveMat(M) == MPos);
  matrix MI = M;
  MI->myRowMul(2, RingElem(QQ,3));
  TEST_ASSERT(NewIntegerOrdMat(M) == MI);

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
