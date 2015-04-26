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
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/DenseUPolyRing.H"
#include "CoCoA/FractionField.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/MatrixOperations.H"
#include "CoCoA/RingDenseUPolyClean.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/random.H"
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
// Test for matrix
// functions: DeleteRow, DeleteCol
// environments: DenseMatrix over Z and Q
//----------------------------------------------------------------------

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


void TestMatrix(ConstMatrixView origM, ConstMatrixView Ma,
                ConstMatrixView Mb, ConstMatrixView Mc)
{
  matrix M(NewDenseMat(origM));
  //  cout << "M is " << M << endl;
  DeleteCol(M,2);
  //  cout << "after DeleteCol(M,2) is " << M << endl;
  //  cout << " " << Ma << endl;
  TEST_ASSERT(M==Ma);
  DeleteRow(M,0);
  //  cout << "after DeleteRow(M,0) is " << M << endl;
  TEST_ASSERT(M==Mb);
  DeleteRow(M,2);
  //  cout << "after DeleteRow(M,2) is " << M << endl;
  TEST_ASSERT(M==Mc);
  //  cout << "---------------------------------------------" << endl;
}


void program()
{
  GlobalManager CoCoAFoundations;

  const ring ZZ = RingZZ();
  const ring QQ = RingQQ();

  int M1[4*4] = {1, 0, 0, 4,
                 0, 3, 2, 0,
                 0, 3, 2, 1,
                 1, 0, 0, 0};
  int M1a[4*3] = {1, 0,  4,    0, 3,  0,    0, 3,  1,    1, 0,  0};
  int M1b[3*3] = {             0, 3,  0,    0, 3,  1,    1, 0,  0};
  int M1c[3*3] = {             0, 3,  0,    0, 3,  1};
  TestMatrix(NewMatrixFromC(ZZ, M1, 4,4),
             NewMatrixFromC(ZZ, M1a, 4,3),
             NewMatrixFromC(ZZ, M1b, 3,3),
             NewMatrixFromC(ZZ, M1c, 2,3));
  TestMatrix(NewMatrixFromC(QQ, M1, 4,4),
             NewMatrixFromC(QQ, M1a, 4,3),
             NewMatrixFromC(QQ, M1b, 3,3),
             NewMatrixFromC(QQ, M1c, 2,3));
  int M2[4*4] = {1, 0, 0, 0,
                 0, 3, 2, 0,
                 0, 3, 2, 1,
                 1, 1, 1, 0};
  int M2a[4*3] = {1, 0,  0,    0, 3,  0,    0, 3,  1,    1, 1,  0};
  int M2b[3*3] = {             0, 3,  0,    0, 3,  1,    1, 1,  0};
  int M2c[3*3] = {             0, 3,  0,    0, 3,  1};
  ring R = NewPolyRing_DUP(ZZ);
  TestMatrix(NewMatrixFromC(R, M2, 4,4),
             NewMatrixFromC(R, M2a, 4,3),
             NewMatrixFromC(R, M2b, 3,3),
             NewMatrixFromC(R, M2c, 2,3));
  //  TestMatrix(NewMatrixFromC(NewPolyRing(QQ,3), M2, 4,4));
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
