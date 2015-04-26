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
#include "CoCoA/FractionField.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/MatrixOperations.H"
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
// functions: det, rank
// environments: DenseMatrix over Z and Q
//----------------------------------------------------------------------

// convention: a function containing a "new" should be named "New.."
matrix NewMatrixFromC(ring K, int cmat[4][4])
{
  const int NumRows = 4;
  const int NumCols = 4;
  matrix M(NewDenseMat(K,NumRows,NumCols));

  for (int i=0; i < NumRows; ++i)
    for (int j=0; j < NumCols; ++j)
      SetEntry(M, i, j, cmat[i][j]);
  return M;
}


void TestMatrix(ConstMatrixView M)
{
  cout << M << endl;
  if (NumRows(M)==NumCols(M))
    cout << "det(M) gives " << det(M) << endl;
  cout << "rank(M) gives " << rank(M) << endl;
  cout << "---------------------------------------------" << endl;
}


void TestMatrixPower(long Nrows, long EntrySize)
{
  ring QQ = RingQQ();
  matrix M = NewDenseMat(QQ,Nrows,Nrows);
  RandomSeqLong RndLong(-EntrySize, EntrySize);
  for (int i=0; i < Nrows; ++i)
    for (int j=0; j < Nrows; ++j)
      SetEntry(M,i,j,NextValue(RndLong));

  for (long N = 1; N < 99; ++N)
  {
    matrix Mpower = power(M,N);

    matrix MMpower = M;
    for (int i=1; i < N; ++i)
      mul(MMpower,MMpower,M);
    if (Mpower != MMpower) cout << "ERROR: Mpower != MMpower!!" << endl;
  }
}


void program()
{
  GlobalManager CoCoAFoundations;

  const ring ZZ = RingZZ();
  const ring QQ = RingQQ();

  int M1[4][4] = {{1, 0, 0, 4},
                  {0, 3, 2, 0},
                  {0, 3, 2, 1},
                  {1, 0, 0, 0}};
  TestMatrix(NewMatrixFromC(ZZ, M1));
  TestMatrix(NewMatrixFromC(QQ, M1));

  int M2[4][4] = {{1, 0, 0, 0},
                  {0, 3, 2, 0},
                  {0, 3, 2, 1},
                  {1, 1, 1, 0}};
  TestMatrix(NewMatrixFromC(ZZ, M2));
  TestMatrix(NewMatrixFromC(QQ, M2));

  TestMatrixPower(1,99);
  TestMatrixPower(5,99);
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
