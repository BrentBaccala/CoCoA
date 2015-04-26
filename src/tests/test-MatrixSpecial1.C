//   Copyright (c)  2011  Anna Bigatti

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
#include "CoCoA/GlobalManager.H"
#include "CoCoA/MatrixOperations.H"
#include "CoCoA/MatrixSpecial.H"
#include "CoCoA/MatrixView.H"
#include "CoCoA/RingFp.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/error.H"
#include "CoCoA/matrix.H"
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
// Test for special matrices:
// functions: jacobian
//----------------------------------------------------------------------

#define TEST_ASSERT(cond)   do { if (!(cond)) CoCoA::AssertionFailed(#cond, __FILE__, __LINE__); } while(0)


void TestJacobian(const PolyRing& P)
{
  std::vector<RingElem> f;
  long n = NumIndets(P);
  
  matrix M = jacobian(indets(P));
  TEST_ASSERT(NumRows(M)==n && NumCols(M)==n);
  TEST_ASSERT(IsOne(M(0,0)) && IsOne(M(1,1)));
  TEST_ASSERT(IsZero(M(0,1)) && IsZero(M(1,0)));

  M = jacobian(indets(P), indets(P));
  TEST_ASSERT(NumRows(M)==n && NumCols(M)==n);
  TEST_ASSERT(IsOne(M(0,0)) && IsOne(M(1,1)));
  TEST_ASSERT(IsZero(M(0,1)) && IsZero(M(1,0)));

  M = jacobian(indets(P), f);
  TEST_ASSERT(NumRows(M)==n && NumCols(M)==0);

  M = jacobian(f, indets(P));
  TEST_ASSERT(NumRows(M)==0 && NumCols(M)==n);
  //  std::cout << jacobian(f, f) << endl; ---> CoCoA_ERROR
  
}

void TestTensorMat(const ring& R)
{
  {
    matrix M = TensorMat(IdentityMat(RingZZ(), 2), IdentityMat(RingZZ(), 2));
    TEST_ASSERT(M == IdentityMat(RingZZ(), 4));
  }
  {
    ConstMatrixView I1x1 = IdentityMat(R, 1);
    ConstMatrixView A = BlockMat2x2(2*I1x1, 3*I1x1, 5*I1x1, 7*I1x1);
    matrix M = TensorMat(A, A);
    std::vector<long> ZeroOne;
    std::vector<long> TwoThree;
    ZeroOne.push_back(0);    ZeroOne.push_back(1);
    TwoThree.push_back(2);   TwoThree.push_back(3);
    TEST_ASSERT(submat(M, ZeroOne, ZeroOne) == 2 * A);
    TEST_ASSERT(submat(M, ZeroOne, TwoThree) == 3 * A);
    TEST_ASSERT(submat(M, TwoThree, ZeroOne) == 5 * A);
    TEST_ASSERT(submat(M, TwoThree, TwoThree) == 7 * A);
  }
  
  
}




void program()
{
  GlobalManager CoCoAFoundations;

  PolyRing P = NewPolyRing(RingQQ(), symbols("x", "y"));
  TestJacobian(P);
  TestTensorMat(P);
  
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
