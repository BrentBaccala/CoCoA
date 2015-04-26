//   Copyright (c)  2014  John Abbott

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
#include "CoCoA/RingHom.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/SparsePolyRing.H"

using namespace CoCoA;

#include <iostream>
using std::cerr;
using std::endl;


//----------------------------------------------------------------------
// Test for RingHom dangling reference (see redmine #190)
//----------------------------------------------------------------------

#define TEST_ASSERT(cond)   do { if (!(cond)) CoCoA::AssertionFailed(#cond, __FILE__, __LINE__); } while(0)

RingHom CoeffHomDanglingRef(const ring& K)
{
  PolyRing P(NewPolyRing(K,1));
  return CoeffEmbeddingHom(P);
}

void program()
{
  GlobalManager CoCoAFoundations;

  ring QQ = RingQQ();
  RingHom phi = CoeffHomDanglingRef(QQ);
  RingElem img = phi(one(QQ));
  ring QQx = owner(img);
  TEST_ASSERT(IsPolyRing(QQx));
  TEST_ASSERT(CoeffRing(QQx) == QQ);
  TEST_ASSERT(NumIndets(QQx) == 1);
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
