//   Copyright (c)  2009  John Abbott

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
#include "CoCoA/NumTheory.H"
#include "CoCoA/IntOperations.H"
#include "CoCoA/error.H"

#include <algorithm>
using std::max;
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

using namespace CoCoA;

void program()
{
  // The call to max used to cause a SEGV -- an undesired "feature" of the proxy class used to reduce wasteful copying upon returning BigInt values.
  // This bug was originally reported by Max Caboara.
  GlobalManager CoCoAFoundations;

  const BigInt a(6);
  const BigInt b(10);
  const BigInt c(15);

  const BigInt d = max(gcd(a,b), gcd(b,c));
  cout << "Result=" << d << endl;
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
