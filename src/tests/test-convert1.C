//   Copyright (c)  2007,2009  John Abbott

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


#include "CoCoA/IntOperations.H"
#include "CoCoA/BigRat.H"
#include "CoCoA/BuildInfo.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/convert.H"

using namespace CoCoA;

#include <limits>
using std::numeric_limits;
#include <string>
using std::string;
#include <iostream>
using std::cerr;
using std::endl;
#include <cmath>
using std::pow;

#define TEST_ASSERT(cond)   do { if (!(cond)) CoCoA::AssertionFailed(#cond, __FILE__, __LINE__); } while(0)


// Tests for the various conversion functions.

void program()
{
  GlobalManager CoCoAFoundations;

  // Conversion from BigInt to machine integers
  BigInt N;
  {
    // Conversion from BigInt to (un)signed long.
    long sl;
    N = numeric_limits<long>::max();
    TEST_ASSERT(IsConvertible(sl, N) && sl == N);
    TEST_ASSERT(!IsConvertible(sl, N+1));
    N = numeric_limits<long>::min();
    TEST_ASSERT(IsConvertible(sl, N) && sl == N);
    TEST_ASSERT(!IsConvertible(sl, N-1));

    unsigned long ul;
    N = numeric_limits<unsigned long>::max();
    TEST_ASSERT(IsConvertible(ul, N) && ul == N);
    TEST_ASSERT(!IsConvertible(ul, N+1));
    N = numeric_limits<unsigned long>::min(); // just zero
    TEST_ASSERT(IsConvertible(ul, N) && ul == N);
    TEST_ASSERT(!IsConvertible(ul, N-1));


    // Conversion from BigInt to (un)signed int.
    int si;
    N = numeric_limits<int>::max();
    TEST_ASSERT(IsConvertible(si, N) && si == N);
    TEST_ASSERT(!IsConvertible(si, N+1));
    N = numeric_limits<int>::min();
    TEST_ASSERT(IsConvertible(si, N) && si == N);
    TEST_ASSERT(!IsConvertible(si, N-1));

    unsigned int ui;
    N = numeric_limits<unsigned int>::max();
    TEST_ASSERT(IsConvertible(ui, N) && ui == N);
    TEST_ASSERT(!IsConvertible(ui, N+1));
    N = numeric_limits<unsigned int>::min(); // just zero
    TEST_ASSERT(IsConvertible(ui, N) && ui == N);
    TEST_ASSERT(!IsConvertible(ui, N-1));


// !!! NYI Conversion from BigInt to (un)signed short.
//   short ss;
//   N = numeric_limits<short>::max();
//   TEST_ASSERT(IsConvertible(ss, N) && ss == N);
//   TEST_ASSERT(!IsConvertible(ss, N+1));
//   N = numeric_limits<short>::min();
//   TEST_ASSERT(IsConvertible(ss, N) && ss == N);
//   TEST_ASSERT(!IsConvertible(ss, N-1));

//   unsigned short us;
//   N = numeric_limits<unsigned short>::max();
//   TEST_ASSERT(IsConvertible(us, N) && us == N);
//   TEST_ASSERT(!IsConvertible(us, N+1));
//   N = numeric_limits<unsigned short>::min(); // just zero
//   TEST_ASSERT(IsConvertible(us, N) && us == N);
//   TEST_ASSERT(!IsConvertible(us, N-1));


// !!! NYI Conversion from BigInt to (un)signed char
//   signed char sc;
//   N = numeric_limits<signed char>::max();
//   TEST_ASSERT(IsConvertible(sc, N) && sc == N);
//   TEST_ASSERT(!IsConvertible(sc, N+1));
//   N = numeric_limits<signed char>::min();
//   TEST_ASSERT(IsConvertible(sc, N) && sc == N);
//   TEST_ASSERT(!IsConvertible(sc, N-1));

//   unsigned char uc;
//   N = numeric_limits<unsigned char>::max();
//   TEST_ASSERT(IsConvertible(uc, N) && uc == N);
//   TEST_ASSERT(!IsConvertible(uc, N+1));
//   N = numeric_limits<unsigned char>::min(); // just zero
//   TEST_ASSERT(IsConvertible(uc, N) && uc == N);
//   TEST_ASSERT(!IsConvertible(uc, N-1));
  }


  // Conversion from BigRat to machine integers
  BigRat Q;
  {
    const BigRat half = BigRat(1,2);
    // Conversion from BigRat to (un)signed long.
    long sl;
    Q = numeric_limits<long>::max();
    TEST_ASSERT(IsConvertible(sl, Q) && sl == Q);
    TEST_ASSERT(!IsConvertible(sl, Q+1));
    TEST_ASSERT(!IsConvertible(sl, Q-half));
    Q = numeric_limits<long>::min();
    TEST_ASSERT(IsConvertible(sl, Q) && sl == Q);
    TEST_ASSERT(!IsConvertible(sl, Q-1));
    TEST_ASSERT(!IsConvertible(sl, Q+half));

    unsigned long ul;
    Q = numeric_limits<unsigned long>::max();
    TEST_ASSERT(IsConvertible(ul, Q) && ul == Q);
    TEST_ASSERT(!IsConvertible(ul, Q+1));
    TEST_ASSERT(!IsConvertible(ul, Q-half));
    Q = numeric_limits<unsigned long>::min(); // just zero
    TEST_ASSERT(IsConvertible(ul, Q) && ul == Q);
    TEST_ASSERT(!IsConvertible(ul, Q-1));
    TEST_ASSERT(!IsConvertible(ul, Q+half));


    // Conversion from BigRat to (un)signed int.
    int si;
    Q = numeric_limits<int>::max();
    TEST_ASSERT(IsConvertible(si, Q) && si == Q);
    TEST_ASSERT(!IsConvertible(si, Q+1));
    TEST_ASSERT(!IsConvertible(si, Q-half));
    Q = numeric_limits<int>::min();
    TEST_ASSERT(IsConvertible(si, Q) && si == Q);
    TEST_ASSERT(!IsConvertible(si, Q-1));
    TEST_ASSERT(!IsConvertible(si, Q+half));

    unsigned int ui;
    Q = numeric_limits<unsigned int>::max();
    TEST_ASSERT(IsConvertible(ui, Q) && ui == Q);
    TEST_ASSERT(!IsConvertible(ui, Q+1));
    TEST_ASSERT(!IsConvertible(ui, Q-half));
    Q = numeric_limits<unsigned int>::min(); // just zero
    TEST_ASSERT(IsConvertible(ui, Q) && ui == Q);
    TEST_ASSERT(!IsConvertible(ui, Q-1));
    TEST_ASSERT(!IsConvertible(ui, Q+half));


// !!! QYI Conversion from BigRat to (un)signed short.
//   short ss;
//   Q = numeric_limits<short>::max();
//   TEST_ASSERT(IsConvertible(ss, Q) && ss == Q);
//   TEST_ASSERT(!IsConvertible(ss, Q+1));
//   TEST_ASSERT(!IsConvertible(ss, Q-half));
//   Q = numeric_limits<short>::min();
//   TEST_ASSERT(IsConvertible(ss, Q) && ss == Q);
//   TEST_ASSERT(!IsConvertible(ss, Q-1));
//   TEST_ASSERT(!IsConvertible(ss, Q+half));

//   unsigned short us;
//   Q = numeric_limits<unsigned short>::max();
//   TEST_ASSERT(IsConvertible(us, Q) && us == Q);
//   TEST_ASSERT(!IsConvertible(us, Q+1));
//   TEST_ASSERT(!IsConvertible(us, Q-half));
//   Q = numeric_limits<unsigned short>::min(); // just zero
//   TEST_ASSERT(IsConvertible(us, Q) && us == Q);
//   TEST_ASSERT(!IsConvertible(us, Q-1));
//   TEST_ASSERT(!IsConvertible(us, Q+half));


// !!! QYI Conversion from BigRat to (un)signed char
//   signed char sc;
//   Q = numeric_limits<signed char>::max();
//   TEST_ASSERT(IsConvertible(sc, Q) && sc == Q);
//   TEST_ASSERT(!IsConvertible(sc, Q+1));
//   TEST_ASSERT(!IsConvertible(sc, Q-half));
//   Q = numeric_limits<signed char>::min();
//   TEST_ASSERT(IsConvertible(sc, Q) && sc == Q);
//   TEST_ASSERT(!IsConvertible(sc, Q-1));
//   TEST_ASSERT(!IsConvertible(sc, Q+half));

//   unsigned char uc;
//   Q = numeric_limits<unsigned char>::max();
//   TEST_ASSERT(IsConvertible(uc, Q) && uc == Q);
//   TEST_ASSERT(!IsConvertible(uc, Q+1));
//   TEST_ASSERT(!IsConvertible(uc, Q-half));
//   Q = numeric_limits<unsigned char>::min(); // just zero
//   TEST_ASSERT(IsConvertible(uc, Q) && uc == Q);
//   TEST_ASSERT(!IsConvertible(uc, Q-1));
//   TEST_ASSERT(!IsConvertible(uc, Q+half));
  }



  // Conversions to and from doubles -- need to be careful as the precision of doubles may be platform dependent
  TEST_ASSERT(!IsConvertible(N, 1.25));
  TEST_ASSERT(IsConvertible(Q, 1.25) && Q == BigRat(5,4));
  TEST_ASSERT(IsConvertible(N, 12345.0) && N == 12345);
  TEST_ASSERT(IsConvertible(Q, 12345.0) && Q == 12345);
  double z = std::pow(2.0, 200) + 1.0;
  TEST_ASSERT(IsConvertible(N, z) && N == power(2,200));
  TEST_ASSERT(IsConvertible(Q, z) && Q == power(2,200));
  z = 1.0/3.0;
  TEST_ASSERT(!IsConvertible(N, z));
  TEST_ASSERT(IsConvertible(Q,z) && Q != BigRat(1,3) && den(Q)%2 == 0);
  TEST_ASSERT(IsConvertible(z, power(2,200)) && z == std::pow(2.0,200));
  TEST_ASSERT(IsConvertible(z, BigRat(5,4)) && z == 1.25);
  TEST_ASSERT(IsConvertible(z, BigRat(4,3)) && std::abs(3*z-4) < 1.0e-15);

  // Checking overflow handling
  N = power(5,100000);
  TEST_ASSERT(!IsConvertible(z, N)); // conversion fails due to overflow
  // Now test conversion of rational with large num & den, but small value.
  TEST_ASSERT(IsConvertible(z, BigRat(5*N+1, 4*N)) && z == 1.25);
  TEST_ASSERT(IsConvertible(z, -BigRat(5*N+1, 4*N)) && z == -1.25);
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
