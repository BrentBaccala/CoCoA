//   Copyright (c)  2010  John Abbott

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
//#include "CoCoA/assert.H"  // for using TEST_ASSERT
#include "CoCoA/convert.H"
#include "CoCoA/error.H"


#include <iostream>
using std::clog;
using std::cerr;
using std::endl;
#include<limits>
using std::numeric_limits;

using namespace CoCoA;

// Handy macro for making assertions.
#define TEST_ASSERT(cond)   do { if (!(cond)) CoCoA::AssertionFailed(#cond, __FILE__, __LINE__); } while(0)

void program()
{
  // This test checks the functions iroot, isqrt, IsSquare & IsPower.
  GlobalManager CoCoAFoundations(UseGMPAllocator);

  BigInt ROOT;   // used as arg to IsExactIroot
  long root; // used as arg to IsExactIroot
  BigInt N;
  for (int a=-99; a <= 99; ++a)
    for (int b=2; b < 99; ++b)
    {
      N = power(a, b);
      TEST_ASSERT(IsExactIroot(ROOT,N,b));
      {
        long n;
        if (IsConvertible(n, N))
        {
          TEST_ASSERT(IsExactIroot(root,n,b));
          TEST_ASSERT(root == ROOT);
        }
      }
      if (a >= 0 || IsOdd(b))
      {
        TEST_ASSERT(a == iroot(N,b));
        TEST_ASSERT(a == ROOT);
      }
      else
      {
        TEST_ASSERT(-a == iroot(N,b));
        TEST_ASSERT(-a == ROOT);
      }
      TEST_ASSERT(IsPower(N));
      if (!((N == -1) || (N == 0) || (N == 8)))
      {
        TEST_ASSERT(!IsPower(N+1));
        TEST_ASSERT(!IsExactIroot(ROOT, N+1,b));
        if (a > 0 || IsEven(b)) TEST_ASSERT(ROOT == abs(a));
        else TEST_ASSERT(ROOT == a+1); // recall that a < 0 here!!
      }
      if (!((N == 0) || (N == 1) || (N == 9)))
      {
        TEST_ASSERT(!IsPower(N-1));
        TEST_ASSERT(!IsExactIroot(ROOT, N-1,b));
        if (a > 0 || IsEven(b)) TEST_ASSERT(ROOT == abs(a)-1);
        else TEST_ASSERT(ROOT == a); // recall that a < 0 here!!
      }
      if (IsEven(b))
      {
        TEST_ASSERT(IsSquare(N));
        TEST_ASSERT(isqrt(N) == power(abs(a),b/2));
      }
      if (N > 0) TEST_ASSERT(abs(a) == iroot(N+1,b));
      if (N < 0) TEST_ASSERT(a == iroot(N-1,b));
      long n;
      if (IsConvertible(n, N))
      {
        TEST_ASSERT(IsPower(n));
        if (!((n == -1) || (n == 0) || (n == 8))) TEST_ASSERT(!IsPower(n+1));
        if (!((n == 0) || (n == 1) || (n == 9) || (n == numeric_limits<long>::min()))) TEST_ASSERT(!IsPower(n-1));
        if (IsEven(b))
        {
          TEST_ASSERT(IsSquare(n));
          TEST_ASSERT(isqrt(n) == power(abs(a),b/2));
        }
        TEST_ASSERT(iroot(n,b) == iroot(N,b));
        if (n > 0) TEST_ASSERT(abs(a) == iroot(n+1,b)); // n+1 cannot overflow (unless 2^k-1 is a perfect power)
        if (n < 0 && n != numeric_limits<long>::min()) TEST_ASSERT(a == iroot(n-1,b));
      }
      try { IsExactIroot(ROOT,N,0); TEST_ASSERT(!"NEVER GET HERE!"); }
      catch (const ErrorInfo& err) { TEST_ASSERT(err == ERR::BadArg);}
      try { IsExactIroot(ROOT,N,-b); TEST_ASSERT(!"NEVER GET HERE!"); }
      catch (const ErrorInfo& err) { TEST_ASSERT(err == ERR::BadArg);}
      if (IsEven(b) && N > 0)
        try { IsExactIroot(ROOT,-N,b); TEST_ASSERT(!"NEVER GET HERE!"); }
        catch (const ErrorInfo& err) { TEST_ASSERT(err == ERR::BadArg);}
    }
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
