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
#include "CoCoA/GlobalManager.H"
#include "CoCoA/IntOperations.H"
#include "CoCoA/NumTheory.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/error.H"
#include "CoCoA/ring.H"

using namespace CoCoA;

#include <iostream>
using std::cerr;
using std::endl;

// This macro is largely copied from CoCoA/config.H
#define TEST_ASSERT(cond)   do { if (!(cond)) CoCoA::AssertionFailed(#cond, __FILE__, __LINE__); } while(0)


void test(ring ZZmodN)
{

  TEST_ASSERT(IsQuotientRing(ZZmodN));
  const QuotientRing QR = ZZmodN;
  TEST_ASSERT(IsZZ(BaseRing(QR)));

  BigInt p = characteristic(ZZmodN);
  TEST_ASSERT(IsProbPrime(p) == IsField(ZZmodN));
  TEST_ASSERT(IsZero(RingElem(ZZmodN, p)));

  RingElem n(ZZmodN);

  TEST_ASSERT(IsZero(n));
  TEST_ASSERT(!IsInvertible(n));
  TEST_ASSERT(n == 0);
  TEST_ASSERT(0 == n);
  TEST_ASSERT(n == zero(ZZmodN));

  // We cannot do inequalities: (not even in ZZ/(0))
  try { n > 0; TEST_ASSERT(!"NEVER GET HERE!"); }
  catch (const CoCoA::ErrorInfo& err) { TEST_ASSERT(err == ERR::NotOrdDom); }

  BigInt N = BigInt(1);
  n = N;
  TEST_ASSERT(IsInvertible(n));
  TEST_ASSERT(IsOne(n));
  TEST_ASSERT(IsMinusOne(-n));
  TEST_ASSERT(n == 1);
  TEST_ASSERT(1 == n);
  TEST_ASSERT(n == N);
  TEST_ASSERT(N == n);
  n = N;
  TEST_ASSERT(IsOne(n));

  // Some simple arithmetic
  n = n+1; TEST_ASSERT(n == 2);
  n = n-1; TEST_ASSERT(n == 1);
  n = n*1; TEST_ASSERT(n == 1);
  n = n/1; TEST_ASSERT(n == 1);
  n = 1+n; TEST_ASSERT(n == 2);
  n = 1-n; TEST_ASSERT(n == -1);
  n = 1*n; TEST_ASSERT(n == -1);
  n = 1/n; TEST_ASSERT(n == -1);
  n = n+n; TEST_ASSERT(n == -2);
  n = n-n; TEST_ASSERT(n == 0);
  n = n*n; TEST_ASSERT(n == 0);
  n = n/(n+1); TEST_ASSERT(n == 0);

  // assignment arithmetic
  n = p-1;
  RingElem m = n;

  n += m;
  n -= m;
  n *= m;
  n /= m;

  n = m; n += n; TEST_ASSERT(n == 2*m);
  n = m; n -= n; TEST_ASSERT(IsZero(n));
  n = m; n *= n; TEST_ASSERT(n == m*m);
  n = m; n /= n; TEST_ASSERT(IsOne(n));
}




void program()
{
  GlobalManager CoCoAFoundations;

  test(NewZZmod(2));
  test(NewZZmod(-2));
  test(NewZZmod(3));
  test(NewZZmod(-3));
  test(NewZZmod(29641));
  test(NewZZmod(-29641));
  test(NewZZmod(10000019));
  test(NewZZmod(-10000019));
  test(NewZZmod(32768)); // not a field

  const ring ZZ = RingZZ();
  ideal I(RingElem(ZZ, 100000007));
  test(NewQuotientRing(ZZ, I));  // field

  I = ideal(power(RingElem(ZZ, 10), 20));
  test(NewQuotientRing(ZZ, I));  // not a field

  // Now check the three forbidden cases.

  // Should this really be forbidden?
  try { test(NewZZmod(0)); TEST_ASSERT(!"NEVER GET HERE!"); }
  catch (const CoCoA::ErrorInfo& err) { TEST_ASSERT(err == ERR::BadQuotRing); }

  try { test(NewZZmod(1)); TEST_ASSERT(!"NEVER GET HERE!"); }
  catch (const CoCoA::ErrorInfo& err) { TEST_ASSERT(err == ERR::BadQuotRing); }

  try { test(NewZZmod(-1)); TEST_ASSERT(!"NEVER GET HERE!"); }
  catch (const CoCoA::ErrorInfo& err) { TEST_ASSERT(err == ERR::BadQuotRing); }
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
