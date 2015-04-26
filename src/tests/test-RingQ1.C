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
#include "CoCoA/FractionField.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/IntOperations.H"
#include "CoCoA/error.H"
#include "CoCoA/ring.H"

using namespace CoCoA;

#include <iostream>
using std::cerr;
using std::endl;

#define TEST_ASSERT(cond)   do { if (!(cond)) CoCoA::AssertionFailed(#cond, __FILE__, __LINE__); } while(0)


void tautology(const RingElem& n)
{
  TEST_ASSERT(IsZero(n) ^ IsInvertible(n));
  BigInt N;
  TEST_ASSERT(IsInteger(N, n));
  TEST_ASSERT(n == N);
  TEST_ASSERT(N == n);
  RingElem n2(owner(n), N);
  TEST_ASSERT(n == n2);
  TEST_ASSERT(n2 == n);
}


void greater(const RingElem& n, int m)
{
  TEST_ASSERT(!(n == m));
  TEST_ASSERT(n != m);
  TEST_ASSERT(n >= m);
  TEST_ASSERT(n > m);
  TEST_ASSERT(!(n <= m));
  TEST_ASSERT(!(n < m));
  TEST_ASSERT(!(m == n));
  TEST_ASSERT(m != n);
  TEST_ASSERT(!(m >= n));
  TEST_ASSERT(!(m > n));
  TEST_ASSERT(m <= n);
  TEST_ASSERT(m < n);
}

void equal(const RingElem& n, int m)
{
  TEST_ASSERT(n == m);
  TEST_ASSERT(!(n != m));
  TEST_ASSERT(n >= m);
  TEST_ASSERT(!(n > m));
  TEST_ASSERT(n <= m);
  TEST_ASSERT(!(n < m));
  TEST_ASSERT(m == n);
  TEST_ASSERT(!(m != n));
  TEST_ASSERT(m >= n);
  TEST_ASSERT(!(m > n));
  TEST_ASSERT(m <= n);
  TEST_ASSERT(!(m < n));
}

void less(const RingElem& n, int m)
{
  TEST_ASSERT(!(n == m));
  TEST_ASSERT(n != m);
  TEST_ASSERT(!(n >= m));
  TEST_ASSERT(!(n > m));
  TEST_ASSERT(n <= m);
  TEST_ASSERT(n < m);
  TEST_ASSERT(!(m == n));
  TEST_ASSERT(m != n);
  TEST_ASSERT(m >= n);
  TEST_ASSERT(m > n);
  TEST_ASSERT(!(m <= n));
  TEST_ASSERT(!(m < n));
}


void greater(const RingElem& n, const BigInt& m)
{
  TEST_ASSERT(!(n == m));
  TEST_ASSERT(n != m);
  TEST_ASSERT(n >= m);
  TEST_ASSERT(n > m);
  TEST_ASSERT(!(n <= m));
  TEST_ASSERT(!(n < m));
  TEST_ASSERT(!(m == n));
  TEST_ASSERT(m != n);
  TEST_ASSERT(!(m >= n));
  TEST_ASSERT(!(m > n));
  TEST_ASSERT(m <= n);
  TEST_ASSERT(m < n);
}

void equal(const RingElem& n, const BigInt& m)
{
  TEST_ASSERT(n == m);
  TEST_ASSERT(!(n != m));
  TEST_ASSERT(n >= m);
  TEST_ASSERT(!(n > m));
  TEST_ASSERT(n <= m);
  TEST_ASSERT(!(n < m));
  TEST_ASSERT(m == n);
  TEST_ASSERT(!(m != n));
  TEST_ASSERT(m >= n);
  TEST_ASSERT(!(m > n));
  TEST_ASSERT(m <= n);
  TEST_ASSERT(!(m < n));
}

void less(const RingElem& n, const BigInt& m)
{
  TEST_ASSERT(!(n == m));
  TEST_ASSERT(n != m);
  TEST_ASSERT(!(n >= m));
  TEST_ASSERT(!(n > m));
  TEST_ASSERT(n <= m);
  TEST_ASSERT(n < m);
  TEST_ASSERT(!(m == n));
  TEST_ASSERT(m != n);
  TEST_ASSERT(m >= n);
  TEST_ASSERT(m > n);
  TEST_ASSERT(!(m <= n));
  TEST_ASSERT(!(m < n));
}


void arithmetic(RingElem n, const RingElem& m)
{
  BigInt N;
  TEST_ASSERT(IsInteger(N, n));
  BigInt M;
  TEST_ASSERT(IsInteger(M, m));
  TEST_ASSERT(n == N);
  TEST_ASSERT(m == M);

  // Some simple arithmetic
  n = n+m; N=N+M; TEST_ASSERT(n == N);
  n = n-m; N=N-M; TEST_ASSERT(n == N);
  n = n*m; N=N*M; TEST_ASSERT(n == N);
  n = n/m; N=N/M; TEST_ASSERT(n == N);

  n = n+M; N=N+M; TEST_ASSERT(n == N);
  n = n-M; N=N-M; TEST_ASSERT(n == N);
  n = n*M; N=N*M; TEST_ASSERT(n == N);
  n = n/M; N=N/M; TEST_ASSERT(n == N);

  n = N+m; N=N+M; TEST_ASSERT(n == N);
  n = N-m; N=N-M; TEST_ASSERT(n == N);
  n = N*m; N=N*M; TEST_ASSERT(n == N);
  n = N/m; N=N/M; TEST_ASSERT(n == N);

  // Checks for aliasing problems
  n = n+n; N=N+N; TEST_ASSERT(n == N);
  n = n*n; N=N*N; TEST_ASSERT(n == N);
  n = n/n; N=N/N; TEST_ASSERT(n == N);
  n = n-n; N=N-N; TEST_ASSERT(n == N);

  // assignment arithmetic
  N = 123456789;  n = N; TEST_ASSERT(n == N);
  n += m; N+=M; TEST_ASSERT(n == N);
  n -= m; N-=M; TEST_ASSERT(n == N);
  n *= m; N*=M; TEST_ASSERT(n == N);
  n /= m; N/=M; TEST_ASSERT(n == N);

  // Checks for aliasing problems
  n = m; n += n; N=M; N+=M; TEST_ASSERT(n == N);
  n = m; n -= n; N=M; N-=M; TEST_ASSERT(n == N);
  n = m; n *= n; N=M; N*=M; TEST_ASSERT(n == N);
  n = m; n /= n; N=M; N/=M; TEST_ASSERT(n == N);
}

void program()
{
  GlobalManager CoCoAFoundations;

  const ring QQ = RingQQ();

  TEST_ASSERT(IsFractionField(QQ));
  TEST_ASSERT(IsZZ(BaseRing(QQ)));
  TEST_ASSERT(BaseRing(QQ) == RingZZ());

  const BigInt Zero;
  RingElem n(QQ);
  // Check an implementation detail
  RingElem m = n; TEST_ASSERT(m == n); TEST_ASSERT(raw(m).myRawPtr() != raw(n).myRawPtr());

  TEST_ASSERT(IsZero(n));
  TEST_ASSERT(!IsInvertible(n));
  TEST_ASSERT(!IsOne(n));
  TEST_ASSERT(!IsMinusOne(n));
  equal(n, 0);
  equal(n, Zero);
  greater(n, -1);
  greater(n, BigInt(-1));
  less(n, 1);
  less(n, BigInt(1));

  n = 1;
  TEST_ASSERT(!IsZero(n));
  TEST_ASSERT(IsInvertible(n));
  TEST_ASSERT(IsOne(n));
  TEST_ASSERT(!IsMinusOne(n));
  equal(n, 1);
  equal(n, BigInt(1));
  greater(n, -1);
  greater(n, BigInt(-1));
  less(n, 2);
  less(n, BigInt(2));

  n = -1;
  TEST_ASSERT(!IsZero(n));
  TEST_ASSERT(IsInvertible(n));
  TEST_ASSERT(!IsOne(n));
  TEST_ASSERT(IsMinusOne(n));
  equal(n, -1);
  equal(n, BigInt(-1));
  greater(n, -2);
  greater(n, BigInt(-2));
  less(n, 1);
  less(n, BigInt(1));

  n = power(RingElem(QQ,10), 1000);
  TEST_ASSERT(IsInvertible(n));
  TEST_ASSERT(n == power(BigInt(10), 1000));

  tautology(zero(QQ));
  tautology(RingElem(QQ, 1));
  tautology(RingElem(QQ, -1));
  tautology(RingElem(QQ, 32003));
  tautology(RingElem(QQ, -32003));
  tautology(RingElem(QQ, power(BigInt(17),99)));
  tautology(RingElem(QQ, power(BigInt(-17),99)));

  // Check correct arithmetic.  Cases are small/large, positive/negative.
  arithmetic(RingElem(QQ, 23), RingElem(QQ, 37));
  arithmetic(RingElem(QQ, 100), RingElem(QQ, power(BigInt(10), 100)));
  arithmetic(RingElem(QQ, power(BigInt(10), 100)), RingElem(QQ, 100));
  arithmetic(RingElem(QQ, power(BigInt(10), 100)), RingElem(QQ, power(BigInt(2), 301)));

  arithmetic(RingElem(QQ, 23), -RingElem(QQ, 37));
  arithmetic(RingElem(QQ, 100), -RingElem(QQ, power(BigInt(10), 100)));
  arithmetic(RingElem(QQ, power(BigInt(10), 100)), -RingElem(QQ, 100));
  arithmetic(RingElem(QQ, power(BigInt(10), 100)), -RingElem(QQ, power(BigInt(2), 301)));

  arithmetic(-RingElem(QQ, 23), RingElem(QQ, 37));
  arithmetic(-RingElem(QQ, 100), RingElem(QQ, power(BigInt(10), 100)));
  arithmetic(-RingElem(QQ, power(BigInt(10), 100)), RingElem(QQ, 100));
  arithmetic(-RingElem(QQ, power(BigInt(10), 100)), RingElem(QQ, power(BigInt(2), 301)));

  arithmetic(-RingElem(QQ, 23), -RingElem(QQ, 37));
  arithmetic(-RingElem(QQ, 100), -RingElem(QQ, power(BigInt(10), 100)));
  arithmetic(-RingElem(QQ, power(BigInt(10), 100)), -RingElem(QQ, 100));
  arithmetic(-RingElem(QQ, power(BigInt(10), 100)), -RingElem(QQ, power(BigInt(2), 301)));

  // assignment arithmetic
  n = 355; n /= 113; // Aside: 355/113 is a famous approximation for pi.
  m = n;

  n += m;
  n -= m;
  n *= m;
  n /= m;

  n = m; n += n; TEST_ASSERT(n == 2*m);
  n = m; n -= n; TEST_ASSERT(IsZero(n));
  n = m; n *= n; TEST_ASSERT(n == m*m);
  n = m; n /= n; TEST_ASSERT(IsOne(n));

  // compute a large number
  n = power(m, 100);
  TEST_ASSERT(IsInvertible(n));
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
    cerr << "***ERROR***  UNCAUGHT CoCoA::Error";
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
