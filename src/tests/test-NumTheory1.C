//   Copyright (c)  2009-2010  John Abbott

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
#include "CoCoA/error.H"
#include "CoCoA/VectorOperations.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/NumTheory.H"
#include "CoCoA/IntOperations.H"

#include <cmath>
using std::log;
#include <cstdlib>
using std::abs;
#include <iomanip>
using std::flush;
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <limits>
using std::numeric_limits;
#include <vector>
using std::vector;

using namespace CoCoA;

#define TEST_ASSERT(cond)   do { if (!(cond)) CoCoA::AssertionFailed(#cond, __FILE__, __LINE__); } while(0)




void program()
{
  GlobalManager CoCoAFoundations(UseGMPAllocator);

  cout << "Testing gcd and lcm on small numbers." << endl;
  const int Nmax = 257;
  for (int a=-Nmax; a < Nmax; ++a)
    for (int b=-Nmax; b < Nmax; ++b)
    {
      const BigInt A(a);
      const BigInt B(b);
      const int gcdab = gcd(a,b);
      TEST_ASSERT(gcdab >= 0);
      TEST_ASSERT((a == 0 && b == 0) || gcdab > 0);
      TEST_ASSERT(gcdab == 0 || (a%gcdab == 0 && b%gcdab == 0));
      TEST_ASSERT(gcdab == gcd(a,B));
      TEST_ASSERT(gcdab == gcd(A,b));
      TEST_ASSERT(gcdab == gcd(A,B));

      const int lcmab = lcm(a,b);
      TEST_ASSERT(gcdab*lcmab == abs(a*b));
      TEST_ASSERT(lcmab == lcm(a,B));
      TEST_ASSERT(lcmab == lcm(A,b));
      TEST_ASSERT(lcmab == lcm(A,B));
    }

  cout << "\nTesting ExtGcd on small numbers" << endl;
  for (int a=-Nmax; a < Nmax; ++a)
    for (int b=-Nmax; b < Nmax; ++b)
    {
      if (a == 0 && b == 0) continue;
      long CofacA, CofacB;
      const int g = ExtGcd(CofacA,CofacB,a,b);
      const int g2 = gcd(a,b);
      TEST_ASSERT(g == g2);
      TEST_ASSERT(CofacA*a+CofacB*b == g);
      if (abs(a) == abs(b))
      {
        TEST_ASSERT((CofacA == 0 && CofacB == sign(b)) || (CofacB == 0 && CofacA == sign(a)));
        continue;
      }
      if (abs(a) == g) { TEST_ASSERT(CofacA == sign(a) && CofacB == 0); continue; }
      if (abs(b) == g) { TEST_ASSERT(CofacB == sign(b) && CofacA == 0); continue; }
      // General case
      TEST_ASSERT(abs(CofacA) <= (abs(b)/g)/2  && abs(CofacB) <= (abs(a)/g)/2);
      if (2*g == abs(b))
        TEST_ASSERT(abs(CofacA) == 1);
      else
        TEST_ASSERT(2*abs(CofacA) < abs(b)/g);
      if (2*g == abs(a))
        TEST_ASSERT(abs(CofacB) == 1);
      else
        TEST_ASSERT(2*abs(CofacB) < abs(a)/g);
    }


  cout << "\nTesting ExtGcd on BigInt numbers" << endl;
  for (int a=-Nmax; a < Nmax; ++a)
    for (int b=-Nmax; b < Nmax; ++b)
    {
      if (a == 0 && b == 0) continue;
      BigInt A(a);
      BigInt B(b);
      BigInt CofacA, CofacB;
      const BigInt G = ExtGcd(CofacA,CofacB,A,B);
      const int g = gcd(a,b);
      TEST_ASSERT(G == g);
      TEST_ASSERT(CofacA*a+CofacB*b == g);
      if (abs(a) == abs(b))
      {
        TEST_ASSERT((CofacA == 0 && CofacB == sign(b)) || (CofacB == 0 && CofacA == sign(a)));
        continue;
      }
      if (abs(a) == g) { TEST_ASSERT(CofacA == sign(a) && CofacB == 0); continue; }
      if (abs(b) == g) { TEST_ASSERT(CofacB == sign(b) && CofacA == 0); continue; }
      // General case
      TEST_ASSERT(abs(CofacA) <= (abs(b)/g)/2  && abs(CofacB) <= (abs(a)/g)/2);
      if (2*g == abs(b))
        TEST_ASSERT(abs(CofacA) == 1);
      else
        TEST_ASSERT(2*abs(CofacA) < abs(b)/g);
      if (2*g == abs(a))
        TEST_ASSERT(abs(CofacB) == 1);
      else
        TEST_ASSERT(2*abs(CofacB) < abs(a)/g);
    }


  cout << "\nTesting InvMod" << endl;
  for (int r=-Nmax; r < Nmax; ++r)
    for (int m=2; m < Nmax; ++m)
    {
      BigInt R(r);
      BigInt M(m);
      const int invr = InvMod(r,m);
      TEST_ASSERT(gcd(r,m)==1 || invr==0);
      if (invr==0) continue;
      if (r > 0) TEST_ASSERT((r*invr)%m == 1);
      if (r < 0) TEST_ASSERT(((-r)*invr)%m == m-1);
      TEST_ASSERT(invr == InvMod(R,m));
      TEST_ASSERT(invr == InvMod(r,M));
      TEST_ASSERT(invr == InvMod(R,M));
    }

  cout << "\nTesting PowerMod" << endl;
  for (int r=-Nmax; r < Nmax; ++r)
    for (int m=2; m < Nmax; ++m)
    {
      BigInt R(r);
      BigInt M(m);
      TEST_ASSERT(PowerMod(r,0,m) == 1);
      TEST_ASSERT(PowerMod(r,0,M) == 1);
      TEST_ASSERT(PowerMod(R,0,m) == 1);
      TEST_ASSERT(PowerMod(R,0,M) == 1);

      long RmodM;
      if (r >= 0) RmodM = r%m; else RmodM = (m-((-r)%m))%m;
      TEST_ASSERT(PowerMod(r,1,m) == RmodM);
      TEST_ASSERT(PowerMod(r,1,M) == RmodM);
      TEST_ASSERT(PowerMod(R,1,m) == RmodM);
      TEST_ASSERT(PowerMod(R,1,M) == RmodM);

      if (gcd(r,m) == 1)
      {
        long InvR = InvMod(r,m);
        TEST_ASSERT(PowerMod(r,-1,m) == InvR);
        TEST_ASSERT(PowerMod(r,-1,M) == InvR);
        TEST_ASSERT(PowerMod(R,-1,m) == InvR);
        TEST_ASSERT(PowerMod(R,-1,M) == InvR);
      }
      for (int e=0; e < 9; ++e)
      {
        TEST_ASSERT(PowerMod(PowerMod(r,e,m), 2, m) == PowerMod(r, 2*e, m));
        if (gcd(r,m) == 1)
          TEST_ASSERT((PowerMod(r,e,m)*PowerMod(r,-e,m))%m == 1);
      }
    }



  cout << "\nTesting IsPrime for small numbers (incl negative ones)." << endl;
  int PrimeCounter = 0;
  for (int n=1; n <= 1000; ++n)
  {
    if (IsPrime(n)) ++PrimeCounter;
    TEST_ASSERT(IsPrime(n) == IsPrime(BigInt(n)));
  }
  TEST_ASSERT(PrimeCounter == 168);

  // ATTN!  this next test depends on the platform (32-bit or 64-bit)
  const int LongBits = numeric_limits<long>::digits;
  cout << "\nTesting IsPrime for numbers just greater than max long." << endl;
  const BigInt StartingPoint = power(2,LongBits);
  for (int delta=1; delta < 99; delta+=2)
  {
    const BigInt N = StartingPoint+delta;
    if (!IsProbPrime(N)) continue;
    if (IsPrime(N)) continue;
    cout << "Surprise: " << N << " is a composite pseudo-prime." << endl;
  }


  cout << "\nTesting IsProbPrime: lengths of repunit (probable) primes (up to 500):" << flush;
  for (int n=1; n <= 500; ++n)
  {
    if (!IsPrime(n)) continue; // no need to try composite n
    if (IsProbPrime((power(10,n)-1)/9))
      cout << " " << n << flush;
  }
  cout << endl;


  cout << "\nTesting NextPrime & PrevPrime" << endl;
  const long MaxLong = numeric_limits<long>::max();
  TEST_ASSERT(NextPrime(1) == 2);
  TEST_ASSERT(NextPrime(2) == 3);
  TEST_ASSERT(NextPrime(3) == 5);
  TEST_ASSERT(NextPrime(4) == 5);

  TEST_ASSERT(NextPrime(MaxLong) == 0);
  for (long n=MaxLong; !IsPrime(n--); )
    TEST_ASSERT(NextPrime(n) == 0);

  const long MaxPrime = IsPrime(MaxLong)?MaxLong:PrevPrime(MaxLong);
  TEST_ASSERT(IsPrime(MaxPrime));
  TEST_ASSERT(NextPrime(MaxPrime) == 0);
  for (long n=MaxPrime; !IsPrime(n--); )
    TEST_ASSERT(NextPrime(n) == MaxPrime);

  TEST_ASSERT(PrevPrime(4) == 3);
  TEST_ASSERT(PrevPrime(3) == 2);
  TEST_ASSERT(PrevPrime(2) == 0);
  TEST_ASSERT(PrevPrime(1) == 0);
  for (int n=3; n < 1000; n += 2)
  {
    if (IsPrime(n))
    {
      TEST_ASSERT(PrevPrime(NextPrime(n)) == n);
      TEST_ASSERT(NextPrime(PrevPrime(n)) == n);
    }
    else
    {
      TEST_ASSERT(NextPrime(PrevPrime(n)) == NextPrime(n));
      TEST_ASSERT(PrevPrime(NextPrime(n)) == PrevPrime(n));
    }
  }

  cout << "\nTesting PrimitiveRoot" << endl;
  const int MaxHist = 20;
  vector<int> histogram(MaxHist);
  int p = 1;
  while (p < 65535)
  {
    p = NextPrime(p);
    const int g = PrimitiveRoot(p);
    if (g < MaxHist) ++histogram[g];
//     double ratio = g/(std::log(p)*std::log(std::log(p)));
//     if (ratio > 2.0) cout << "ratio=" << ratio << "  for p=" << p << " and g=" << g << endl;
  }
  cout << "Histogram of least positive primitive roots: " << histogram << endl;


  cout << "\nTesting factorize on smaller & larger integers:" << endl;
  for (int i=1; i <= 16; ++i)
  {
    const BigInt N = power(10,i)-1;
    cout << "Factorization of " << N << ": " << factor(N) << endl;
  }

  cout << "\nTesting EulerPhi on small numbers." << endl;
  cout << "Here are the most frequent values of EulerPhi on numbers up to 100000:\n"
                 << "Value  Freq\n";
  vector<int> freq(100000);
  for (int i=2; i < 100000; ++i)
    ++freq[EulerPhi(i)];
  for (int i=2; i < 100000; ++i)
    if (freq[i] > 200) 
      cout << i << "  " << freq[i] << endl;

  cout << "\nTesting EulerPhi on larger numbers." << endl
                 << "Looking for smallest number with EulerPhi(N) < N/10." << endl;
  BigInt BigNum;
  BigNum = 2;
  int LastPrime = -1;
  for (int p=3; 10*EulerPhi(BigNum) > BigNum; p=NextPrime(p))
  {
    BigNum *= p;
    LastPrime = p;
  }
  BigInt PhiBigNum = EulerPhi(BigNum);
  cout << "Result is  N=" << BigNum << endl
                 << "Its EulerPhi=" << PhiBigNum << endl
                 << "In fact N is just the product of all primes up to " << LastPrime << endl;


  // SmoothFactor
  cout << "\nTesting SmoothFactor" << endl;
  factorization<BigInt> FactoredFactorial = SmoothFactor(factorial(101), 101);
  cout << "101-smooth factors of factorial of 101: " << FactoredFactorial << endl;
  FactoredFactorial = SmoothFactor(factorial(101), 50);
  cout << "50-smooth factors of factorial of 101: " << FactoredFactorial << endl;
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
