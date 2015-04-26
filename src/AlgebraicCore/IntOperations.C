//   Copyright (c)  2012  John Abbott

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
#include "CoCoA/convert.H"
#include "CoCoA/error.H"
#include "CoCoA/utils.H"

#include <cmath>
//using std::log;
//using std::atan;
#include <limits>
using std::numeric_limits;



namespace CoCoA
{

  namespace //anonymous
  {
    // Compute  round(n/d)  rounding halves up.
    // NB cannot simply do  (n+d/2)/d  because of possible overflow!
    unsigned long uround_half_up(unsigned long n, unsigned long d)
    {
      CoCoA_ASSERT(d != 0);
      const unsigned long q = n/d;
      if (n-q*d <= (d-1)/2) return q;
      return q+1; // cannot overflow
    }

    // Compute  round(n/d)  rounding halves down.
    // NB cannot simply do  (n+(d-1)/2)/d  because of possible overflow!
    unsigned long uround_half_down(unsigned long n, unsigned long d)
    {
      CoCoA_ASSERT(d != 0);
      const unsigned long q = n/d;
      if (n-q*d <= d/2) return q;
      return q+1; // cannot overflow
    }

    // These two defns are used inside cmp(MachineInt, MachineInt)
    inline int cmp(long a, long b)
    {
      if (a < b) return -1;
      if (b < a) return 1;
      return 0;
    }

    inline int cmp(unsigned long a, unsigned long b)
    {
      if (a < b) return -1;
      if (b < a) return 1;
      return 0;
    }

  } // end of anonymous namespace


  const BigInt abs(const BigInt& N)
  {
    BigInt ans;
    mpz_abs(mpzref(ans), mpzref(N));
    return ans;
  }


  const BigInt operator-(const BigInt& N)
  {
    BigInt ans(N);
    negate(ans);
    return ans;
  }


  const BigInt operator+(const BigInt& N1, const BigInt& N2)
  {
    BigInt ans;
    mpz_add(mpzref(ans), mpzref(N1), mpzref(N2));
    return ans;
  }

  const BigInt operator+(const BigInt& N1, const MachineInt& n2)
  {
    BigInt ans;
    if (IsNegative(n2))
      mpz_sub_ui(mpzref(ans), mpzref(N1), abs(n2));
    else
      mpz_add_ui(mpzref(ans), mpzref(N1), AsUnsignedLong(n2));
    return ans;
  }

  const BigInt operator+(const MachineInt& n1, const BigInt& N2)
  {
    BigInt ans;
    if (IsNegative(n1))
      mpz_sub_ui(mpzref(ans), mpzref(N2), abs(n1));
    else
      mpz_add_ui(mpzref(ans), mpzref(N2), AsUnsignedLong(n1));
    return ans;
  }


  const BigInt operator-(const BigInt& N1, const BigInt& N2)
  {
    BigInt ans;
    mpz_sub(mpzref(ans), mpzref(N1), mpzref(N2));
    return ans;
  }

  const BigInt operator-(const BigInt& N1, const MachineInt& n2)
  {
    BigInt ans;
    if (IsNegative(n2))
      mpz_add_ui(mpzref(ans), mpzref(N1), abs(n2));
    else
      mpz_sub_ui(mpzref(ans), mpzref(N1), AsUnsignedLong(n2));
    return ans;
  }

  const BigInt operator-(const MachineInt& n1, const BigInt& N2)
  {
    BigInt ans;
    if (IsNegative(n1))
      mpz_add_ui(mpzref(ans), mpzref(N2), abs(n1));
    else
      mpz_sub_ui(mpzref(ans), mpzref(N2), AsUnsignedLong(n1));
    negate(ans);
    return ans;
  }


  const BigInt operator*(const BigInt& N1, const BigInt& N2)
  {
    BigInt ans;
    mpz_mul(mpzref(ans), mpzref(N1), mpzref(N2));
    return ans;
  }

  const BigInt operator*(const BigInt& N1, const MachineInt& n2)
  {
    BigInt ans;
    if (IsNegative(n2))
      mpz_mul_si(mpzref(ans), mpzref(N1), AsSignedLong(n2));
    else
      mpz_mul_ui(mpzref(ans), mpzref(N1), AsUnsignedLong(n2));
    return ans;
  }

  const BigInt operator*(const MachineInt& n1, const BigInt& N2)
  {
    BigInt ans;
    if (IsNegative(n1))
      mpz_mul_si(mpzref(ans), mpzref(N2), AsSignedLong(n1));
    else
      mpz_mul_ui(mpzref(ans), mpzref(N2), AsUnsignedLong(n1));
    return ans;
  }


  const BigInt operator/(const BigInt& N1, const BigInt& N2)
  {
    if (IsZero(N2))
      CoCoA_ERROR(ERR::DivByZero, "BigInt / BigInt");
    BigInt ans;
    mpz_tdiv_q(mpzref(ans), mpzref(N1), mpzref(N2));
    return ans;
  }

  const BigInt operator/(const BigInt& N1, const MachineInt& n2)
  {
    if (IsZero(n2))
      CoCoA_ERROR(ERR::DivByZero, "BigInt / MachineInt");
    BigInt ans;
    if (IsNegative(n2))
    {
      mpz_tdiv_q_ui(mpzref(ans), mpzref(N1), abs(n2));
      negate(ans);
    }
    else
      mpz_tdiv_q_ui(mpzref(ans), mpzref(N1), AsUnsignedLong(n2));
    return ans;
  }

  // Impl is simple rather than fast.
  const BigInt operator/(const MachineInt& n1, const BigInt& N2)
  {
    if (IsZero(N2))
      CoCoA_ERROR(ERR::DivByZero, "MachineInt / BigInt");
    BigInt ans(n1);
    ans /= N2;
    return ans;
  }


  const BigInt operator%(const BigInt& N1, const BigInt& N2)
  {
    if (IsZero(N2))
      CoCoA_ERROR(ERR::ZeroModulus, "BigInt % BigInt");
    BigInt ans;
    mpz_tdiv_r(mpzref(ans), mpzref(N1), mpzref(N2));
    return ans;
  }

  // This impl is valid for any second arg, but would then
  // have to return an unsigned long.  I have decided to forbid
  // 2nd args which are too big to fit into signed long; this
  // guarantees that the result will fit into signed long.
  // An alternative would be simply to NumericCast the result to
  // long (which will throw if the result doesn't fit).
  long operator%(const BigInt& N1, const MachineInt& n2)
  {
    if (IsZero(n2))
      CoCoA_ERROR(ERR::DivByZero, "BigInt % MachineInt");
    if (!IsSignedLong(n2))
      CoCoA_ERROR(ERR::ArgTooBig, "BigInt % MachineInt");
    const long ans = mpz_tdiv_ui(mpzref(N1), abs(n2)); // abs value of remainder
    if (N1 < 0) return -ans; // CAREFUL: read GMP doc for mpz_tdiv_ui
    return ans;
  }

  // Impl is simple rather than fast.
  const BigInt operator%(const MachineInt& n1, const BigInt& N2)
  {
    if (IsZero(N2))
      CoCoA_ERROR(ERR::DivByZero, "MachineInt % BigInt");
    return BigInt(n1)%N2;
//    BigInt ans(n1);
//    ans %= N2;
//    return ans;
  }


  const BigInt LeastNNegRemainder(const BigInt& N1, const BigInt& N2)
  {
    if (IsZero(N2))
      CoCoA_ERROR(ERR::ZeroModulus, "LeastNNegRemainder(BigInt, BigInt)");
    BigInt ans;
    mpz_mod(mpzref(ans), mpzref(N1), mpzref(N2));
    return ans;
  }

  const BigInt LeastNNegRemainder(const BigInt& N1, const MachineInt& n2)
  {
    return BigInt(mpz_fdiv_ui(mpzref(N1), abs(n2)));
  }

  const BigInt LeastNNegRemainder(const MachineInt& n1, const BigInt& N2)
  {
    // Simple rather than fast;
    return LeastNNegRemainder(BigInt(n1), N2);
  }

  const BigInt LeastNNegRemainder(const MachineInt& n1, const MachineInt& n2)
  {
    const unsigned long N2 = abs(n2);
    const unsigned long ans = abs(n1)%N2;
    if (ans != 0 && IsNegative(n1)) return BigInt(N2-ans);
    return BigInt(ans);
  }


  const BigInt SymmRemainder(const BigInt& N1, const BigInt& N2)
  {
    if (IsZero(N2))
      CoCoA_ERROR(ERR::ZeroModulus, "SymmRemainder(BigInt, BigInt)");
    BigInt ans;
    mpz_tdiv_r(mpzref(ans), mpzref(N1), mpzref(N2));
    BigInt HalfN2; mpz_tdiv_q_2exp(mpzref(HalfN2), mpzref(N2), 1); mpz_abs(mpzref(HalfN2), mpzref(HalfN2));
    if (CmpAbs(ans, HalfN2) < 0) return ans;
    if (ans == HalfN2) return ans;
    if (sign(ans) == sign(N2)) ans -= N2;
    else ans += N2;
    return ans;
  }

  long SymmRemainder(const BigInt& N1, const MachineInt& n2)
  {
    const unsigned long N2 = abs(n2);
    unsigned long ans = mpz_fdiv_ui(mpzref(N1), N2);
    if (ans > N2/2) return -NumericCast<long>(N2-ans); // cannot overflow (on a binary computer)
    return ans;
  }

  const BigInt SymmRemainder(const MachineInt& n1, const BigInt& N2)
  {
    // Simple rather than fast.
    return SymmRemainder(BigInt(n1), N2);
  }


  long SymmRemainder(const MachineInt& n1, const MachineInt& n2)
  {
    const unsigned long N2 = abs(n2);
    unsigned long ans = abs(n1)%N2;
    if (IsNegative(n1)) ans = N2-ans;
    if (ans > N2/2) return -NumericCast<long>(N2-ans); // cannot overflow (on a binary computer)
    return ans;
  }


  const BigInt power(const BigInt& base, const BigInt& exponent)
  {
    if (exponent < 0)
      CoCoA_ERROR(ERR::NegExp, "power(BigInt, BigInt)");
    // BUG??? Could allow high powers of 0 and +/-1.  Is it worth it?
    unsigned long exp;
    if (!IsConvertible(exp, exponent))
      CoCoA_ERROR(ERR::ExpTooBig, "power(BigInt, BigInt)");
    BigInt ans;
    mpz_pow_ui(mpzref(ans), mpzref(base), exp);
    return ans;
  }

  const BigInt power(const BigInt& base, const MachineInt& exponent)
  {
    if (IsNegative(exponent))
      CoCoA_ERROR(ERR::NegExp, "power(BigInt, MachineInt)");
    BigInt ans;
    mpz_pow_ui(mpzref(ans), mpzref(base), AsUnsignedLong(exponent));
    return ans;
  }

  const BigInt power(const MachineInt& base, const BigInt& exponent)
  {
    const char* const FnName = "power(MachineInt, BigInt)";
    if (exponent < 0)
      CoCoA_ERROR(ERR::NegExp, FnName);
    // BUG???  Could allow high powers of 0 and +/-1.  Is it worth it?
    unsigned long exp;
    if (!IsConvertible(exp, exponent))
      CoCoA_ERROR(ERR::ExpTooBig, FnName);
    const bool NegateResult = IsNegative(base) & IsOdd(exp);
    const unsigned long B = abs(base);
    BigInt ans;
    mpz_ui_pow_ui(mpzref(ans), B, exp);
    if (NegateResult)
      negate(ans);
    return ans;
  }

  const BigInt power(const MachineInt& base, const MachineInt& exponent)
  {
    const char* const FnName = "power(MachineInt, MachineInt)";
    if (IsNegative(exponent))
      CoCoA_ERROR(ERR::NegExp, FnName);

    const bool NegateResult = IsNegative(base) & IsOdd(AsUnsignedLong(exponent));
    const unsigned long B = abs(base);
    BigInt ans;
    mpz_ui_pow_ui(mpzref(ans), B, AsUnsignedLong(exponent));
    if (NegateResult)
      negate(ans);
    return ans;
  }

  long SmallPower(const MachineInt& base, const MachineInt& exponent)
  {
    if (IsNegative(exponent)) CoCoA_ERROR(ERR::NegExp, "SmallPower");
    if (!IsSignedLong(base)) CoCoA_ERROR(ERR::BadArg, "SmallPower: base is too big to fit into a long");
    if (IsZero(exponent)) return 1; // Order is important for  0^0
    if (IsZero(base)) return 0;     // ...
    unsigned long e = abs(exponent);
    long b = abs(base);
    const bool negative = (e&1) && IsNegative(base);
    long ans = 1;
    while (e != 1)
    {
      if (e&1) ans *= b;
      e >>= 1;
      b *= b;
    }
    ans *= b;
    if (negative) return -ans;
    return ans;
  }


  void add(BigInt& lhs, const BigInt& N1, const BigInt& N2)
  {
    mpz_add(mpzref(lhs), mpzref(N1), mpzref(N2));
  }

  void add(BigInt& lhs, const BigInt& N1, const MachineInt& n2)
  {
    if (IsNegative(n2))
      mpz_sub_ui(mpzref(lhs), mpzref(N1), abs(n2));
    else
      mpz_add_ui(mpzref(lhs), mpzref(N1), AsUnsignedLong(n2));
  }

  void add(BigInt& lhs, const MachineInt& n1, const BigInt& N2)
  {
    add(lhs, N2, n1);
  }


  void sub(BigInt& lhs, const BigInt& N1, const BigInt& N2)
  {
    mpz_sub(mpzref(lhs), mpzref(N1), mpzref(N2));
  }

  void sub(BigInt& lhs, const BigInt& N1, const MachineInt& n2)
  {
    if (IsNegative(n2))
      mpz_add_ui(mpzref(lhs), mpzref(N1), abs(n2));
    else
      mpz_sub_ui(mpzref(lhs), mpzref(N1), AsUnsignedLong(n2));
  }

  void sub(BigInt& lhs, const MachineInt& n1, const BigInt& N2)
  {
    sub(lhs, N2, n1);
    negate(lhs);
  }


  void mul(BigInt& lhs, const BigInt& N1, const BigInt& N2)
  {
    mpz_mul(mpzref(lhs), mpzref(N1), mpzref(N2));
  }

  void mul(BigInt& lhs, const BigInt& N1, const MachineInt& n2)
  {
    if (IsNegative(n2))
      mpz_mul_si(mpzref(lhs), mpzref(N1), AsSignedLong(n2));
    else
      mpz_mul_ui(mpzref(lhs), mpzref(N1), AsUnsignedLong(n2));
  }

  void mul(BigInt& lhs, const MachineInt& n1, const BigInt& N2)
  {
    mul(lhs, N2, n1);
  }


  void div(BigInt& lhs, const BigInt& N1, const BigInt& N2)
  {
    if (IsZero(N2))
      CoCoA_ERROR(ERR::DivByZero, "div(BigInt, BigInt, BigInt)");
    mpz_tdiv_q(mpzref(lhs), mpzref(N1), mpzref(N2));
  }

  void div(BigInt& lhs, const BigInt& N1, const MachineInt& n2)
  {
    if (IsZero(n2))
      CoCoA_ERROR(ERR::DivByZero, "div(BigInt, BigInt, MachineInt)");
    if (IsNegative(n2))
    {
      mpz_tdiv_q_ui(mpzref(lhs), mpzref(N1), abs(n2));
      negate(lhs);
    }
    else
      mpz_tdiv_q_ui(mpzref(lhs), mpzref(N1), AsUnsignedLong(n2));
  }

  // Impl is simple rather than fast.
  void div(BigInt& lhs, const MachineInt& n1, const BigInt& N2)
  {
    if (N2 == 0)
      CoCoA_ERROR(ERR::DivByZero, "div(BigInt, MachineInt, BigInt)");
    div(lhs, BigInt(n1), N2);
  }


  void mod(BigInt& lhs, const BigInt& N1, const BigInt& N2)
  {
    if (IsZero(N2))
      CoCoA_ERROR(ERR::ZeroModulus, "mod(BigInt, BigInt, BigInt)");
    mpz_tdiv_r(mpzref(lhs), mpzref(N1), mpzref(N2));
  }

  void mod(BigInt& lhs, const BigInt& N1, const MachineInt& n2)
  {
    if (IsZero(n2))
      CoCoA_ERROR(ERR::ZeroModulus, "mod(BigInt, BigInt, MachineInt)");
    mpz_tdiv_r_ui(mpzref(lhs), mpzref(N1), abs(n2));
  }

  // Impl is simple rather than fast.
  void mod(BigInt& lhs, const MachineInt& n1, const BigInt& N2)
  {
    if (IsZero(N2))
      CoCoA_ERROR(ERR::ZeroModulus, "mod(BigInt, MachineInt, BigInt)");
    mod(lhs, BigInt(n1), N2);
  }


  void quorem(BigInt& quo, BigInt& rem, const BigInt& num, const BigInt& den)
  {
    if (IsZero(den))
      CoCoA_ERROR(ERR::DivByZero, "quorem(BigInt, BigInt, BigInt, BigInt)");
    mpz_tdiv_qr(mpzref(quo), mpzref(rem), mpzref(num), mpzref(den));
  }

  void quorem(BigInt& quo, BigInt& rem, const BigInt& num, const MachineInt& den)
  {
    if (IsZero(den))
      CoCoA_ERROR(ERR::DivByZero, "quorem(BigInt, BigInt, BigInt, MachineInt)");
    if (IsNegative(den))
    {
      mpz_tdiv_qr_ui(mpzref(quo), mpzref(rem), mpzref(num), abs(den));
      negate(quo);
    }
    else
      mpz_tdiv_qr_ui(mpzref(quo), mpzref(rem), mpzref(num), AsUnsignedLong(den));
  }

  // Impl is simple rather than fast.
  void quorem(BigInt& quo, BigInt& rem, const MachineInt& num, const BigInt& den)
  {
    if (IsZero(den))
      CoCoA_ERROR(ERR::DivByZero, "quorem(BigInt, BigInt, MachineInt, BigInt)");
    quorem(quo, rem, BigInt(num), den);
  }


  void divexact(BigInt& lhs, const BigInt& N1, const BigInt& N2)
  {
    CoCoA_ASSERT(N1%N2 == 0);
    mpz_divexact(mpzref(lhs), mpzref(N1), mpzref(N2));
  }

//???  void divexact(BigInt& lhs, const MachineInt& N1, const BigInt& N2);
//???  void divexact(BigInt& lhs, const BigInt& N1, const MachineInt& N2);


  bool IsZero(const BigInt& N)
  {
    return mpz_sgn(mpzref(N))==0;
  }

  bool IsOne(const BigInt& N)
  {
    return mpz_cmp_ui(mpzref(N), 1)==0;
  }

  bool IsMinusOne(const BigInt& N)
  {
    return mpz_cmp_si(mpzref(N),-1)==0;
  }


  bool IsPowerOf2(const BigInt& N)
  {
    if (N <= 0) return false;
    const unsigned long LogN = mpz_sizeinbase(mpzref(N), 2)-1;
    const unsigned long Last1Bit = mpz_scan1(mpzref(N), 0);
    return (Last1Bit == LogN);
  }

  bool IsPowerOf2(const MachineInt& n)
  {
    if (IsZero(n) || IsNegative(n)) return false;
    const unsigned long N = AsUnsignedLong(n);
    return !(N & (N-1));
  }


  bool IsDivisible(const MachineInt& N, const MachineInt& D) // is N divisibile by D?
  {
    if (IsZero(D)) CoCoA_ERROR(ERR::DivByZero, "IsDivisible");
    return abs(N)%abs(D) == 0;
  }

  bool IsDivisible(const MachineInt& N, const BigInt& D) // is N divisibile by D?
  {
    if (IsZero(D)) CoCoA_ERROR(ERR::DivByZero, "IsDivisible");
    return IsDivisible(BigInt(N), D);
  }

  bool IsDivisible(const BigInt& N, const MachineInt& D) // is N divisibile by D?
  {
    if (IsZero(D)) CoCoA_ERROR(ERR::DivByZero, "IsDivisible");
    return mpz_divisible_ui_p(mpzref(N), abs(D));
  }

  bool IsDivisible(const BigInt& N, const BigInt& D) // is N divisibile by D?
  {
    if (IsZero(D)) CoCoA_ERROR(ERR::DivByZero, "IsDivisible");
    return mpz_divisible_p(mpzref(N), mpzref(D));
  }


  int cmp(const BigInt& N1, const BigInt& N2)
  {
    return sign(mpz_cmp(mpzref(N1), mpzref(N2)));
  }

  int cmp(const BigInt& N1, const MachineInt& n2)
  {
    if (IsNegative(n2))
      return sign(mpz_cmp_si(mpzref(N1), AsSignedLong(n2)));
    else
      return sign(mpz_cmp_ui(mpzref(N1), AsUnsignedLong(n2)));
  }

  int cmp(const MachineInt& n1, const BigInt& N2)
  {
    return -cmp(N2, n1);
  }

  int cmp(const MachineInt& n1, const MachineInt& n2)
  {
    if (IsNegative(n1))
    {
      if (IsNegative(n2))
        return cmp(AsSignedLong(n1), AsSignedLong(n2));
      else
        return -1;
    }
    // n1 is non-negative
    if (IsNegative(n2)) return 1;
    // Both n1 and n2 are non-negative
    return cmp(AsUnsignedLong(n1), AsUnsignedLong(n2));
  }


  int CmpAbs(const BigInt& N1, const BigInt& N2)
  {
    return sign(mpz_cmpabs(mpzref(N1), mpzref(N2)));
  }

  int CmpAbs(const BigInt& N1, const MachineInt& n2)
  {
    return sign(mpz_cmpabs_ui(mpzref(N1), abs(n2)));
  }

  int CmpAbs(const MachineInt& n1, const BigInt& N2)
  {
    return -sign(mpz_cmpabs_ui(mpzref(N2), abs(n1)));
  }

  int CmpAbs(const MachineInt& n1, const MachineInt& n2)
  {
    const unsigned long a1 = abs(n1);
    const unsigned long a2 = abs(n2);
    if (a1 > a2) return 1;
    if (a1 < a2) return -1;
    return 0;
  }



  bool operator==(const BigInt& N1, const BigInt& N2)
  {
    return cmp(N1, N2) == 0;
  }

  bool operator==(const BigInt& N1, const MachineInt& n2)
  {
    return cmp(N1, n2) == 0;
  }

  bool operator==(const MachineInt& n1, const BigInt& N2)
  {
    return cmp(n1, N2) == 0;
  }


  bool operator!=(const BigInt& N1, const BigInt& N2)
  {
    return cmp(N1, N2) != 0;
  }

  bool operator!=(const BigInt& N1, const MachineInt& n2)
  {
    return cmp(N1, n2) != 0;
  }

  bool operator!=(const MachineInt& n1, const BigInt& N2)
  {
    return cmp(n1, N2) != 0;
  }


  bool operator> (const BigInt& N1, const BigInt& N2)
  {
    return cmp(N1, N2) > 0;
  }

  bool operator> (const BigInt& N1, const MachineInt& n2)
  {
    return cmp(N1, n2) > 0;
  }

  bool operator> (const MachineInt& n1, const BigInt& N2)
  {
    return cmp(n1, N2) > 0;
  }


  bool operator>=(const BigInt& N1, const BigInt& N2)
  {
    return cmp(N1, N2) >= 0;
  }

  bool operator>=(const BigInt& N1, const MachineInt& n2)
  {
    return cmp(N1, n2) >= 0;
  }

  bool operator>=(const MachineInt& n1, const BigInt& N2)
  {
    return cmp(n1, N2) >= 0;
  }


  bool operator< (const BigInt& N1, const BigInt& N2)
  {
    return cmp(N1, N2) < 0;
  }

  bool operator< (const BigInt& N1, const MachineInt& n2)
  {
    return cmp(N1, n2) < 0;
  }

  bool operator< (const MachineInt& n1, const BigInt& N2)
  {
    return cmp(n1, N2) < 0;
  }


  bool operator<=(const BigInt& N1, const BigInt& N2)
  {
    return cmp(N1, N2) <= 0;
  }

  bool operator<=(const BigInt& N1, const MachineInt& n2)
  {
    return cmp(N1, n2) <= 0;
  }
			
  bool operator<=(const MachineInt& n1, const BigInt& N2)
  {
    return cmp(n1, N2) <= 0;
  }
			

  double mantissa(const BigInt& N)
  {
    long DiscardExponent;
    return mpz_get_d_2exp(&DiscardExponent, mpzref(N));
  }


  long exponent(const BigInt& N)
  {
    return NumericCast<long>(mpz_sizeinbase(mpzref(N), 2));
  }


  long NumDigits(const BigInt& N, long base)
  {
    if (base < 2 || base > 36)
      CoCoA_ERROR(ERR::BadNumBase, "NumDigits(BigInt,base)");
    return NumericCast<long>(mpz_sizeinbase(mpzref(N), base));
  }


  double log(const BigInt& N)
  {
    if (IsZero(N))
      CoCoA_ERROR(ERR::LogZero, "log(BigInt)");
    if (exponent(N) < numeric_limits<double>::max_exponent)
      return std::log(std::abs(ConvertTo<double>(N)));
    long exp;
    const double mantissa = mpz_get_d_2exp(&exp, mpzref(N));
    const static double log2 = std::log(2.0);
    return std::log(std::abs(mantissa)) + exp*log2;
  }


  long ILogBase(const MachineInt& n, const MachineInt& base)
  {
    return ILogBase(BigInt(n), BigInt(base));
  }

  long ILogBase(const MachineInt& n, const BigInt& base)
  {
    return ILogBase(BigInt(n), base);
  }

  long ILogBase(const BigInt& N, const MachineInt& base)
  {
    return ILogBase(N, BigInt(base));
  }

  long ILogBase(const BigInt& N, const BigInt& base)
  {
    if (abs(base) < 2) CoCoA_ERROR(ERR::BadArg, "ILogBase: base must be at least 2");
    if (IsZero(N)) CoCoA_ERROR(ERR::LogZero, "ILogBase");
    const double ApproxLog = log(N)/log(base);
    // Check whether ApproxLog is far from integer.
    const double delta = 5* ApproxLog * numeric_limits<double>::epsilon();
    const long candidate = static_cast<long>(floor(ApproxLog+delta)); // probably right but could be too big by 1
    if (std::abs(ApproxLog - candidate) > delta)
      return candidate;
    // ApproxLog was close to integer, so we make a slow but sure check.
    const BigInt pwr = power(base, candidate);
    const int test = cmp(abs(N), pwr);
    if (test == 0) return candidate;
    if (test < 0) return candidate-1;
    if (abs(N) >= base*pwr) return candidate+1;
    return candidate;
  }


  // This fn comes out a long horrible mess because I cannot return a
  // machine int as a BigInt.  There must be a design error somewhere.
  BigInt RangeFactorial(const MachineInt& lo, const MachineInt& hi)
  {
    BigInt ans;
    // If zero is in the range, result is zero
    if ((IsZero(lo) || IsNegative(lo)) && !IsNegative(hi)) return ans;
    ans = 1;
    // First some checks for empty products.
    if (IsNegative(hi) && !IsNegative(lo)) return ans;
    // We now know that hi and lo have the same sign.
    // Two more checks for empty products.
    if (IsNegative(lo) && AsSignedLong(lo) > AsSignedLong(hi)) return ans;
    if (!IsNegative(lo) && AsUnsignedLong(lo) > AsUnsignedLong(hi)) return ans;
    // We now know that the product is not empty & does not span 0.
    unsigned long l = abs(lo);
    unsigned long h = abs(hi);
    if (IsNegative(hi))
    {
      std::swap(l, h);
      if (((h^l)&1) == 0) negate(ans);
    }

    if (h-l > 15)
    {
      const unsigned long mid = l + (h-l)/2; // equiv to (l+h)/2 but avoids overflow
      return RangeFactorial(l,mid) * RangeFactorial(1+mid,h);
    }
    // From here on 0 < l < h.
    for (; l <= h; ++l)
    {
      ans *= l;
//      mpz_mul_ui(mpzref(ans), mpzref(ans), l);
    }
    return ans;
  }

  const BigInt factorial(const BigInt& N)
  {
    if (N < 0)
      CoCoA_ERROR(ERR::NotNonNegative, "factorial(BigInt)");
    unsigned long n;
    if (!IsConvertible(n, N))
      CoCoA_ERROR(ERR::ArgTooBig, "factorial(BigInt)");
//    return factorial(n);
    BigInt ans;
    mpz_fac_ui(mpzref(ans), n);
    return ans;
  }


  const BigInt factorial(const MachineInt& n)
  {
    if (IsNegative(n))
      CoCoA_ERROR(ERR::NotNonNegative, "factorial(MachineInt)");
    BigInt ans;
    mpz_fac_ui(mpzref(ans), AsUnsignedLong(n));
    return ans;
  }


// From Wikipedia; apparently due to Ramanujan
// Experimentally max error is 4.6*10^(-8) at n=24.
double LogFactorial(const MachineInt& N)
{
  using std::log;
  using std::atan;
  if (IsNegative(N)) CoCoA_ERROR(ERR::NotNonNegative, "LogFactorial");
  const long n = AsSignedLong(N);
  const static double LogSqrtPi = log(4*atan(1.0))/2;
  if (n > 23)
    return n*log(double(n))-n + log(n*(1+n*(4.0+8*n)))/6 + LogSqrtPi;
  double fact = 1.0;
  for (long k=2; k <= n; ++k)
    fact *= k;
  return log(fact);
}

double LogFactorial(const BigInt& N)
{
  if (N < 0) CoCoA_ERROR(ERR::NotNonNegative, "LogFactorial");
  double n;
  if (!IsConvertible(n, N)) CoCoA_ERROR(ERR::ArgTooBig, "LogFactorial");
  {
    long SmallN;
    if (IsConvertible(SmallN, N)) return LogFactorial(SmallN);
  }
  using std::log;
  using std::atan;
  const static double LogSqrtPi = log(4*atan(1.0))/2;
  return n*log(n)-n + log(n*(1+n*(4.0+8*n)))/6 + LogSqrtPi;
}



  const BigInt binomial(const BigInt& N, const BigInt& R)
  {
//    if (R < 0) return BigInt(0);
    if (R < 0)
      CoCoA_ERROR(ERR::BadArg, "binomial(BigInt,BigInt): 2nd arg must be non-neg");
    if (N < 0)
    {
      if (IsEven(R)) return binomial(R-N-1,R);
      else return -binomial(R-N-1,R);
    }
    if (R > N)
      return BigInt(0);

    // General case:  N >= R >= 0
    unsigned long r;
    const bool RIsSmall = (2*R < N);
    if ((RIsSmall && !IsConvertible(r, R)) ||
        (!RIsSmall && !IsConvertible(r, N-R)))
      CoCoA_ERROR(ERR::ArgTooBig, "binomial(BigInt, BigInt)");

    BigInt ans;
    mpz_bin_ui(mpzref(ans), mpzref(N), r); // we know that r >= 0.
    return ans;
  }


  const BigInt binomial(const MachineInt& n, const BigInt& R)
  {
    if (R < 0)
      CoCoA_ERROR(ERR::BadArg, "binomial(MachineInt,BigInt): 2nd arg must be non-neg");
    return binomial(BigInt(n), R);
  }


  const BigInt binomial(const BigInt& N, const MachineInt& r)
  {
//    if (IsNegative(r)) return BigInt(0);
    if (IsNegative(r))
      CoCoA_ERROR(ERR::BadArg, "binomial(BigInt,MachineInt): 2nd arg must be non-neg");
    BigInt ans;
    mpz_bin_ui(mpzref(ans), mpzref(N), AsUnsignedLong(r));
    return ans;
  }


  const BigInt binomial(const MachineInt& n, const MachineInt& r)
  {
//    if (IsNegative(r)) return BigInt(0);
    if (IsNegative(r))
      CoCoA_ERROR(ERR::BadArg, "binomial(MachineInt,MachineInt): 2nd arg must be non-neg");

    BigInt ans;
    if (IsNegative(n))
      mpz_bin_ui(mpzref(ans), mpzref(BigInt(n)), AsUnsignedLong(r));
    else
      mpz_bin_uiui(mpzref(ans), AsUnsignedLong(n), AsUnsignedLong(r));
    return ans;
  }


  const BigInt fibonacci(const BigInt& N)
  {
    unsigned long n;
    if (!IsConvertible(n, N))
      CoCoA_ERROR(ERR::ArgTooBig, "fibonacci(BigInt)");
    return fibonacci(n);
  }


  const BigInt fibonacci(const MachineInt& n)
  {
    const bool EvenNegative = IsNegative(n) && IsEven(n);
    BigInt ans;
    unsigned long arg = abs(n);
    mpz_fib_ui(mpzref(ans), arg);
    if (EvenNegative)
      negate(ans);
    return ans;
  }


  // NB halves are rounded away from zero.
  const BigInt RoundDiv(const BigInt& num, const BigInt& den)
  {
    if (IsZero(den))
      CoCoA_ERROR(ERR::DivByZero, "RoundDiv(BigInt,BigInt)");
    BigInt q,r;
    quorem(q, r, num, den);
    if (CmpAbs(2*r,den) >= 0) // ROUND AWAY FROM ZERO
//    if (CmpAbs(2*r,den) > 0) // ROUND TOWARDS ZERO
      q += sign(r)*sign(den);
    return q;
  }

  const BigInt RoundDiv(const BigInt& num, const MachineInt& den)
  {
    if (IsZero(den))
      CoCoA_ERROR(ERR::DivByZero, "RoundDiv(BigInt,MachInt)");
    return RoundDiv(num, BigInt(den));
  }

  const BigInt RoundDiv(const MachineInt& num, const BigInt& den)
  {
    if (IsZero(den))
      CoCoA_ERROR(ERR::DivByZero, "RoundDiv(MachInt,BigInt)");
    return RoundDiv(BigInt(num), den);
  }

  // NB Halves are rounded away from zero.
  long RoundDiv(const MachineInt& n, const MachineInt& d)
  {
    if (IsZero(d))
      CoCoA_ERROR(ERR::DivByZero, "RoundDiv(n,d)");
    if (IsZero(n)) return 0;
    // Henceforth n and d are both non-zero.
    const long s = sign(n)*sign(d);
    const unsigned long q = uround_half_up(abs(n), abs(d)); // ROUND AWAY FROM ZERO
//    const unsigned long q = uround_half_down(abs(n), abs(d)); // ROUND TOWARDS ZERO
    if (s > 0) return NumericCast<long>(q);
    static const long MinLong = numeric_limits<long>::min();
    if (q + MinLong == 0) return MinLong; // to avoid overflow if result is MinLong.
    return -NumericCast<long>(q);

  }


  const BigInt iroot(const MachineInt& n, const MachineInt& r)
  {
    return iroot(BigInt(n), r);
  }

  const BigInt iroot(const MachineInt& n, const BigInt& R)
  {
    if (R < 1) CoCoA_ERROR(ERR::BadArg, "iroot: 2nd arg must be at least 1");
    unsigned long r;
    if (!IsConvertible(r, R))
      CoCoA_ERROR(ERR::ArgTooBig, "iroot: 2nd arg is too big");
    return iroot(BigInt(n), r);
  }

  const BigInt iroot(const BigInt& N, const MachineInt& r)
  {
    if (IsNegative(r) || IsZero(r))
      CoCoA_ERROR(ERR::BadArg, "iroot: 2nd arg must be at least 1");
    if (N < 0 && IsEven(r))
      CoCoA_ERROR(ERR::BadArg, "iroot: even root of negative number");
    BigInt ans;
    mpz_root(mpzref(ans), mpzref(N), AsUnsignedLong(r));
    return ans;
  }

  const BigInt iroot(const BigInt& N, const BigInt& R)
  {
    if (R < 1) CoCoA_ERROR(ERR::BadArg, "iroot: 2nd arg must be at least 1");
    unsigned long r;
    if (!IsConvertible(r, R))
      CoCoA_ERROR(ERR::ArgTooBig, "iroot: 2nd arg is too big");
    return iroot(N, r);
  }


  bool IsExactIroot(long& ans, const MachineInt& n, const MachineInt& r)
  {
    BigInt root;
    const bool IsExact = IsExactIroot(root, BigInt(n), r);
    if (!IsConvertible(ans, root))
      CoCoA_ERROR(ERR::ArgTooBig, "IsExactIroot");
    return IsExact;
  }

  bool IsExactIroot(BigInt& ans, const MachineInt& n, const MachineInt& r)
  {
    return IsExactIroot(ans, BigInt(n), r);
  }

  bool IsExactIroot(long& ans, const MachineInt& n, const BigInt& R)
  {
    unsigned long r;
    if (!IsConvertible(r, R))
      CoCoA_ERROR(ERR::ArgTooBig, "IsExactIroot: root index must be >=1 & fit into a ulong");
    BigInt root;
    const bool IsExact = IsExactIroot(root, BigInt(n), r);
    if (!IsConvertible(ans, root))
      CoCoA_ERROR(ERR::ArgTooBig, "IsExactIroot");
    return IsExact;
  }

  bool IsExactIroot(BigInt& ans, const MachineInt& n, const BigInt& R)
  {
    unsigned long r;
    if (!IsConvertible(r, R))
      CoCoA_ERROR(ERR::ArgTooBig, "IsExactIroot: root index must be >=1 & fit into a ulong");
    return IsExactIroot(ans, BigInt(n), r);
  }

  bool IsExactIroot(BigInt& ans, const BigInt& N, const MachineInt& r)
  {
    if (IsNegative(r) || IsZero(r))
      CoCoA_ERROR(ERR::BadArg, "IsExactIroot: root index must be at least 1");
    if (N < 0 && IsEven(AsUnsignedLong(r)))
      CoCoA_ERROR(ERR::BadArg, "IsExactIroot: even root of negative number");
    const bool IsExact = mpz_root(mpzref(ans), mpzref(N), AsUnsignedLong(r));
    return IsExact;
  }

  bool IsExactIroot(BigInt& ans, const BigInt& N, const BigInt& R)
  {
    unsigned long r;
    if (!IsConvertible(r, R))
      CoCoA_ERROR(ERR::ArgTooBig, "IsExactIroot: root index must be >=1 & fit into a ulong");
    return IsExactIroot(ans, N, r);
  }


  long isqrt(const MachineInt& n)
  {
    if (IsNegative(n)) CoCoA_ERROR(ERR::NotNonNegative, "isqrt(n)");
    long ans;
    IsConvertible(ans, isqrt(BigInt(n))); // conversion will succeed.
    return ans;
  }

  const BigInt isqrt(const BigInt& N)
  {
    if (N < 0) CoCoA_ERROR(ERR::NotNonNegative, "isqrt(N)");
    BigInt ans;
    mpz_sqrt(mpzref(ans), mpzref(N));
    return ans;
  }


  bool IsSquare(const MachineInt& n)
  {
    return IsSquare(BigInt(n));
  }

  bool IsSquare(const BigInt& N)
  {
    return mpz_perfect_square_p(mpzref(N));
  }


  bool IsPower(const MachineInt& n)
  {
    return IsPower(BigInt(n));
  }

  bool IsPower(const BigInt& N)
  {
    return mpz_perfect_power_p(mpzref(N));
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/IntOperations.C,v 1.11 2014/05/16 14:15:29 abbott Exp $
// $Log: IntOperations.C,v $
// Revision 1.11  2014/05/16 14:15:29  abbott
// Summary: Rewrote DoundDiv for MachineInts to handle rounding away from zero correctly when result is MinLong
// Author: JAA
//
// Revision 1.10  2014/05/16 12:26:11  abbott
// Summary: Changed rounding (in RoundDiv): it is now "away from zero"
// Author: JAA
//
// Revision 1.9  2014/05/05 18:21:00  abbott
// Removed another M_PI (replaced by 4*atan(1.0))
//
// Revision 1.8  2014/04/30 16:08:42  abbott
// Summary: Removed dependency on nonstd M_PI; now uses atan; cast arg of log to double
// Author: JAA
//
// Revision 1.7  2014/04/10 15:30:33  abbott
// Summary: Corrected impl of IsPowerOf2
// Author: JAA
//
// Revision 1.6  2014/04/08 13:04:13  abbott
// Summary: Added new fn IsPowerOf2
// Author: JAA
//
// Revision 1.5  2014/03/06 15:49:30  abbott
// Summary: Zero to power zero now gives 1
// Author: JAA
//
// Revision 1.4  2013/05/20 15:48:13  abbott
// Added new fn LogFactorial (placed in IntOperations).
//
// Revision 1.3  2013/05/14 14:23:00  abbott
// Minor improvement to log(BigInt)
//
// Revision 1.2  2012/12/04 09:55:47  abbott
// Added new fns LeastNNegRemainder and SymmRemainder (with doc & new test).
// Some minor corrections to the doc (for operator/ and operator%).
//
// Revision 1.1  2012/05/28 09:18:21  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
//
