//   Copyright (c)  2009-2010,2013  John Abbott

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


#include "CoCoA/BigRat.H"

#include "CoCoA/IntOperations.H"
#include "CoCoA/NumTheory.H"
#include "CoCoA/OpenMath.H"
#include "CoCoA/error.H"
#include "CoCoA/utils_gmp.H"

#include <cmath>
using std::abs;
using std::floor;
#include <iostream>
using std::ostream;
using std::istream;
#include <limits>
using std::numeric_limits;
#include <string>
using std::string;
#include <vector>
using std::vector;

namespace CoCoA
{
  
//????????  BigRat::BigRat(MP_RAT to_be_owned); // for efficient return values in operator+ etc.

  BigRat::BigRat()
  {
    mpq_init(myRep);
  }

  BigRat::BigRat(const mpq_t q)
  {
    mpq_init(myRep);
    mpq_set(myRep, q);
  }

//   BigRat::BigRat(const MachineInt& n)
//   {
//     mpq_init(myRep);
//     if (IsNegative(n))
//       mpq_set_si(myRep, AsSignedLong(n), 1);
//     else
//       mpq_set_ui(myRep, AsUnsignedLong(n), 1);
//   }

//   BigRat::BigRat(const BigInt& N)
//   {
//     mpq_init(myRep);
//     mpq_set_z(myRep, mpzref(N));
//   }

  BigRat::BigRat(const MachineInt& n1, const MachineInt& n2, ReduceFlag status)
  {
//    if (IsZero(n2))
//      CoCoA_ERROR(ERR::DivByZero, "BigRat(n1,n2)");
    mpq_init(myRep);
    myAssign(BigInt(n1), BigInt(n2), status);
//    BELOW IS ORIGINAL CODE -- slightly better 'cos does not create temporaries.
//     const bool IsNegativeFraction = IsNegative(n1) ^ IsNegative(n2);
//     mpq_set_ui(myRep, abs(n1), abs(n2));
//     if (status == NotReduced)
//       mpq_canonicalize(myRep);
//     else
//       CoCoA_ASSERT(gcd(n1,n2) == 1);
//     if (IsNegativeFraction)
//       mpq_neg(myRep, myRep);
  }

  BigRat::BigRat(const MachineInt& n1, const BigInt& N2, ReduceFlag status)
  {
    if (IsZero(N2))
      CoCoA_ERROR(ERR::DivByZero, "BigRat(n1,N2)");
    mpq_init(myRep);
    myAssign(BigInt(n1), N2, status);
  }

  BigRat::BigRat(const BigInt& N1, const MachineInt& n2, ReduceFlag status)
  {
    if (IsZero(n2))
      CoCoA_ERROR(ERR::DivByZero, "BigRat(N1,n2)");
    mpq_init(myRep);
    myAssign(N1, BigInt(n2), status);
  }

  BigRat::BigRat(const BigInt& N1, const BigInt& N2, ReduceFlag status)
  {
    if (IsZero(N2))
      CoCoA_ERROR(ERR::DivByZero, "BigRat(N1,N2)");
    mpq_init(myRep);
    myAssign(N1, N2, status);
  }


  BigRat::BigRat(const std::string& str, ReduceFlag status)
  {
    mpq_init(myRep);
//     if (base != 0 && (base < 2 || base > 36))
//       CoCoA_ERROR(ERR::BadNumBase, "BigRat(string,int)");
    if (mpq_set_str(myRep, str.c_str(), 10) != 0)
      CoCoA_ERROR(ERR::BadArg, "BigRat(string)");
    if (status == NotReduced)
      mpq_canonicalize(myRep);
  }


  BigRat::BigRat(const BigRat& from)
  {
    mpq_init(myRep);
    mpq_set(myRep, from.myRep);
  }


// BigRat& clone() const;

  BigRat::~BigRat()
  {
    mpq_clear(myRep);
  }


  // NOTE: This is NOT EXCEPTION CLEAN if the GMP fns can throw.
  void BigRat::myAssign(const BigInt& N1, const BigInt& N2, ReduceFlag status/*=NotReduced*/)
  {
    CoCoA_ASSERT(!IsZero(N2) || status == AlreadyReduced);
    const bool IsNegativeFraction = (N1 < 0) ^ (N2 < 0);
    mpz_abs(mpq_numref(myRep), mpzref(N1));
    mpz_abs(mpq_denref(myRep), mpzref(N2));
    if (status == NotReduced)
      mpq_canonicalize(myRep);
    else
      CoCoA_ASSERT(gcd(N1,N2) == 1);
    if (IsNegativeFraction)
      mpq_neg(myRep, myRep);
  }


  BigRat& BigRat::operator=(const BigRat& rhs)
  {
    mpq_set(myRep, rhs.myRep);
    return *this;
  }

    // -------- functions that modify at least one argument or `*this' ----------

  BigRat& BigRat::operator+=(const BigRat& rhs)
  {
    mpq_add(myRep, myRep, rhs.myRep);
    return *this;
  }

  BigRat& BigRat::operator-=(const BigRat& rhs)
  {
    mpq_sub(myRep, myRep, rhs.myRep);
    return *this;
  }

  BigRat& BigRat::operator*=(const BigRat& rhs)
  {
//     if (mpz_sgn(mpq_numref(myRep)) == 0) return *this;
//     if (mpz_sgn(mpq_numref(rhs.myRep)) == 0) return operator=(rhs);
    mpq_mul(myRep, myRep, rhs.myRep);
    return *this;
  }

  BigRat& BigRat::operator/=(const BigRat& rhs)
  {
    if (mpz_sgn(mpq_numref(rhs.myRep)) == 0)
      CoCoA_ERROR(ERR::DivByZero, "q1 /= q2");
    mpq_div(myRep, myRep, rhs.myRep);
    return *this;
  }
                        
    // Same but with RHS a BigInt...
  BigRat& BigRat::operator=(const BigInt& rhs)
  {
    mpq_set_z(myRep, mpzref(rhs));
    return *this;
  }


  // SLUG: impl is needlessly inefficient: makes useless copy in D
  BigRat& BigRat::operator+=(const BigInt& rhs)
  {
    const BigInt D(mpq_denref(myRep));
    const BigInt tmp = rhs*D;
    mpz_add(mpq_numref(myRep), mpq_numref(myRep), mpzref(tmp));
    // no need to call mpq_canonicalize
    return *this;
  }

  // SLUG: impl is needlessly inefficient: makes useless copy in D
  BigRat& BigRat::operator-=(const BigInt& rhs)
  {
    const BigInt D(mpq_denref(myRep));
    const BigInt tmp = rhs*D;
    mpz_sub(mpq_numref(myRep), mpq_numref(myRep), mpzref(tmp));
    // no need to call mpq_canonicalize
    return *this;
  }

  BigRat& BigRat::operator*=(const BigInt& rhs)
  {
    return operator*=(BigRat(rhs,1));
  }

  BigRat& BigRat::operator/=(const BigInt& rhs)
  {
    if (IsZero(rhs))
      CoCoA_ERROR(ERR::DivByZero, "Q /= N");
    // Could be more efficient if "*this" is 0.
    return operator/=(BigRat(rhs,1));
  }
                        
    // Same but with RHS a MachinInteger...
  BigRat& BigRat::operator= (const MachineInt& rhs)
  {
    if (IsNegative(rhs))
      mpq_set_si(myRep, AsSignedLong(rhs), 1);
    else
      mpq_set_ui(myRep, AsUnsignedLong(rhs), 1);
    return *this;
  }

  BigRat& BigRat::operator+=(const MachineInt& rhs)
  {
    return operator+=(BigInt(rhs));
  }

  BigRat& BigRat::operator-=(const MachineInt& rhs)
  {
    return operator-=(BigInt(rhs));
  }

  BigRat& BigRat::operator*=(const MachineInt& rhs)
  {
    return operator*=(BigInt(rhs));
  }

  BigRat& BigRat::operator/=(const MachineInt& rhs)
  {
    if (IsZero(rhs))
      CoCoA_ERROR(ERR::DivByZero, "Q /= n");
    return operator/=(BigInt(rhs));
  }


  const BigRat& BigRat::operator++()
  {
    mpz_add(mpq_numref(myRep), mpq_numref(myRep), mpq_denref(myRep)); // no need to reduce
    return *this;
  }

  const BigRat& BigRat::operator--()
  {
    mpz_sub(mpq_numref(myRep), mpq_numref(myRep), mpq_denref(myRep));
    return *this;
  }

  const BigRat BigRat::operator++(int) // INEFFICIENT
  {
    BigRat ans(*this);
    operator++();
    return ans;
  }

  const BigRat BigRat::operator--(int) // INEFFICIENT
  {
    BigRat ans(*this);
    operator--();
    return ans;
  }



  // I/O FUNCTIONS

  string ConvertToString(const BigRat& src, int base/*=10*/)
  {
    if (base < 2 || base > 36)
      CoCoA_ERROR(ERR::BadNumBase, "IsConvertible(string,BigRat,int)");
    const long digits = NumDigits(num(src),base) + NumDigits(den(src),base);
    vector<char> buffer(digits+3); // +2 to allow for minus sign, "/" character and terminating NUL
    mpq_get_str(&buffer[0], base, mpqref(src));
    return &buffer[0];
  }


  std::ostream& operator<<(std::ostream& out, const BigRat& Q)
  {
    out << num(Q);
    if (!IsOneDen(Q))
      out << "/" << den(Q);
    return out;
  }

  std::istream& operator>>(std::istream& in, BigRat& Q)
  {
    if (!in.good()) return in;
    BigInt N;
    in >> N;
    if (!in.good()) return in;
    const char slash = in.peek();
    if (slash != '/')
    {
      Q = N;
      return in;
    }
    in.ignore();
    BigInt D;
    in >> D;
    if (!in.good()) return in;
    BigRat ans(N,D);
    swap(Q, ans); // really an assignment
    return in;
  }


  OpenMathOutput& operator<<(OpenMathOutput& OMOut, const BigRat& Q)
  {
    OMOut->mySendApplyStart();
    OMOut->mySendSymbol("nums1","rational");
    OMOut->mySend(num(Q));
    OMOut->mySend(den(Q));
    OMOut->mySendApplyEnd();
    return OMOut;
  }


  OpenMathInput& operator>>(OpenMathInput& OMIn, BigRat& /*Q*/)
  {
    CoCoA_ERROR(ERR::NYI, "OpenMathInput fn for BigRat");
    return OMIn;
  }


  // STANDARD ARITHMETIC OPERATIONS

  void swap(BigRat& a, BigRat& b)
  {
    mpq_swap(mpqref(a), mpqref(b));
  }

  const BigInt num(const BigRat& Q)
  {
    return BigInt(mpq_numref(mpqref(Q)));
  }

  const BigInt den(const BigRat& Q)
  {
    return BigInt(mpq_denref(mpqref(Q)));
  }

  const BigRat abs(const BigRat& Q)
  {
    BigRat ans;
    mpq_abs(mpqref(ans), mpqref(Q));
    return ans;
  }

  const BigRat operator-(const BigRat& Q)
  {
    BigRat ans;
    mpq_neg(mpqref(ans), mpqref(Q));
    return ans;
  }

  const BigRat operator+(const BigRat& Q1, const BigRat& Q2)
  {
    BigRat ans;
    mpq_add(mpqref(ans), mpqref(Q1), mpqref(Q2));
    return ans;
  }

  const BigRat operator-(const BigRat& Q1, const BigRat& Q2)
  {
    BigRat ans;
    mpq_sub(mpqref(ans), mpqref(Q1), mpqref(Q2));
    return ans;
  }

  const BigRat operator*(const BigRat& Q1, const BigRat& Q2)
  {
    BigRat ans;
    mpq_mul(mpqref(ans), mpqref(Q1), mpqref(Q2));
    return ans;
  }

  const BigRat operator/(const BigRat& Q1, const BigRat& Q2)
  {
    if (IsZero(Q2))
      CoCoA_ERROR(ERR::DivByZero, "Q1/Q2");
    BigRat ans;
    mpq_div(mpqref(ans), mpqref(Q1), mpqref(Q2));
    return ans;
  }


  const BigRat operator+(const BigRat& Q, const BigInt& N)
  {
    BigRat ans = Q;
    return ans += N;
//   THE LINES BELOW SHOULD BE MORE EFFICIENT (but they're not very readable).
//     BigRat ans;
//     mpz_mul(mpq_denref(ans), mpzref(N), mpq_denref(mpqref(Q)));
//     mpz_add(mpq_numref(mpqref(ans)), mpq_denref(mpqref(ans)), mpq_numref(mpqref(Q)));
//     mpz_set(mpq_denref(mpqref(ans)), mpq_denref(mpqref(Q)));
//     return ans;
  }

  const BigRat operator-(const BigRat& Q, const BigInt& N)
  {
    BigRat ans(Q);
    return ans -= N;
  }

  const BigRat operator*(const BigRat& Q, const BigInt& N)
  {
    BigRat ans(Q);
    return ans *= N;
  }

  const BigRat operator/(const BigRat& Q, const BigInt& N)
  {
    if (IsZero(N))
      CoCoA_ERROR(ERR::DivByZero, "Q/N");
    BigRat ans(Q);
    return ans /= N;
  }

  const BigRat operator+(const BigInt& N, const BigRat& Q)
  {
    return Q+N;
  }

  const BigRat operator-(const BigInt& N, const BigRat& Q)
  {
    BigRat ans = Q-N;
    return -ans;
  }

  const BigRat operator*(const BigInt& N, const BigRat& Q)
  {
    return Q*N;
  }

  const BigRat operator/(const BigInt& N, const BigRat& Q)
  {
    if (IsZero(Q))
      CoCoA_ERROR(ERR::DivByZero, "N/Q");
    return BigRat(N,1)/Q;
  }


  const BigRat operator+(const BigRat& Q, const MachineInt& n)
  {
    return Q + BigRat(n,1);
  }

  const BigRat operator-(const BigRat& Q, const MachineInt& n)
  {
    return Q - BigRat(n,1);
  }

  const BigRat operator*(const BigRat& Q, const MachineInt& n)
  {
    return Q * BigRat(n,1);
  }

  const BigRat operator/(const BigRat& Q, const MachineInt& n)
  {
    if (IsZero(n))
      CoCoA_ERROR(ERR::DivByZero, "Q/n");
    return Q / BigRat(n,1);
  }


  const BigRat operator+(const MachineInt& n, const BigRat& Q)
  {
    return BigRat(n,1) + Q;
  }

  const BigRat operator-(const MachineInt& n, const BigRat& Q)
  {
    return BigRat(n,1) - Q;
  }

  const BigRat operator*(const MachineInt& n, const BigRat& Q)
  {
    return BigRat(n,1) * Q;
  }

  const BigRat operator/(const MachineInt& n, const BigRat& Q)
  {
    if (IsZero(Q))
      CoCoA_ERROR(ERR::DivByZero, "n/Q");
    return BigRat(n,1) / Q;
  }

  const BigRat power(const BigRat& base, const BigInt& exponent)
  {
    if (exponent >= 0)
      return BigRat(power(num(base), exponent), power(den(base), exponent), BigRat::AlreadyReduced);
    if (IsZero(base))
      CoCoA_ERROR(ERR::BadPwrZero, "power(BigRat,BigInt)");
    return BigRat(power(den(base), -exponent), power(num(base), -exponent), BigRat::AlreadyReduced);
  }

  const BigRat power(const BigRat& base, const MachineInt& exponent)
  {
    if (!IsNegative(exponent))
      return BigRat(power(num(base), exponent), power(den(base), exponent), BigRat::AlreadyReduced);
    if (IsZero(base))
      CoCoA_ERROR(ERR::BadPwrZero, "power(BigRat,MachineInt)");
    return BigRat(power(den(base), abs(exponent)), power(num(base), abs(exponent)), BigRat::AlreadyReduced);
  }



  bool IsZero(const BigRat& Q)
  {
    return IsZero(num(Q));
  }


  bool IsOne(const BigRat& Q)
  {
    return mpq_cmp_ui(mpqref(Q), 1,1) == 0;
  }


  bool IsMinusOne(const BigRat& Q)
  {
    return mpq_cmp_si(mpqref(Q), -1,1) == 0;
  }


  bool IsOneNum(const BigRat& Q)
  {
    return mpz_cmp_ui(mpq_numref(mpqref(Q)), 1) == 0;
  }


  bool IsOneDen(const BigRat& Q)
  {
    return mpz_cmp_ui(mpq_denref(mpqref(Q)), 1) == 0;
  }


  int sign(const BigRat& Q)
  {
    return mpq_sgn(mpqref(Q));
  }


  // COMPARISON FUNCTIONS

  int cmp(const BigRat& Q1, const BigRat& Q2)
  {
    return sign(mpq_cmp(mpqref(Q1), mpqref(Q2)));
  }


  int cmp(const BigRat& Q, const BigInt& N)
  {
    return cmp(num(Q), N*den(Q));
  }


  int cmp(const BigInt& N, const BigRat& Q)
  {
    return cmp(N*den(Q), num(Q));
  }


  int cmp(const BigRat& Q, const MachineInt& n)
  {
    if (IsNegative(n))
      return sign(mpq_cmp_si(mpqref(Q), AsSignedLong(n),1));
    return sign(mpq_cmp_ui(mpqref(Q), AsUnsignedLong(n),1));
  }


  int cmp(const MachineInt& n, const BigRat& Q)
  {
    return -cmp(Q, n);
  }


  int CmpAbs(const BigRat& Q1, const BigRat& Q2)
  {
    return mpq_cmpabs(mpqref(Q1), mpqref(Q2));
  }

  // The next 4 prefer simplicity over speed.
  int CmpAbs(const BigRat& Q, const BigInt& N)
  {
    return CmpAbs(Q, BigRat(N,1));
  }

  int CmpAbs(const BigInt& N, const BigRat& Q)
  {
    return CmpAbs(BigRat(N,1), Q);
  }

  int CmpAbs(const BigRat& Q, const MachineInt& n)
  {
    return CmpAbs(Q, BigRat(n,1));
  }

  int CmpAbs(const MachineInt& n, const BigRat& Q)
  {
    return CmpAbs(BigRat(n,1), Q);
  }



  bool operator==(const BigRat& Q1, const BigRat& Q2)
  {
    return mpq_equal(mpqref(Q1), mpqref(Q2));
  }

  bool operator!=(const BigRat& Q1, const BigRat& Q2)
  {
    return !(Q1 == Q2);
  }

  bool operator> (const BigRat& Q1, const BigRat& Q2)
  {
    return cmp(Q1,Q2) > 0;
  }

  bool operator>=(const BigRat& Q1, const BigRat& Q2)
  {
    return cmp(Q1,Q2) >= 0;
  }

  bool operator< (const BigRat& Q1, const BigRat& Q2)
  {
    return cmp(Q1,Q2) < 0;
  }

  bool operator<=(const BigRat& Q1, const BigRat& Q2)
  {
    return cmp(Q1,Q2) <= 0;
  }

                        
  bool operator==(const BigRat& Q, const BigInt& N)
  {
    return IsOneDen(Q) && num(Q) == N;
  }

  bool operator!=(const BigRat& Q, const BigInt& N)
  {
    return !(Q == N);
  }

  bool operator> (const BigRat& Q, const BigInt& N)
  {
    return cmp(Q,N) > 0;
  }

  bool operator>=(const BigRat& Q, const BigInt& N)
  {
    return cmp(Q,N) >= 0;
  }

  bool operator< (const BigRat& Q, const BigInt& N)
  {
    return cmp(Q,N) < 0;
  }

  bool operator<=(const BigRat& Q, const BigInt& N)
  {
    return cmp(Q,N) <= 0;
  }

                        
  bool operator==(const BigInt& N, const BigRat& Q)
  {
    return Q == N;
  }

  bool operator!=(const BigInt& N, const BigRat& Q)
  {
    return !(Q == N);
  }

  bool operator> (const BigInt& N, const BigRat& Q)
  {
    return cmp(N,Q) > 0;
  }

  bool operator>=(const BigInt& N, const BigRat& Q)
  {
    return cmp(N,Q) >= 0;
  }

  bool operator< (const BigInt& N, const BigRat& Q)
  {
    return cmp(N,Q) < 0;
  }

  bool operator<=(const BigInt& N, const BigRat& Q)
  {
    return cmp(N,Q) <= 0;
  }

                        
  bool operator==(const BigRat& Q, const MachineInt& n)
  {
    return IsOneDen(Q) && num(Q) == n;
  }

  bool operator!=(const BigRat& Q, const MachineInt& n)
  {
    return !(Q == n);
  }

  bool operator> (const BigRat& Q, const MachineInt& n)
  {
    return cmp(Q,n) > 0;
  }

  bool operator>=(const BigRat& Q, const MachineInt& n)
  {
    return cmp(Q,n) >= 0;
  }

  bool operator< (const BigRat& Q, const MachineInt& n)
  {
    return cmp(Q,n) < 0;
  }

  bool operator<=(const BigRat& Q, const MachineInt& n)
  {
    return cmp(Q,n) <= 0;
  }

                
  bool operator==(const MachineInt& n, const BigRat& Q)
  {
    return IsOneDen(Q) && num(Q) == n;
  }

  bool operator!=(const MachineInt& n, const BigRat& Q)
  {
    return !(n == Q);
  }

  bool operator> (const MachineInt& n, const BigRat& Q)
  {
    return cmp(n,Q) > 0;
  }

  bool operator>=(const MachineInt& n, const BigRat& Q)
  {
    return cmp(n,Q) >= 0;
  }

  bool operator< (const MachineInt& n, const BigRat& Q)
  {
    return cmp(n,Q) < 0;
  }

  bool operator<=(const MachineInt& n, const BigRat& Q)
  {
    return cmp(n,Q) <= 0;
  }

                        

  // MISCELLANEOUS FUNCTIONS

  double mantissa(const BigRat& /*Q*/)
  {
    CoCoA_ERROR(ERR::NYI, "mantissa(BigRat)");
    return 0.0;
  }

  long exponent(const BigRat& /*Q*/)
  {
    CoCoA_ERROR(ERR::NYI, "exponent(BigRat)");
    return 0;
  }

  // BUG BUG BUG Is this log(Q)  or log(abs(Q))???
  double log(const BigRat& Q)
  {
    if (IsZero(Q))
      CoCoA_ERROR(ERR::BadArg, "log(Q)");
    return log(num(Q)) - log(den(Q));
  }

  long ILogBase(const BigRat& Q, const MachineInt& base)
  {
    return ILogBase(Q, BigInt(base));
  }

  // BUGLY: Can we merge the two almost identical sections?
  long ILogBase(const BigRat& Q, BigInt base)
  {
    using std::abs;
    base = abs(base);
    if (base < 2) CoCoA_ERROR(ERR::BadArg, "ILogBase: base must be at least 2");
    if (IsZero(Q)) CoCoA_ERROR(ERR::BadArg, "ILogBase: cannot compute log(0)");
    const double ApproxLog = log(Q)/log(base);
    const double delta = 5 * abs(ApproxLog) * numeric_limits<double>::epsilon();
///???BUG not yet fully implemented    if (ApproxLog > numeric_limits<long>::max()) CoCoA_ERROR(ERR::ArgTooBig, "ILogBase");
    const long candidate = static_cast<long>(std::floor(ApproxLog+delta)); // probably right but could be too big by 1
    if (std::abs(ApproxLog - candidate) > delta)
      return candidate;
    if (candidate >= 0)
    {
      const BigInt pwr = power(base, candidate);
      const int test = cmp(abs(Q), pwr);
      if (test == 0) return candidate;
      if (test < 0) return candidate-1;
      if (abs(Q) >= base*pwr) return candidate+1;
      return candidate;
    }
    // candidate < 0
    const BigInt pwr = power(base, -candidate);
    const BigRat shifted = abs(Q)*pwr;
    const int test = cmp(shifted, 1);
    if (test == 0) return candidate;
    if (test < 0) return candidate-1;
    if (shifted >= base) return candidate+1;
    return candidate;
  }


  BigInt floor(const BigRat& Q)
  {
    BigInt ans;
    mpz_fdiv_q(mpzref(ans), mpq_numref(mpqref(Q)), mpq_denref(mpqref(Q)));
    return ans;
  }

  BigInt ceil(const BigRat& Q)
  {
    BigInt ans;
    mpz_cdiv_q(mpzref(ans), mpq_numref(mpqref(Q)), mpq_denref(mpqref(Q)));
    return ans;
  }


  // This impl clearly guarantees that rounding is compatible with RoundDiv!
  BigInt round(const BigRat& Q)
  {
    if (IsOneDen(Q)) return num(Q);
    return RoundDiv(num(Q), den(Q));
  }


} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/BigRat.C,v 1.12 2014/07/09 11:41:35 abbott Exp $
// $Log: BigRat.C,v $
// Revision 1.12  2014/07/09 11:41:35  abbott
// Summary: Corrected (embarassing) bug in ILogBase
// Author: JAA
//
// Revision 1.11  2014/07/07 12:11:12  abbott
// Summary: Corrected operator>> (forgot to ignore the '/')
// Author: JAA
//
// Revision 1.10  2014/06/14 19:25:11  abbott
// Summary: Added new fn CmpAbs (for BigRat)
// Author: JAA
//
// Revision 1.9  2014/05/16 12:02:28  abbott
// Summary: Changed comment about fn "round"
// Author: JAA
//
// Revision 1.8  2014/01/28 09:58:30  abbott
// Revised impl of ctor from std::string so that it accepts and respects 2nd arg saying whether the fraction should be canonicalized.
//
// Revision 1.7  2013/05/20 15:50:20  abbott
// Added new ctor for BigRat from mpq_t.
//
// Revision 1.6  2013/03/26 14:56:06  abbott
// Updated the conversion fns (in ptic removed procedure "convert");
// numerous consequential changes.
//
// Revision 1.5  2012/12/12 10:38:35  abbott
// Changed assertion to allow creation of 1/0 if marked as AlreadyReduced.
//
// Revision 1.4  2012/12/04 20:14:49  abbott
// Modified BigRat ctor to allow one to create 1/0 (if specified as AleadyReduced).
//
// Revision 1.3  2012/05/28 09:18:21  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.2  2011/11/09 14:03:40  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.1  2011/09/23 13:20:35  bigatti
// -- QQ.C renamed into BigRat.C
//
// Revision 1.18  2011/09/06 15:21:53  abbott
// Changed "cmp" functions so that the return value is in {-1,0,+1}.
//
// Revision 1.17  2011/08/24 10:28:49  bigatti
// -- renamed QQ --> BigRat
// -- sorted #include
//
// Revision 1.16  2011/08/23 16:18:38  abbott
// Simplified defn of round; added comment about rounding halves towards zero.
//
// Revision 1.15  2011/08/17 11:57:39  abbott
// Added static_cast to keep compiler quiet.
//
// Revision 1.14  2011/08/14 15:52:17  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.13  2011/06/23 16:01:07  abbott
// Removed single arg ctor QQ(MachineInteger), & consequential changes.
//
// Revision 1.12  2011/03/01 15:26:10  abbott
// Improved impl of ILogBase -- faster in most cases.
//
// Revision 1.11  2011/02/25 12:06:51  abbott
// Added new fn IsOneNum; also some minor code cleaning in QQ.C
//
// Revision 1.10  2011/01/14 17:23:19  abbott
// Fixed a minor bug in power.
//
// Revision 1.9  2010/12/26 13:03:16  abbott
// Added ILogBase function (to ZZ & QQ).
//
// Revision 1.8  2010/05/07 14:57:52  abbott
// Two main changes:
//   power(QQ,ZZ) now allows negative exponent
//   renamed QQ::AlreadyNormalized to QQ::AlreadyReduced
//           (and allowed denoms to be negative; the ctor then makes them positive).
//
// Revision 1.7  2010/03/22 11:49:28  abbott
// Added ctor from a string.
//
// Revision 1.6  2010/03/18 16:40:42  abbott
// Added missing include directive.
//
// Revision 1.5  2010/03/18 16:34:10  abbott
// Added new pseudo-ctors for QQ with optional flag to indicate that value is already normalized.
// Added OpenMath I/O operators.
//
// Revision 1.4  2009/10/26 15:39:24  bigatti
// -- added CopyFromMPZ in ZZ ctor
//
// Revision 1.3  2009/07/08 12:26:53  abbott
// Added floor and ceil functions for QQs.
// Added example program for QQs.
// Minor correction to convert.C; minor cleaning to ex-ZZ1.C
//
// Revision 1.2  2009/07/06 12:31:26  abbott
// Commented out two unused function arguments (to keep compiler quiet
// when compiling a debugging version).
//
// Revision 1.1  2009/07/02 16:29:42  abbott
// Added new class QQ to represent rational numbers.
// Consequent change to the Makefile.
//
//
