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

#include "CoCoA/FloatApprox.H"

#include "CoCoA/BigRat.H"
#include "CoCoA/IntOperations.H"
#include "CoCoA/error.H"


#include <iostream>
using std::ostream;

namespace CoCoA
{

  const int MantExp2::ourDefaultMantBits = 53; // same as IEEE "double"

  MantExp2 MantissaAndExponent2(const MachineInt& n, const MachineInt& MantWidth)
  {
    return MantissaAndExponent2(BigRat(n,1), MantWidth);
  }


  // let rational version do the work so that halves are rounded consistently!
  MantExp2 MantissaAndExponent2(const BigInt& N, const MachineInt& MantWidth)
  {
    return MantissaAndExponent2(BigRat(N,1), MantWidth);
  }


  // Simple/compact rather than fast;  is speed so important here?
  MantExp2 MantissaAndExponent2(const BigRat& q, const MachineInt& MantWidth)
  {
    if (IsNegative(MantWidth) || !IsSignedLong(MantWidth) || abs(MantWidth) < 2)
      CoCoA_ERROR(ERR::BadArg, "MantissaAndExponent2");
    if (IsZero(q)) return MantExp2(0,0,BigInt(0),0);
    const long MantBits = AsSignedLong(MantWidth);
    BigInt N = abs(num(q));
    BigInt D = den(q);
    const int SignQ = sign(q);
    const long LogQ = ILogBase(q,2);
    const long exp = LogQ-MantBits+1;  // NB exp-1, exp+1, and -exp  will not overflow!
    if (exp <= 0)
      mpz_mul_2exp(mpzref(N), mpzref(N), -exp); // N *= 2^|exp|
    else
      mpz_mul_2exp(mpzref(D), mpzref(D), exp);  // D *= 2^exp

    N = RoundDiv(N,D);
    if (ILogBase(N,2) == MantBits) // true iff mantissa has "overflowed"
      return MantExp2(SignQ, 1+LogQ, N/2, MantBits);
    return MantExp2(SignQ, LogQ, N, MantBits);
  }


  std::ostream& operator<<(std::ostream& out, const MantExp2& ME)
  {
    out << "MantExp2(sign=" << ME.mySign << ", exp=" << ME.myExponent << ", mant=" << ME.myMantissa << ", NumDigits=" << ME.myNumDigits << ")";
    return out;
  }


  //------------------------------------------------------------------
  BigRat FloatApprox(const MachineInt& n, const MachineInt& MantBits)
  {
    const MantExp2 ME = MantissaAndExponent2(n, MantBits);
    return (ME.mySign * ME.myMantissa) * power(BigRat(2,1), ME.myExponent-ME.myNumDigits+1);
  }

  BigRat FloatApprox(const BigInt& n, const MachineInt& MantBits)
  {
    const MantExp2 ME = MantissaAndExponent2(n, MantBits);
    return (ME.mySign * ME.myMantissa) * power(BigRat(2,1), ME.myExponent-ME.myNumDigits+1);
  }

  BigRat FloatApprox(const BigRat& n, const MachineInt& MantBits)
  {
    const MantExp2 ME = MantissaAndExponent2(n, MantBits);
    return (ME.mySign * ME.myMantissa) * power(BigRat(2,1), ME.myExponent-ME.myNumDigits+1);
  }


  //------------------------------------------------------------------
  // Decimal "floating point" representation

  const int MantExp10::ourDefaultSigFig = 5;


  std::ostream& operator<<(std::ostream& out, const MantExp10& ME)
  {
    out << "MantExp10(sign=" << ME.mySign << ", exp=" << ME.myExponent << ", mant=" << ME.myMantissa << ", NumDigits=" << ME.myNumDigits << ")";
    return out;
  }


  MantExp10 MantissaAndExponent10(const BigInt& N, const MachineInt& SigFig)
  {
    if (IsNegative(SigFig) || !IsSignedLong(SigFig) || IsZero(SigFig))
      CoCoA_ERROR(ERR::BadArg, "MantissaAndExponent10");
    if (IsZero(N)) return MantExp10(0,0,BigInt(0),0);
    const long ndigits = AsSignedLong(SigFig);
    const int SignN = sign(N);
    const long e = ILogBase(N,10); // overflow???
    if (e < ndigits)
      return MantExp10(SignN, e, abs(N)*power(10, ndigits-e-1), ndigits);
    const BigInt HalfULP = 5*power(10, e-ndigits);
    const BigInt digits = (1+abs(N)/HalfULP)/2;
    // Must check whether digits has overflowed...
    if (abs(digits) == power(10, ndigits))
      return MantExp10(SignN, e+1, digits/10, ndigits);
    return MantExp10(SignN, e, digits, ndigits);
  }


  MantExp10 MantissaAndExponent10(const BigRat& q, const MachineInt& SigFig)
  {
    if (IsNegative(SigFig) || !IsSignedLong(SigFig) || IsZero(SigFig))
      CoCoA_ERROR(ERR::BadArg, "MantissaAndExponent10");
    if (IsZero(q)) return MantExp10(0,0,BigInt(0),0);
    if (IsOneDen(q)) return MantissaAndExponent10(num(q), SigFig);
    const long ndigits = AsSignedLong(SigFig);
    const int signq = sign(q);
    const long e = ILogBase(q,10); // overflow???
    BigInt digits;
    if (e < ndigits)
      digits = round(abs(q)*power(10,ndigits-e-1));
    else
      digits = round(abs(q)/power(10,1+e-ndigits)); 
    // Must check whether digits has overflowed...
    if (abs(digits) == power(10, ndigits))
      return MantExp10(signq, e+1, digits/10, ndigits);
    return MantExp10(signq, e, digits, ndigits);
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/FloatApprox.C,v 1.5 2014/05/14 13:18:11 abbott Exp $
// $Log: FloatApprox.C,v $
// Revision 1.5  2014/05/14 13:18:11  abbott
// Summary: Updated impls of FloatApprox to follow new defn of MantissaAndExponent2
// Author: JAA
//
// Revision 1.4  2014/05/14 10:51:16  abbott
// Summary: Added new field myNumDigits to MantExp2 and MantExp10
// Author: JAA
//
// Revision 1.3  2014/05/13 11:13:02  abbott
// Summary: MantissaAndExponent2 now accepts NumBits from 2 onwards
// Author: JAA
//
// Revision 1.2  2014/04/11 13:33:03  abbott
// Summary: Added MantissaAndExponent2 and MantissaAndExponent10
// Author: JAA
//
// Revision 1.1  2014/04/10 15:32:10  abbott
// Summary: New fn FloatApprox
// Author: JAA
//
//
