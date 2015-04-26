//   Copyright (c)  2003-2011  John Abbott

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


#include "CoCoA/BigInt.H"
#include "CoCoA/IntOperations.H"
#include "CoCoA/OpenMath.H"
#include "CoCoA/error.H"


#include <iostream>
using std::ostream;
using std::istream;
using std::ios;

#include <string>
using std::string;

#include <vector>
using std::vector;

namespace CoCoA
{

  BigInt::BigInt()
  {
    mpz_init(myRepr);
  }


  // CAREFUL!  BigInt(0) goes to another ctor below (from mpz_t)
  BigInt::BigInt(const MachineInt& n)
  {
    if (IsNegative(n))
      mpz_init_set_si(myRepr, AsSignedLong(n));
    else
      mpz_init_set_ui(myRepr, AsUnsignedLong(n));
  }


  BigInt::BigInt(const std::string& str)
  {
    mpz_init(myRepr);
//     if (base != 0 && (base < 2 || base > 36))
//       CoCoA_ERROR(ERR::BadNumBase, "BigInt(string)");
    if (mpz_set_str(myRepr, str.c_str(), 10) != 0)
      CoCoA_ERROR(ERR::BadArg, "BigInt(string)");
  }


  // CAREFUL!!  BigInt(0) comes here!!
  BigInt::BigInt(const mpz_t N)
  {
    if (N == 0/*nullptr*/)
    {
      mpz_init(myRepr);
      return;
    }
    mpz_init_set(myRepr, N);
  }

//   BigInt::BigInt(MP_INT to_be_owned)
//   {
//     myRepr[0] = to_be_owned;
//   }


//   BigInt::BigInt(CopyFromMPZMarker, const mpz_t to_be_copied)
//   {
//     mpz_init_set(myRepr, to_be_copied);
//   }

// COPY CONSTRUCTOR

  BigInt::BigInt(const BigInt& from)
  {
    mpz_init_set(myRepr, from.myRepr);
  }


  BigInt& BigInt::operator=(const BigInt& rhs)
  {
    if (this == &rhs) return *this;
    mpz_set(myRepr, rhs.myRepr);
    return *this;
  }



  // BUG/SLUG slow and wastes space -- better to use ostringstream?
  std::string ConvertToString(const BigInt& src, int base)
  {
    if (base < 2 || base > 36)
      CoCoA_ERROR(ERR::BadNumBase, "ConvertToString(BigInt,int)");
    const long digits = NumDigits(src, base);
    string ans; ans.reserve(digits+1); // +1 to allow for minus sign
    vector<char> buffer(digits+2); // +2 to allow for minus sign and terminating NUL
    mpz_get_str(&buffer[0], base, mpzref(src));
    ans = &buffer[0]; // Won't throw as space is already reserved.
    return ans;
  }


  ostream& operator<<(ostream& out, const BigInt& N)
  {
//    Using GMPXX next two lines should do it all...
//    mpz_class Ncopy(N.myRepr);
//    return out << Ncopy;

    // Dispose of the exceptional case n=0 immediately.
    if (IsZero(N)) return out << '0';
    const int base = (out.flags() & ios::oct) ? 8 :
                     (out.flags() & ios::hex) ? 16 : 10;

    // prefix will contain sign and/or base specification
    string prefix;

    // Firstly put the sign in the prefix...
    if (N < 0) prefix += '-';
    else if (out.flags() & ios::showpos) prefix += '+';

    // Next put base indication in the prefix if required...
    if (out.flags() & ios::showbase)
    {
      if (base == 8 || base == 16) prefix += '0';
      if (base == 16) prefix += 'x';
    }

    out << prefix << ConvertToString(abs(N), base);
    return out;
  }


  namespace  // anonymous namespace for file local functions
  {

    // Returns true if c is a valid digit for the given base
    // (both upper and lower cases are allowed for hexadecimal).
    // Note: base must be 8, 10 or 16 (otherwise always gives false)
    bool IsDigitBase(char c, int base) // base = 8, 10 or 16
    {
      switch (base)
      {
          case 10:
            return isdigit(c);
          case 16:
            return isxdigit(c);
          case 8:
            return (isdigit(c) &&  c != '8' && c != '9');
          default:
            return false;
      }
    }
  } // end of anonymous namespace

  istream& operator>>(istream& in, BigInt& N)
  {
    if (!in.good()) return in;
//     mpz_class buf;
//     in >> buf;
//     mpz_set(N.myRepr, buf.get_mpz_t());

    const int base = (in.flags() & ios::oct) ? 8 : (in.flags() & ios::hex) ? 16 : 10;
    ws(in);

    // Look for sign of number.
    const char FirstChar = in.peek(); // this may set eofbit
    if (!in.good()) return in;
    if (FirstChar == '-' || FirstChar == '+')
      in.ignore();

    // Read in as many digits as possible.
    string digits;
    while (true)
    {
      const char ch = in.peek(); // this may set eofbit
      if (!in.good() || !IsDigitBase(ch, base)) break;
      in.ignore();
      digits += ch;
    }
    if (digits.empty())
    {
      if (FirstChar == '-' || FirstChar == '+')
        in.putback(FirstChar);
      in.setstate(ios::failbit); // set failbit (perhaps in addition to eofbit)
      return in;
    }
    // We found some digits, so convert them into a number.
    N = BigInt(digits); // could throw if number is huge.
    if (FirstChar == '-') negate(N);
    return in;
  }


  OpenMathOutput& operator<<(OpenMathOutput& OMOut, const BigInt& N)
  {
    OMOut->mySend(N);
    return OMOut;
  }


  OpenMathInput& operator>>(OpenMathInput& OMIn, BigInt& /*N*/)
  {
    CoCoA_ERROR(ERR::NYI, "OpenMathInput fns");
    return OMIn;
  }


  //---------------------------------------------------------------------------
  // Assignment and assignment arithmetic functions

  BigInt& BigInt::operator+=(const BigInt& rhs)
  {
    mpz_add(myRepr, myRepr, rhs.myRepr);
    return *this;
  }


  BigInt& BigInt::operator-=(const BigInt& rhs)
  {
    mpz_sub(myRepr, myRepr, rhs.myRepr);
    return *this;
  }


  BigInt& BigInt::operator*=(const BigInt& rhs)
  {
    mpz_mul(myRepr, myRepr, rhs.myRepr);
    return *this;
  }


  BigInt& BigInt::operator/=(const BigInt& rhs)
  {
    if (IsZero(rhs))
      CoCoA_ERROR(ERR::DivByZero, "BigInt /= BigInt");
//???    if (rhs < 0 && !mpz_divisible_p(myRepr, mpzref(rhs)))
//???      CoCoA_ERROR(ERR::IntDivByNeg, "BigInt /= BigInt");
    mpz_tdiv_q(myRepr, myRepr, rhs.myRepr);
    return *this;
  }


  BigInt& BigInt::operator%=(const BigInt& rhs)
  {
    if (IsZero(rhs))
      CoCoA_ERROR(ERR::ZeroModulus, "BigInt %= BigInt");
//???    if (rhs < 0 && !mpz_divisible_p(myRepr, mpzref(rhs)))
//???      CoCoA_ERROR(ERR::IntDivByNeg, "BigInt %= BigInt");
    mpz_tdiv_r(myRepr, myRepr, rhs.myRepr);
    return *this;
  }


  //---------------------------------------------------------------------------
  // Assignment and assignment arithmetic with rhs a MachineInt

  BigInt& BigInt::operator= (const MachineInt& rhs)
  {
    if (IsNegative(rhs))
      mpz_set_si(myRepr, AsSignedLong(rhs));
    else
      mpz_set_ui(myRepr, AsUnsignedLong(rhs));

    return *this;
  }


  BigInt& BigInt::operator+=(const MachineInt& rhs)
  {
    if (IsNegative(rhs))
      mpz_sub_ui(myRepr, myRepr, abs(rhs));
    else
      mpz_add_ui(myRepr, myRepr, AsUnsignedLong(rhs));
    return *this;
  }


  BigInt& BigInt::operator-=(const MachineInt& rhs)
  {
    if (IsNegative(rhs))
      mpz_add_ui(myRepr, myRepr, abs(rhs));
    else
      mpz_sub_ui(myRepr, myRepr, AsUnsignedLong(rhs));
    return *this;
  }


  BigInt& BigInt::operator*=(const MachineInt& rhs)
  {
    if (IsNegative(rhs))
      mpz_mul_si(myRepr, myRepr, AsSignedLong(rhs));
    else
      mpz_mul_ui(myRepr, myRepr, AsUnsignedLong(rhs));
    return *this;
  }


  BigInt& BigInt::operator/=(const MachineInt& rhs)
  {
    if (IsZero(rhs))
      CoCoA_ERROR(ERR::DivByZero, "BigInt /= MachineInt");
    if (IsNegative(rhs))
    {
//???      if (!mpz_divisible_ui_p(myRepr, abs(rhs)))
//???        CoCoA_ERROR(ERR::IntDivByNeg, "BigInt /= MachineInt");
      mpz_neg(myRepr, myRepr);
      mpz_tdiv_q_ui(myRepr, myRepr, abs(rhs));
    }
    else
      mpz_tdiv_q_ui(myRepr, myRepr, AsUnsignedLong(rhs));
    return *this;
  }


  BigInt& BigInt::operator%=(const MachineInt& rhs)
  {
    if (IsZero(rhs))
      CoCoA_ERROR(ERR::ZeroModulus, "BigInt %= MachineInt");
    if (IsNegative(rhs))
    {
//???      if (!mpz_divisible_ui_p(myRepr, abs(rhs)))
//???        CoCoA_ERROR(ERR::IntDivByNeg, "BigInt %= MachineInt");
      mpz_tdiv_r_ui(myRepr, myRepr, abs(rhs));
    }
    else
      mpz_tdiv_r_ui(myRepr, myRepr, AsUnsignedLong(rhs));
    return *this;
  }



  //---------------------------------------------------------------------------
  // increment and decrement

  BigInt& BigInt::operator++()
  {
    mpz_add_ui(myRepr, myRepr, 1);
    return *this;
  }


  BigInt& BigInt::operator--()
  {
    mpz_sub_ui(myRepr, myRepr, 1);
    return *this;
  }


  const BigInt BigInt::operator++(int)  // INEFFICIENT
  {
    BigInt ans(*this);
    ++*this;
    return ans;
  }


  const BigInt BigInt::operator--(int)  // INEFFICIENT
  {
    BigInt ans(*this);
    --*this;
    return ans;
  }


  //---------------------------------------------------------------------------

} // end of namespace CoCoA



// The next few lines contain RCS header/log information.
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/BigInt.C,v 1.14 2013/03/26 14:56:06 abbott Exp $
// $Log: BigInt.C,v $
// Revision 1.14  2013/03/26 14:56:06  abbott
// Updated the conversion fns (in ptic removed procedure "convert");
// numerous consequential changes.
//
// Revision 1.13  2013/02/14 15:35:48  abbott
// Added a cosmetic space.
//
// Revision 1.12  2012/05/28 09:18:21  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.11  2012/05/25 13:01:23  abbott
// Added fn IsDivisible.
//
// Revision 1.10  2012/02/07 10:47:12  bigatti
// -- updated comment to "BigInt"
//
// Revision 1.9  2011/11/09 14:03:40  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.8  2011/09/13 15:49:20  abbott
// Added "static" to log2 local var in log(BigInt).
//
// Revision 1.7  2011/09/06 15:21:53  abbott
// Changed "cmp" functions so that the return value is in {-1,0,+1}.
//
// Revision 1.6  2011/08/25 06:30:08  bigatti
// -- added log from ZZ.C
//
// Revision 1.5  2011/08/23 16:16:29  abbott
// Corrected defn of RoundDiv; corrected comment too.
//
// Revision 1.4  2011/08/17 11:56:57  abbott
// Added static_cast to keep compiler quiet.
//
// Revision 1.3  2011/08/14 15:52:17  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.2  2011/08/12 16:31:03  abbott
// COMMENTED OUT SOME FNS SO THAT BigInt CAN EXIST ALONGSIDE ZZ
// FOR THE TIME BEING.
//
// Revision 1.1  2011/08/12 15:21:26  abbott
// Added BigInt impl (derived from ZZ); not used by anyone yet.
//
//------ log with old name ZZ.C -------------------------------
// Revision 1.28  2011/03/21 11:13:25  abbott
// Changed return type for operator%(ZZ,MachineInteger) from
// unsigned long to long to fit in with the new coding conventions.
//
// Revision 1.27  2011/03/14 10:32:22  abbott
// Changed size_t into long (return types of NumDigits & exponent).
//
// Revision 1.26  2011/03/01 15:27:23  abbott
// Improved impl of ILogBase -- faster in most cases.
// Changed the error objects thrown in certain cases (e.g. ERR::BadArg --> ERR::LogZero).
//
// Revision 1.25  2011/01/14 17:21:15  abbott
// Added isqrt, iroot, IsExactIroot, IsSquare, IsPower.
//
// Revision 1.24  2010/12/26 13:03:15  abbott
// Added ILogBase function (to BigInt & BigRat).
//
// Revision 1.23  2010/04/23 13:54:44  abbott
// binomial function now triggers error if 2nd arg is negative.
//
// Revision 1.22  2010/03/23 11:59:35  abbott
// Cleaned mantissa -- GMP already has a fn for this!
// Removed use of shared static "temporary" ZZ value -- caused GMPAllocator to report an error!
//
// Revision 1.21  2010/03/22 21:01:31  abbott
// Cleaned impl of binomial (& fixed a couple of silly bugs).
//
// Revision 1.20  2010/03/22 11:50:31  abbott
// Added ctor from a string.
// Fixed stupid bug in operator-.
//
// Revision 1.19  2010/03/18 16:41:08  abbott
// Commented out unused arg (in NYI fn).
//
// Revision 1.18  2010/03/18 13:54:19  abbott
// Added openmath output fns (moved from OpenMath files).
//
// Revision 1.17  2009/12/29 22:44:32  abbott
// Removed buggy proxy class ZZ::rtn.
// Consequent changes for function prototypes also in NumTheory.
// Corrected some minor buglets in NumTheory.
//
// Revision 1.16  2009/12/11 11:42:16  abbott
// Changed convert into IsConvertible.
// Cleaned impl of operator>>(istream&,ZZ&).
//
// Revision 1.15  2009/11/26 16:18:00  bigatti
// -- including string ZZ.C instead of ZZ.H
//
// Revision 1.14  2009/10/08 13:39:47  abbott
// Renamed "round" into "RoundDiv".
// Added some new versions of "RoundDiv".
// Added a test for "RoundDiv".
//
// Revision 1.13  2009/06/11 14:05:29  abbott
// CLEANING: Removed several functions which are now gathered in NumTheory.H/C
//           (for example: gcd, lcm, PowerMod, InvMod).
//
// Revision 1.12  2009/06/05 12:08:27  abbott
// Changed return type of operator%(ZZ,MachineInteger); it is now unsigned long
// instead of ZZ.
//
// Revision 1.11  2009/01/12 14:34:10  abbott
// Corrected IsDigitBase (following a bug report from anonymous)
// for octal and hexadecimal digits.  Added 3 lines of comment too.
//
// Revision 1.10  2008/12/16 21:14:47  abbott
// In functions taking a machine integer changed arg type from MachineInteger
// to const-ref-MachineInteger.
//
// Revision 1.9  2008/11/18 10:25:48  abbott
// Added function round.
//
// Revision 1.8  2008/04/22 13:09:16  abbott
// Removed IsPositive and IsNegative functions for ZZs.
// Comparison between RingElem and 0 is now handled specially (specially fast?).
//
// Revision 1.7  2007/10/30 17:14:07  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.6  2007/06/01 13:56:08  abbott
// Now gives error whenever you try to compute 0^0 -- previously some cases were not checked.
//
// Revision 1.5  2007/05/22 22:51:39  abbott
// Changed name of fn ndigits to NumDigits.
// Changed return type of exponent and NumDigits.
// Changed some exceptions from nonstandard to the appropriate standard one.
//
// Revision 1.4  2007/05/21 12:57:28  abbott
// New class for passing machine integers as args; includes some simple
// operations on machine integers (cmp, gcd, IsNegative,...).  Operations
// between ZZ and machine integers modified to use the new class.  Inexact
// integer division (of a ZZ) by a negative value now triggers an error;
// new error for reporting inexact integer division by a negative value.
//
// Revision 1.3  2007/03/23 18:38:42  abbott
// Separated the "convert" functions (and CheckedCast) into their own files.
// Many consequential changes.  Also corrected conversion to doubles.
//
// Revision 1.2  2007/03/16 17:43:05  abbott
// Added new convert function (from ZZ to double).
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.9  2007/03/08 18:22:28  cocoa
// Just whitespace cleaning.
//
// Revision 1.8  2007/01/13 14:14:34  cocoa
// Overhaul of RingHom code: it nows uses SmartPtrIRC, and printing is more logical.
// Have not yet updated the documentation.
//
// Revision 1.7  2007/01/09 15:52:08  cocoa
// Changed QBGenerator to use std::vector instead of std::list for the result.
// Minor mod to configure script.
//
// Revision 1.6  2006/11/29 11:59:35  cocoa
// -- fixed: convert(double& z, const ZZ& num, const ZZ& den) now returns
//    bool (was void) and does not throw errors
//
// Revision 1.5  2006/11/27 13:06:22  cocoa
// Anna and Michael made me check without writing a proper message.
//
// Revision 1.4  2006/10/16 23:18:59  cocoa
// Corrected use of std::swap and various special swap functions.
// Improved myApply memfn for homs of RingDistrMPolyInlPP.
//
// Revision 1.3  2006/10/06 09:36:13  cocoa
// Fixed two minor bugs introduced at last check in.
//
// Revision 1.2  2006/09/29 11:46:54  cocoa
// Corrected bug in convert(ZZ, ZZ, double) -- now correct and simpler.
// Previously went into infinite loop on negative doubles.
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.5  2006/05/12 16:10:58  cocoa
// Added OpenMathFwd.H, and tidied OpenMath.H.
// Many consequential but trivial changes.
//
// Revision 1.4  2006/03/27 12:21:25  cocoa
// Minor silly changes to reduce number of complaints from some compiler or other.
//
// Revision 1.3  2006/03/14 15:01:49  cocoa
// Improved the implementation of ring member fns for computing powers.
// Should keep Intel C++ compiler quieter too.
//
// Revision 1.2  2006/03/12 21:28:33  cocoa
// Major check in after many changes
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.2  2005/08/08 16:36:32  cocoa
// Just checking in before going on holiday.
// Don't really recall what changes have been made.
// Added IsIndet function for RingElem, PPMonoidElem,
// and a member function of OrdvArith.
// Improved the way failed assertions are handled.
//
// Revision 1.1.1.1  2005/05/03 15:47:31  cocoa
// Imported files
//
// Revision 1.3  2005/04/20 15:40:47  cocoa
// Major change: modified the standard way errors are to be signalled
// (now via a macro which records filename and line number).  Updated
// documentation in error.txt accordingly.
//
// Improved the documentation in matrix.txt (still more work to be done).
//
// Revision 1.2  2005/04/19 14:06:03  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.8  2004/11/12 15:49:29  cocoa
// Tidying prior to 0.90 release.
// (a) documentation improved (or marked as poor)
// (b) sundry minor improvements to the code
//
// Revision 1.7  2004/11/05 15:36:57  cocoa
// Removed call to new(nothrow) inside a convert function.  Replaced
// it by use of a std::vector -- code is certainly cleaner and clearer
// now, but possibly with an inconsequential reduction in run-time
// performance.
//
// Revision 1.6  2004/11/04 18:47:42  cocoa
// (1) Ring member functions which previously expected mpz_t args
//     now expect ZZ args.  Numerous minor consequential changes.
// (2) Renamed function which gives access to the mpz_t value inside
//     a ZZ object: previously was raw(...), now is mpzref(...).
//     Plenty of calls had to be altered.
//
// Revision 1.5  2004/05/27 16:14:02  cocoa
// Minor revision for new coding conventions.
//
// Revision 1.4  2004/05/24 15:52:13  cocoa
// Major update:
//   new error mechanism
//   many fixes
//   RingHoms almost work now
//   RingFloat much improved
//
// Revision 1.3  2004/02/03 16:16:20  cocoa
// Removed pointless IamGCDDomain functions from several concrete rings.
// Added IamOrderedDomain functions where appropriate.
// Tidied ctors for the small finite fields.
//
// Revision 1.2  2003/10/17 10:51:06  cocoa
// Major cleaning, and new naming convention.
//
// Revision 1.1.1.1  2003/09/24 12:55:43  cocoa
// Imported files
//
// Revision 1.2  2003/06/23 16:13:01  abbott
// Minor cleaning prior to public release.
//
// Revision 1.1  2003/05/14 17:12:40  abbott
// Initial revision
//
