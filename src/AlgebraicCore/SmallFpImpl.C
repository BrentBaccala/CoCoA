//   Copyright (c)  2004-2007,2009,2011,2013  John Abbott

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

#include "CoCoA/SmallFpImpl.H"

#include "CoCoA/BigRat.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/NumTheory.H"
#include "CoCoA/error.H"

#include <cmath>
using std::floor;
using std::sqrt;
#include <limits>
using std::numeric_limits;

namespace CoCoA
{

  namespace // anonymous namespace
  {
    // JAA believes that this routine always returns a reduced value,
    // but does not have a proof of this currently.
    SmallFpImpl::value_t exgcd(SmallFpImpl::value_t r, SmallFpImpl::value_t m)
    {
      CoCoA_ASSERT(r != 0);  CoCoA_ASSERT(m != 0);
      SmallFpImpl::value_t p = m;
      SmallFpImpl::value_t cr = 1;
      SmallFpImpl::value_t cm = 0; // this is minus what you think it is!
      while (r != 0)
      {
	SmallFpImpl::value_t q = m/r;
	m -= q*r;
	cm += q*cr;
	if (m == 0) break;
	q = r/m;
	r -= q*m;
	cr += q*cm;
      }
      //???  if (r+m != 1) error("cannot invert zero-divisor");
      if (r == 0) return p-cm;
      return cr;
    }


  } // end of anonymous namespace


  // These two inline fns are used only in the ctors.
  inline SmallFpImpl::value_t SmallFpImpl::ourCalcResidueUPB(value_t p)
  {
    const value_t MAX = numeric_limits<value_t>::max();
    return p*(MAX/p/2); // Largest multiple of p <= MAX/2; exploits integer division.
  }

  inline long SmallFpImpl::ourCalcIterLimit(unsigned long p)
  {
    const value_t MAX = numeric_limits<value_t>::max();
    const value_t MaxIters = (MAX/(p-1))/(p-1)/2; // Max no. of unreduced products you can sum without exceeding MAX/2.
    const unsigned long MaxLong = numeric_limits<long>::max();
    if (MaxIters > MaxLong) return MaxLong; // JAA reckons this'll never happen.
    return MaxIters;
  }


  SmallFpImpl::SmallFpImpl(const MachineInt& p, GlobalSettings::ResidueSetting ResidueChoice):
      myModulusValue(ourCheckCtorArg(p)),
      myResiduesAreSymm(ResidueChoice==GlobalSettings::SymmResidues),
      myResidueUPBValue(ourCalcResidueUPB(myModulusValue)),
      myIterLimit(ourCalcIterLimit(myModulusValue))
  {}


  bool SmallFpImpl::IsGoodCtorArg(const MachineInt& n)
  {
    if (IsNegative(n) || !IsSignedLong(n)) return false;
    const long N = AsSignedLong(n);
    return N <= ourMaxModulus() && IsPrime(N);
  }

  long SmallFpImpl::ourMaxModulus()
  {
    const double HalfMaxIntVal = numeric_limits<value_t>::max()/2; // may harmlessly round up
    const long ans = PrevPrime(static_cast<long>(std::floor(sqrt(HalfMaxIntVal))));
    CoCoA_ASSERT(ourCalcIterLimit(ans) > 0); // check that 2*ans^2 <= MAXLONG
    return ans;
  }


  SmallFpImpl::value_t SmallFpImpl::myReduce(const MachineInt& n) const
  {
    const value_t ans =  abs(n)%myModulusValue;
    if (!IsNegative(n) || ans == 0) return ans;
    return myModulusValue - ans;
  }

  SmallFpImpl::value_t SmallFpImpl::myReduce(const BigInt& N) const
  {
    return mpz_fdiv_ui(mpzref(N), myModulusValue);
  }

  SmallFpImpl::value_t SmallFpImpl::myReduce(const BigRat& q) const
  {
    const value_t D = mpz_fdiv_ui(mpq_denref(mpqref(q)), myModulusValue);
    CoCoA_ASSERT(D != 0);
    const value_t N = mpz_fdiv_ui(mpq_numref(mpqref(q)), myModulusValue);
    return myDiv(N, D);
  }


  SmallFpImpl::value_t SmallFpImpl::myDiv(value_t x, value_t y) const
  {
    CoCoA_ASSERT(x == myNormalize(x));
    CoCoA_ASSERT(y == myNormalize(y));
    CoCoA_ASSERT(y != 0);
    const value_t yrecip = exgcd(y, myModulusValue);
    return myMul(x, yrecip);
  }


  SmallFpImpl::value_t SmallFpImpl::myPower(value_t x, long n) const
  {
    CoCoA_ASSERT(x == myNormalize(x));
    CoCoA_ASSERT(x != 0 || n > 0); // complain about any non-positive power of 0
    if (x == 0) { return 0; }
    if (x == 1) { return 1; }
    n %= myModulusValue - 1; // OK by fermat's little theorem.
    if (n < 0) n += myModulusValue-1;
    if (n == 0) { return 1; }
    if (n == 1) { return x; }
    // Here we know that n >= 2 and x is not 0 or 1.
    // Below is an iterative version of repeated squaring.
    unsigned long mask = 1;
    unsigned long quartern = n/4;
    while (mask <= quartern) mask <<= 1;
    value_t ans = x;
    for (; mask != 0; mask >>= 1)
    {
      ans = myMul(ans, ans);
      if (n & mask) ans = myMul(ans, x);
    }
    return ans;
  }



  // If p is a small prime, return p as a value_t.
  // Otherwise throw an exception.
  SmallFpImpl::value_t SmallFpImpl::ourCheckCtorArg(const MachineInt& n)
  {
    if (!IsGoodCtorArg(n))
      CoCoA_ERROR(ERR::BadSmallFpChar, "SmallFpImpl ctor");
    return AsUnsignedLong(n);
  }


  std::ostream& operator<<(std::ostream& out, const SmallFpImpl& arith)
  {
    return out << "SmallFpImpl(" << arith.myModulus() << ")";
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/SmallFpImpl.C,v 1.25 2014/08/26 12:56:36 abbott Exp $
// $Log: SmallFpImpl.C,v $
// Revision 1.25  2014/08/26 12:56:36  abbott
// Summary: Corrected a comment
// Author: JAA
//
// Revision 1.24  2014/05/06 13:13:57  abbott
// Summary: Removed unused fn NumBits
// Author: JAA
//
// Revision 1.23  2014/05/02 13:56:26  abbott
// Summary: Added explanatopry comment about a CoCoA_ASSERT
// Author: JAA
//
// Revision 1.22  2014/04/30 16:12:57  abbott
// Summary: Changes arg of sqrt to be of type double
// Author: JAA
//
// Revision 1.21  2013/05/27 13:01:23  abbott
// Simplified IsGoodCtorArg (following hint on redmine).
// Some minor cosmetic changes.
//
// Revision 1.20  2013/04/29 09:05:24  abbott
// Changed local variable to unsigned to avoid compiler warning.
//
// Revision 1.19  2013/03/25 17:04:19  abbott
// Major clean-up of interface to SmallFpImpl/SmallFpLogImpl/SmallFpDoubleImpl
// (underlying impl remains much the same).  Removed lots of cruft.
// Consequential changes to RingFp* classes; small change to SparsePolyRing.
//
// Revision 1.18  2012/11/23 17:29:47  abbott
// Changed names of a private data field and a private mem fn.
//
// Revision 1.17  2012/09/07 15:21:13  abbott
// First stage of revision of SmallFpImpl interface (and SmallFpLog, SmallFpDouble).
//
// Revision 1.16  2012/05/28 09:18:20  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.15  2012/01/30 23:27:57  abbott
// Added print function.
//
// Revision 1.14  2011/11/09 14:29:37  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.13  2011/08/24 10:29:55  bigatti
// -- renamed QQ --> BigRat
// -- sorted #include
//
// Revision 1.12  2011/08/14 15:52:16  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.11  2011/05/20 19:26:05  abbott
// Updated SmallFp*Impl: removed all output-related fns (must use myExport instead).
//
// Revision 1.10  2011/05/19 14:38:27  abbott
// Updated small prime finite field impls to allow user to specify
// separately for each whether to use symmetric or non-negative
// residues for export operations (myExport and printing).
//
// Revision 1.9  2011/05/03 10:04:35  abbott
// Hacked code slightly to shut the *!@% compiler up  >-{
//
// Revision 1.8  2011/03/22 22:38:41  abbott
// Corrected copyright years.
//
// Revision 1.7  2011/03/22 20:06:13  abbott
// Added static mem fn IsGoodCtorArg (called by RingFp pseudo-ctors).
// Commented out ctors which take ZZ arg -- seems useless.
//
// Revision 1.6  2011/03/14 10:28:14  abbott
// Changed unsigned long into long (and unsigned short into short).
//
// Revision 1.5  2009/12/11 11:46:32  abbott
// Changed fn  convert  into  IsConvertible.
// Added template procedure  convert.
// New version because change is not backward compatible.
//
// Revision 1.4  2009/05/14 09:39:29  abbott
// Added possibility to specify "symmetric" or "non-negative" residues
// in quotients of ZZ.  Affects printing of elements in quotients of ZZ
// (also changed printing of elements in general quotient rings).
// Consequent changes in several tests.
//
// Revision 1.3  2007/10/30 17:14:07  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.2  2007/03/23 18:38:42  abbott
// Separated the "convert" functions (and CheckedCast) into their own files.
// Many consequential changes.  Also corrected conversion to doubles.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.4  2007/03/08 18:22:29  cocoa
// Just whitespace cleaning.
//
// Revision 1.3  2006/11/27 13:06:22  cocoa
// Anna and Michael made me check without writing a proper message.
//
// Revision 1.2  2006/11/02 13:25:43  cocoa
// Simplification of header files: the OpenMath classes have been renamed.
// Many minor consequential changes.
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.3  2006/05/12 16:10:58  cocoa
// Added OpenMathFwd.H, and tidied OpenMath.H.
// Many consequential but trivial changes.
//
// Revision 1.2  2006/03/12 21:28:33  cocoa
// Major check in after many changes
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.4  2005/10/14 15:25:07  cocoa
// Major tidying and cleaning to small prime finite fields.
// Several consequential changes.  Improved their documentation.
//
// Added Makefile and script to include/CoCoA/ directory to
// keep library.H up to date.
//
// Revision 1.3  2005/10/12 15:52:09  cocoa
// Completed test-RingFp1 and corrected/cleaned the SmallFp*
// and RingFp* files.
//
// Some minor tidying elsewhere.
//
// Revision 1.2  2005/10/11 16:37:30  cocoa
// Added new small prime finite field class (see RingFpDouble).
//
// Cleaned makefiles and configuration script.
//
// Tidied PPMonoid code (to eliminate compiler warnings).
//
// Fixed bug in RingFloat::myIsInteger.
//
// Revision 1.1.1.1  2005/05/03 15:47:31  cocoa
// Imported files
//
// Revision 1.3  2005/04/20 15:40:48  cocoa
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
// Revision 1.4  2004/11/12 15:49:29  cocoa
// Tidying prior to 0.90 release.
// (a) documentation improved (or marked as poor)
// (b) sundry minor improvements to the code
//
// Revision 1.3  2004/07/20 15:04:06  cocoa
// The next step in the new "ring element" conversion process:
// handling the case of creating a "const RefRingElem" object
// (since C++ refuses to do this properly itself).
//
// Revision 1.2  2004/07/14 16:40:42  cocoa
// Separated RingFpLog from its implementation which now resides in
// a new class: SmallFpLogImpl.  This is analogous to the change made
// to RingFp yesterday.
//
// Some tidying and other sundry minor changes.
//
// Revision 1.1  2004/07/13 16:32:26  cocoa
// First stage of major revamp of ring elements.
// Implementation of RingFp has been split into "ring interface"
// and "algorithms plus data structures".
//
//
