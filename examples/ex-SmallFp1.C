// Copyright (c) 2013  John Abbott
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"
using namespace CoCoA;
using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "**ADVANCED**  This example program is for advanced CoCoALib users. \n"
  "It shows how to use SmallFpImpl for efficient arithmetic in a small\n"
  "prime finite field.  Not for the faint-hearted!                    \n";

const string LongDescription =
  "**ADVANCED**  This example program is for advanced CoCoALib users.  \n"
  "SmallFpImpl enables you to perform arithmetic efficiently in a small\n"
  "prime finite field.  The catch is that correct use is not as simple \n"
  "as using RingElems directly.  We take as a specific illustrative    \n"
  "example the computation of an inner product.                        \n";
//----------------------------------------------------------------------

// The impl of inner product for RingElem is pretty simple and clear:
RingElem InnerProd(const vector<RingElem>& u, const vector<RingElem>& v)
{
  ring Fp = owner(u[0]);
  const long n = len(u);
  RingElem ans(Fp);
  for (long i=0; i < n; ++i)
    ans += u[i]*v[i];
  return ans;
}

void TimeTrialRingElem(long p)
{
  ring Fp = NewZZmod(p);
  vector<RingElem> u; u.reserve(p*p);
  vector<RingElem> v; v.reserve(p*p);
  for (long i=0; i < p; ++i)
    for (long j=0; j < p; ++j)
    {
      u.push_back(RingElem(Fp,i));
      v.push_back(RingElem(Fp,j));
    }

  // Timing the computation of the inner product:
  const double StartTime = CpuTime();
  RingElem InProd = InnerProd(u, v);
  const double EndTime = CpuTime();
  cout << "Ans is " << InProd << endl;
  cout << "Using ring ZZmod(" << p << ") time is " << EndTime-StartTime << endl;
}


// We present 2 fast impls:
//   (A) is slightly clearer, while
//   (B) is slightly faster.


// Handy typedef to make reading/writing the code simpler!!!
typedef SmallFpImpl::value_t FpElem;


// Impl (A) for "fast" inner product;
// it is much fiddlier than the RingElem implementation above!
FpElem InnerProd_A(const SmallFpImpl& Fp, const vector<FpElem>& u, const vector<FpElem>& v)
{
  const long n = len(u);
  const long MaxSafeIters = Fp.myMaxIters();
  long IterCount = 0;
  FpElem ans = 0;
  for (long i=0; i < n; ++i)
  {
    ans += u[i]*v[i];  // <--- the actual computation
    ++IterCount;
    if (IterCount < MaxSafeIters) continue;
    IterCount = 0;
    ans = Fp.myHalfNormalize(ans);
  }
  return Fp.myNormalize(ans);
}

// Impl (B) for "fast" inner product;
// it is harder to understand than (A), but is a bit faster (on my computer).
FpElem InnerProd_B(const SmallFpImpl& Fp, const vector<FpElem>& u, const vector<FpElem>& v)
{
  const long n = len(u);
  const long MaxSafeIters = Fp.myMaxIters();
  long i = 0;
  long NextNormalize = 0;
  FpElem ans = 0;
  while (NextNormalize < n)
  {
    NextNormalize += MaxSafeIters;
    if (NextNormalize > n) NextNormalize = n;
    for (; i < NextNormalize; ++i)
    {
      ans += u[i]*v[i];  // <--- the actual computation
    }
    ans = Fp.myHalfNormalize(ans);
  }
  return Fp.myNormalize(ans);
}


// Just for comparison here is an *inefficient* impl using SmallFpImpl
// (every intermediate value is normalized).
FpElem InnerProd_slow(const SmallFpImpl& Fp, const vector<FpElem>& u, const vector<FpElem>& v)
{
  const long n = len(u);
  FpElem ans = 0;
  for (long i=0; i < n; ++i)
    ans = Fp.myAdd(ans, Fp.myMul(u[i],v[i]));
  return ans;
}


void TimeTrialSmallFp(long p)
{
  if (!SmallFpImpl::IsGoodCtorArg(p)) return;//????
  SmallFpImpl Fp(p);
  // Create two vectors to work on
  vector<FpElem> u; u.reserve(p*p);
  vector<FpElem> v; v.reserve(p*p);
  for (long i=0; i < p; ++i)
    for (long j=0; j < p; ++j)
    {
      u.push_back(Fp.myReduce(i));
      v.push_back(Fp.myReduce(j));
    }

  // Timing method (A)
  const double StartTime_A = CpuTime();
  const FpElem InProd_A = InnerProd_A(Fp, u, v);
  const double EndTime_A = CpuTime();
  cout << "Ans is " << InProd_A << endl;
  cout << "Using impl (A) for p=" << p << "  time is " << (EndTime_A - StartTime_A) << endl;

  // Timing method (B)
  const double StartTime_B = CpuTime();
  const FpElem InProd_B = InnerProd_B(Fp, u, v);
  const double EndTime_B = CpuTime();
  cout << "Ans is " << InProd_B << endl;
  cout << "Using impl (B) for p=" << p << "  time is " << (EndTime_B - StartTime_B) << endl;

  // Let's see how much it costs to normalize every result
  const double StartTime_slow = CpuTime();
  const FpElem InProd_slow = InnerProd_slow(Fp, u, v);
  const double EndTime_slow = CpuTime();
  cout << "Ans is " << InProd_slow << endl;
  cout << "Using impl (slow) for p=" << p << "  time is " << (EndTime_slow - StartTime_slow) << endl;
}


void program()
{
  GlobalManager CoCoAFoundations;

  cout << ShortDescription << endl;

  const int p = 2237;
  TimeTrialRingElem(p);
  TimeTrialSmallFp(p);
}

//----------------------------------------------------------------------
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
    cerr << "***ERROR***  UNCAUGHT CoCoA error";
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
  return 1;
}

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-SmallFp1.C,v 1.3 2014/04/03 15:43:07 abbott Exp $
// $Log: ex-SmallFp1.C,v $
// Revision 1.3  2014/04/03 15:43:07  abbott
// Summary: Reduced size of prime p (o/w too slow with debugging on some machines)
// Author: JAA
//
// Revision 1.2  2013/05/27 14:48:18  abbott
// Added typedef for FpElem to make code more readable.
//
// Revision 1.1  2013/05/27 12:55:04  abbott
// Added new example ex-SmallFp1.C
//
// Revision 1.6  2012/11/30 14:04:55  abbott
// Increased visibility of comment saying "put your code here".
//
// Revision 1.5  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.4  2008/10/07 12:12:54  abbott
// Removed useless commented out #include.
//
// Revision 1.3  2007/05/31 16:06:16  bigatti
// -- removed previous unwanted checked-in version
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.9  2007/03/07 11:51:40  bigatti
// -- improved test alignment
//
// Revision 1.8  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.7  2007/03/02 17:46:40  bigatti
// -- unique RingZ and RingQ
// -- requires foundations.H ;  foundations blah;  (thik of a better name)
//
// Revision 1.6  2007/03/02 10:47:53  cocoa
// First stage of RingZ modifications -- tests do not compile currently, Anna will fix this.
//
// Revision 1.5  2007/03/01 13:52:59  bigatti
// -- minor: fixed typo
//
// Revision 1.4  2007/02/28 15:15:56  bigatti
// -- minor: removed quotes in description
//
// Revision 1.3  2007/02/12 16:27:43  bigatti
// -- added strings ShortDescription and LongDescription for indexing
//
// Revision 1.2  2007/02/10 18:44:03  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.1  2006/03/12 21:28:34  cocoa
// Major check in after many changes
//
