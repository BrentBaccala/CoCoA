#include "CoCoA/BuildInfo.H"
#include "CoCoA/DenseUPolyRing.H"
#include "CoCoA/FractionField.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingDenseUPolyClean.H"
#include "CoCoA/RingTwinFloat.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/VectorOperations.H"
#include "CoCoA/symbol.H"

using namespace CoCoA;

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <vector>
using std::vector;

//----------------------------------------------------------------------
// Test for DenseUPolyRing functions on RingElem
// functions: +, -f, f-g, *, StdDeg, gcd, TidyGens, %, IsIrred
// environments: DenseUPolyClean (Fp, Z)
//----------------------------------------------------------------------

#define TEST_ASSERT(cond)   do { if (!(cond)) CoCoA::AssertionFailed(#cond, __FILE__, __LINE__); } while(0)


void TestDenseUPolyRing(DenseUPolyRing P)
{
  cout << "TESTING: " << P << endl << endl;

  if (!IsExact(P))
  {
    RingElem tf(NewRingTwinFloat(1024), 4);
    RingElem four(P);
    P->myRecvTwinFloat(raw(four), raw(tf));
    TEST_ASSERT(four == 4);
  }

  const RingElem& x = indet(P, 0);
  RingElem f1 = 3*x +23,
    f2 = IndetPower(P, 0, 4) + x,  // x^4 + x
    f3 = x + 29,
    f4(P);

  P->myAssignNonZeroCoeff(raw(f1), raw(one(CoeffRing(P))), 2); // x^2

  cout << "Given f1 = " << f1 << endl;
  cout << "  and f2 = " << f2 << endl << endl;
  cout << "f1+f2   gives  " << f1+f2 << endl;
  cout << "f2+f1   gives  " << f2+f1 << endl;
  cout << "f2-f1   gives  " << f2-f1 << endl;
  cout << "-f1   gives  "   << -f1 << endl;
  cout << "deg((2*x+1)*(3*x+1)) gives  "   << deg((2*x+1)*(3*x+1)) << endl;
  //   //  cout << "gcd(f1,f2)   gives  " << gcd(f1,f2) << endl;
  cout << "StdDeg(f1)   gives  " << StdDeg(f1) << endl;
  cout << "StdDeg(f2)   gives  " << StdDeg(f2) << endl;

  // if (IsTrueGCDDomain(P)) TEST_ASSERT(gcd(f1*f1f2, f2*f1f2) == f1f2);
  RingElem f1f2 = f1*f2;
  RingElem f2f1 = f2*f1;
  if (IsField(CoeffRing(P)) || IsTrueGCDDomain(CoeffRing(P)))
  {
    TEST_ASSERT( gcd(f1,f2) == 1 );
    TEST_ASSERT( IsInvertible(gcd(f1f2, f1*f3)/f1) );
    TEST_ASSERT( IsInvertible(gcd(f1f2, f2*f3)/f2) );
    TEST_ASSERT( IsInvertible(gcd(f1*f3, f2*f3)/f3) );
    TEST_ASSERT( IsInvertible(gcd(f1*f1f2, f2f1*f2)/(f1f2)) );
    TEST_ASSERT( IsDivisible(f1f2, f1) );
    TEST_ASSERT( IsDivisible(f1f2, f2) );
    TEST_ASSERT( !IsDivisible(f1, f2) );
    TEST_ASSERT( !IsDivisible(f2, f1) );
    TEST_ASSERT( (f1f2)/f1 == f2 );
    TEST_ASSERT( (f1f2)/f2 == f1 );
    TEST_ASSERT( deriv(f1,x) == 2*x+3 );
    try { std::cout << deriv(f1,x+1) << std::endl; }
    catch (const CoCoA::ErrorInfo& err) { if (err != ERR::NotIndet) throw; }
  }
  if (IsField(CoeffRing(P)))
  {
    try
    {
      TEST_ASSERT( IsIrred(f1) );
      TEST_ASSERT( !IsIrred(f2) );
      TEST_ASSERT( !IsIrred(f1f2) );
    }
    catch (const CoCoA::ErrorInfo& err) { if (err != ERR::NYI) throw; }
    RingElem g(P), acof(P), bcof(P);
    P->myExgcd(raw(g), raw(acof), raw(bcof), raw(f1), raw(f2));
    TEST_ASSERT(IsInvertible(g) || (g == gcd(f1,f2)));
    TEST_ASSERT(g == acof*f1 + bcof*f2);
  }
  
  P->myMulByXExp(raw(f1), 5);
  cout << "  P->myMulByXExp(raw(f1), 5) gives " << f1 << endl;
  P->myMulBy1MinusXExp(raw(f1), 2);
  cout << "  P->myMulBy1MinusXExp(raw(f1), 2) gives " << f1 << endl;
  P->myAddMul(raw(f1), raw(LC(x)), 2, raw(f2));
  cout << "  P->myAddMul(raw(f1), raw(LC(x)), 2, raw(f2)) gives " << f1 << endl;
  
  cout << "one(P) = " << one(P) << endl;

  f1f2 = f1*f2;  // previous lines changed f1
  f2f1 = f2*f1;  // previous lines changed f1
  TEST_ASSERT(IsOne(one(P)));
  TEST_ASSERT(!IsOne(f1));
  TEST_ASSERT(!IsOne(f2));
  TEST_ASSERT(P->myIsValid(raw(f1)));
  TEST_ASSERT(P->myIsValid(raw(f2)));
  TEST_ASSERT(P->myIsValid(raw(f3)));
  TEST_ASSERT(P->myIsValid(raw(f4)));
  TEST_ASSERT(P->myIsValid(raw(f1f2)));
  TEST_ASSERT(P->myIsValid(raw(f2f1)));
  if (IsCommutative(P))  TEST_ASSERT(f1f2 == f2f1);
  TEST_ASSERT(f1 != f2);
  TEST_ASSERT(f1+f2 == f2+f1);
  TEST_ASSERT(f1+(-f1) == 0);
  TEST_ASSERT(f1-f1 == 0);
  TEST_ASSERT(f1f2 != f1);
  TEST_ASSERT(f1f2 != f2);
  
  if (IsField(CoeffRing(P)))
  {
    vector<RingElem> g;
    g.push_back(f1f2);
    g.push_back(f1*x);
    g.push_back(zero(P));
    g.push_back(f1*(x-1));
    ideal I = ideal(g);
    cout << "TidyGens(I) = " << TidyGens(I) << endl;
    TEST_ASSERT(!IsZero(I));
    TEST_ASSERT(IsZero(f1 % I));
    TEST_ASSERT(IsZero(f1f2 % I));
    TEST_ASSERT(!IsZero((f1f2+x) % I));
    
    ideal J1 = ideal(f1f2);
    ideal J2 = ideal(f1*(3*x-7));
    TEST_ASSERT(intersect(J1, J2) == ideal(f1f2*(3*x-7))); 
    TEST_ASSERT(colon(J1, J2) == ideal(f2)); 
 
    vector<RingElem> h(3, zero(P));
    TEST_ASSERT(TidyGens(ideal(h)).empty());
    TEST_ASSERT(IsZero(ideal(h)));
  }
  if ( IsFractionFieldOfGCDDomain(CoeffRing(P)) )
  {
    TEST_ASSERT( CommonDenom(x/3+one(P)/2) ==  6);
    TEST_ASSERT( ClearDenom(x/3+one(P)/2) ==  2*x+3);
  }
  TEST_ASSERT(IsZero(f4));
  f4 = 1;
  TEST_ASSERT(IsOne(f4));
  f2f1 = 1;
  TEST_ASSERT(P->myIsValid(raw(f2f1)));
  TEST_ASSERT(IsOne(f2f1));

  cout << "------------------------------------------------" << endl;
}


void program()
{
  GlobalManager CoCoAFoundations(UseNonNegResidues);

  const QuotientRing Fp = NewZZmod(101);
  const QuotientRing F6 = NewZZmod(6);
  const RingTwinFloat RR = NewRingTwinFloat(150);

  TestDenseUPolyRing(NewPolyRing_DUP(RingZZ()));
  TestDenseUPolyRing(NewPolyRing_DUP(Fp));
  TestDenseUPolyRing(NewPolyRing_DUP(F6, symbol("y")));
  TestDenseUPolyRing(NewPolyRing_DUP(RingQQ()));
  TestDenseUPolyRing(NewPolyRing_DUP(RR));
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
