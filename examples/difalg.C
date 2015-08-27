
#include "CoCoA/library.H"
// #include <unordered_map>
#include <algorithm>
#include <functional>
using namespace CoCoA;
using namespace std;

/* PPMonoidRingExpImpl
 *
 * A PPMonoid whose exponents are elements in a ring that is also an
 * ordered domain.  Only "positive" exponents are allowed (positive in
 * the sense that they compare greater than zero).
 */

class PPMonoidRingExpImpl;

class PPMonoidRingExpElem {
  friend PPMonoidRingExpImpl;
  std::vector<RingElem> exponents;

  PPMonoidRingExpElem(long numIndets) {
    exponents.resize(numIndets);
  }
};

class PPMonoidRingExpImpl: public PPMonoidBase
{
  typedef PPMonoidElemRawPtr RawPtr;           // just to save typing
  typedef PPMonoidElemConstRawPtr ConstRawPtr; // just to save typing

public:
  PPMonoidRingExpImpl(const std::vector<symbol>& IndetNames, const PPOrdering& ord, const ring& ExponentRing = RingZZ());
  ~PPMonoidRingExpImpl();
private: // disable copy ctor and assignment
  explicit PPMonoidRingExpImpl(const PPMonoidRingExpImpl& copy);  // NEVER DEFINED -- copy ctor disabled
  PPMonoidRingExpImpl& operator=(const PPMonoidRingExpImpl& rhs); // NEVER DEFINED -- assignment disabled

public:
  void contents() const; // FOR DEBUGGING ONLY

  const std::vector<PPMonoidElem>& myIndets() const;               ///< std::vector whose n-th entry is n-th indet as PPMonoidElem

  // The functions below are operations on power products owned by PPMonoidRingExpImpl
  const PPMonoidElem& myOne() const;
  using PPMonoidBase::myNew;    // disable warnings of overloading
  PPMonoidElemRawPtr myNew() const;                                ///< ctor from nothing
  PPMonoidElemRawPtr myNew(PPMonoidElemConstRawPtr rawpp) const;   ///< ctor by assuming ownership
  PPMonoidElemRawPtr myNew(const std::vector<long>& expv) const;   ///< ctor from exp vector
  void myDelete(RawPtr rawpp) const;                               ///< dtor, frees pp
  void mySwap(RawPtr rawpp1, RawPtr rawpp2) const;                 ///< swap(pp1, pp2);
  void myAssignOne(RawPtr rawpp) const;                            ///< pp = 1
  void myAssign(RawPtr rawpp, ConstRawPtr rawpp1) const;           ///< p = pp1
  void myAssign(RawPtr rawpp, const std::vector<long>& expv) const;///< pp = expv (assign from exp vector, all exps >= 0)

  void myMul(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const;   ///< pp = pp1*pp2
  using PPMonoidBase::myMulIndetPower;    // disable warnings of overloading
  void myMulIndetPower(RawPtr rawpp, long indet, long exp) const;           ///< pp *= indet^exp, assumes exp >= 0
  void myDiv(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const;   ///< pp = pp1/pp2
  void myColon(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const; ///< pp = pp1/gcd(pp1,pp2)
  void myGcd(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const;   ///< pp = gcd(pp1,pp2)
  void myLcm(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const;   ///< pp = lcm(pp1,pp2)
  void myRadical(RawPtr rawpp, ConstRawPtr rawpp1) const;                   ///< pp = radical(pp1)
  void myPowerSmallExp(RawPtr rawpp, ConstRawPtr rawpp1, long exp) const;   ///< pp = pp1^exp, assumes exp >= 0
  void myPowerRingElem(RawPtr rawpp, ConstRawPtr rawpp1, ConstRefRingElem exp) const;   ///< pp = pp1^exp, assumes exp >= 0

  bool myIsOne(ConstRawPtr rawpp) const;                                   ///< is pp = 1?
  bool myIsIndet(long& index, ConstRawPtr rawpp) const;                    ///< true iff pp is an indet
  bool myIsCoprime(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const;          ///< are pp1 & pp2 coprime?
  bool myIsEqual(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const;            ///< is pp1 equal to pp2?
  bool myIsDivisible(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const;        ///< does pp2 divide pp1?
  bool myIsRadical(ConstRawPtr rawpp) const;                               ///< is pp equal to its radical?

  int myCmp(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const;                 ///< -1,0,1 as pp1 < = > pp2
  long myStdDeg(ConstRawPtr rawpp) const;                                  ///< standard degree of pp
  void myWDeg(degree& d, ConstRawPtr rawpp) const;                         ///< d = grading(pp)
  int myCmpWDeg(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const;             ///< <0, =0, >0 as wdeg(pp1) < = > wdeg(pp2)
  int myCmpWDegPartial(ConstRawPtr rawpp1, ConstRawPtr rawpp2, long) const; ///< as myCmpWDeg wrt the first weights
  long myExponent(ConstRawPtr rawpp, long indet) const;                    ///< exponent of indet in pp
  void myBigExponent(BigInt& EXP, ConstRawPtr rawpp, long indet) const;    ///< EXP = exponent of indet in pp
  void myRingElemExponent(RingElem& EXP, ConstRawPtr rawpp, long i) const; ///< EXP = degree of ith indet in pp
  void myExponents(std::vector<long>& expv, ConstRawPtr rawpp) const;      ///< expv[i] = exponent(pp,i)
  void myBigExponents(std::vector<BigInt>& v, ConstRawPtr rawpp) const;    ///< get exponents, SHOULD BE myExponents ???
  void myOutputSelf(std::ostream& out) const;                              ///< out << PPM
  // INHERITED DEFINITION of virtual  void myOutput(std::ostream& out, ConstRawPtr rawpp) const;
  void myDebugPrint(std::ostream& out, ConstRawPtr rawpp) const;           ///< print pp in debugging format???


private: // auxiliary functions
  PPMonoidRingExpElem * myExpv(RawPtr) const;
  const PPMonoidRingExpElem * myExpv(ConstRawPtr) const;

  void myComputeDivMask(DivMask& dm, const DivMaskRule& DivMaskImpl, ConstRawPtr rawpp) const; ///< used by PPWithMask
  bool myCheckExponents(const std::vector<long>& expv) const;

private: // data members
  ring ExponentRing;
  vector<PPMonoidElem> myIndetVector; ///< the indets as PPMonoidElems
  auto_ptr<PPMonoidElem> myOnePtr;
};


// File local inline functions

inline PPMonoidRingExpElem * PPMonoidRingExpImpl::myExpv(RawPtr rawpp) const
{
  return static_cast<PPMonoidRingExpElem *>(rawpp.myRawPtr());
}


inline const PPMonoidRingExpElem * PPMonoidRingExpImpl::myExpv(ConstRawPtr rawpp) const
{
  return static_cast<const PPMonoidRingExpElem *>(rawpp.myRawPtr());
}


bool PPMonoidRingExpImpl::myCheckExponents(const std::vector<long>& expv) const
{
  // Check expv.size == myNumIndets.
  // Check exps are non-neg
  if (len(expv) != myNumIndets) return false;
  for (long i=0; i < myNumIndets; ++i)
    if (expv[i] < 0) return false;
  return true;
}


//----   Constructors & destructor   ----//

PPMonoidRingExpImpl::PPMonoidRingExpImpl(const std::vector<symbol>& IndetNames, const PPOrdering& ord, const ring& ExponentRing):
  PPMonoidBase(ord, IndetNames),
  ExponentRing(ExponentRing),
  myIndetVector()
{
  // XXX do something with PPOrdering ord

  myRefCountInc();  // this is needed for exception cleanliness, in case one of the lines below throws
  myOnePtr.reset(new PPMonoidElem(PPMonoid(this)));
  {
    // IMPORTANT: this block destroys pp *before* the call to myRefCountZero.
    PPMonoidElem pp(PPMonoid(this));
    vector<long> expv(myNumIndets);
    for (long i=0; i < myNumIndets; ++i)
      {
        expv[i] = 1;
        myAssign(raw(pp), expv);
        myIndetVector.push_back(pp);
        expv[i] = 0;
      }
  }
  myRefCountZero();
}


PPMonoidRingExpImpl::~PPMonoidRingExpImpl()
{}

/////////////////////////////////////////////////////////////////////////////



const std::vector<PPMonoidElem>& PPMonoidRingExpImpl::myIndets() const
{
  return myIndetVector;
}


const PPMonoidElem& PPMonoidRingExpImpl::myOne() const
{
  return *myOnePtr;
}


PPMonoidElemRawPtr PPMonoidRingExpImpl::myNew() const
{
  PPMonoidElemRawPtr rawpp(new PPMonoidRingExpElem(myNumIndets));
  myAssignOne(rawpp); // cannot throw
  return rawpp;
}

PPMonoidElemRawPtr PPMonoidRingExpImpl::myNew(PPMonoidElemConstRawPtr rawcopypp) const
{
  PPMonoidElemRawPtr rawpp(new PPMonoidRingExpElem(myNumIndets));
  myAssign(rawpp, rawcopypp); // cannot throw
  return rawpp;
}


PPMonoidElemRawPtr PPMonoidRingExpImpl::myNew(const std::vector<long>& expv) const
{
  CoCoA_ASSERT(myCheckExponents(expv));
  PPMonoidElemRawPtr rawpp(new PPMonoidRingExpElem(myNumIndets));
  myAssign(rawpp, expv); // cannot throw
  return rawpp;
}


void PPMonoidRingExpImpl::myAssignOne(RawPtr rawpp) const
{
  PPMonoidRingExpElem * const expv = myExpv(rawpp);
  for (long i = 0; i < myNumIndets; ++i)
    expv->exponents[i] = zero(ExponentRing);
}


void PPMonoidRingExpImpl::myAssign(RawPtr rawpp, ConstRawPtr rawpp1) const
{
  if (rawpp == rawpp1) return;

  myExpv(rawpp)->exponents = myExpv(rawpp1)->exponents;
}

void PPMonoidRingExpImpl::myAssign(RawPtr rawpp, const vector<long>& expv1) const
{
  CoCoA_ASSERT(myCheckExponents(expv1));

  PPMonoidRingExpElem * const expv = myExpv(rawpp);
  for (long i = 0; i < myNumIndets; ++i)
    expv->exponents[i] = RingElem(ExponentRing, expv1[i]);
}


void PPMonoidRingExpImpl::myDelete(RawPtr rawpp) const
{
  delete myExpv(rawpp);
}


void PPMonoidRingExpImpl::mySwap(RawPtr rawpp1, RawPtr rawpp2) const
{
  if (rawpp1 == rawpp2) return;

  PPMonoidRingExpElem * const expv1 = myExpv(rawpp1);
  PPMonoidRingExpElem * const expv2 = myExpv(rawpp2);

  expv1->exponents.swap(expv2->exponents);
}


void PPMonoidRingExpImpl::myMul(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
{
  // No worries about aliasing.
  PPMonoidRingExpElem * const expv = myExpv(rawpp);
  const PPMonoidRingExpElem * const expv1 = myExpv(rawpp1);
  const PPMonoidRingExpElem * const expv2 = myExpv(rawpp2);
  for (long i=0; i < myNumIndets; ++i)
    {
      expv->exponents[i] = expv1->exponents[i] + expv2->exponents[i];
    }
}


void PPMonoidRingExpImpl::myMulIndetPower(RawPtr rawpp, long indet, long exp) const  // assumes exp >= 0
{
  CoCoA_ASSERT(exp >= 0);
  CoCoA_ASSERT(0 <= indet && indet < myNumIndets);
  PPMonoidRingExpElem * const expv = myExpv(rawpp);
  expv->exponents[indet] += exp;
}


void PPMonoidRingExpImpl::myDiv(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
{
  // No worries about aliasing.
  PPMonoidRingExpElem * const expv = myExpv(rawpp);
  const PPMonoidRingExpElem * const expv1 = myExpv(rawpp1);
  const PPMonoidRingExpElem * const expv2 = myExpv(rawpp2);

  for (long i=0; i < myNumIndets; ++i)
    {
      CoCoA_ASSERT("Exponent Underflow" && expv1->exponents[i] >= expv2->exponents[i]);
      expv->exponents[i] = expv1->exponents[i] - expv2->exponents[i];
    }
}


void PPMonoidRingExpImpl::myColon(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
{
  // No worries about aliasing.
  PPMonoidRingExpElem * const expv = myExpv(rawpp);
  const PPMonoidRingExpElem * const expv1 = myExpv(rawpp1);
  const PPMonoidRingExpElem * const expv2 = myExpv(rawpp2);

  for (long i=0; i < myNumIndets; ++i)
    if (expv1->exponents[i] > expv2->exponents[i])
      expv->exponents[i] = expv1->exponents[i] - expv2->exponents[i];
    else
      expv->exponents[i] = zero(ExponentRing);
}


void PPMonoidRingExpImpl::myGcd(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
{
  // No worries about aliasing.
  PPMonoidRingExpElem * const expv = myExpv(rawpp);
  const PPMonoidRingExpElem * const expv1 = myExpv(rawpp1);
  const PPMonoidRingExpElem * const expv2 = myExpv(rawpp2);

  for (long i = 0; i < myNumIndets; ++i)
    expv->exponents[i] = min(expv1->exponents[i], expv2->exponents[i]);
}


void PPMonoidRingExpImpl::myLcm(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
{
  // No worries about aliasing.
  PPMonoidRingExpElem * const expv = myExpv(rawpp);
  const PPMonoidRingExpElem * const expv1 = myExpv(rawpp1);
  const PPMonoidRingExpElem * const expv2 = myExpv(rawpp2);

  for (long i = 0; i < myNumIndets; ++i)
    expv->exponents[i] = max(expv1->exponents[i], expv2->exponents[i]);
}


void PPMonoidRingExpImpl::myRadical(RawPtr rawpp, ConstRawPtr rawpp1) const
{
  PPMonoidRingExpElem * const expv = myExpv(rawpp);
  const PPMonoidRingExpElem * const expv1 = myExpv(rawpp1);

  for (long i = 0; i < myNumIndets; ++i)
    expv->exponents[i] = (expv1->exponents[i] > 0) ? one(ExponentRing) : zero(ExponentRing);
}


void PPMonoidRingExpImpl::myPowerSmallExp(RawPtr rawpp, ConstRawPtr rawpp1, long LongExp) const  // assumes exp >= 0
{
  CoCoA_ASSERT(LongExp >= 0);

  PPMonoidRingExpElem * const expv = myExpv(rawpp);
  const PPMonoidRingExpElem * const expv1 = myExpv(rawpp1);

  for (long i = 0; i < myNumIndets; ++i)
    expv->exponents[i] = LongExp * expv1->exponents[i];
}


void PPMonoidRingExpImpl::myPowerRingElem(RawPtr rawpp, ConstRawPtr rawpp1, ConstRefRingElem exp) const  // assumes exp >= 0
{
  PPMonoidRingExpElem * const expv = myExpv(rawpp);
  const PPMonoidRingExpElem * const expv1 = myExpv(rawpp1);

  for (long i = 0; i < myNumIndets; ++i)
    expv->exponents[i] = exp * expv1->exponents[i];
}


bool PPMonoidRingExpImpl::myIsOne(ConstRawPtr rawpp) const
{
  const PPMonoidRingExpElem * const expv = myExpv(rawpp);

  for (long i = 0; i < myNumIndets; ++i)
    if (! IsZero(expv->exponents[i])) return false;

  return true;
}


bool PPMonoidRingExpImpl::myIsIndet(long& index, ConstRawPtr rawpp) const
{
  const PPMonoidRingExpElem * const expv = myExpv(rawpp);
  long j = myNumIndets;
  for (long i = 0; i < myNumIndets; ++i)
    {
      if (IsZero(expv->exponents[i])) continue;
      if (j != myNumIndets || ! IsOne(expv->exponents[i])) return false;
      j = i;
    }
  if (j == myNumIndets) return false;
  index = j;
  return true;
}


bool PPMonoidRingExpImpl::myIsCoprime(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
{
  const PPMonoidRingExpElem * const expv1 = myExpv(rawpp1);
  const PPMonoidRingExpElem * const expv2 = myExpv(rawpp2);

  for (long i = 0; i < myNumIndets; ++i)
    if (! IsZero(expv1->exponents[i]) && ! IsZero(expv2->exponents[i])) return false;

  return true;
}


bool PPMonoidRingExpImpl::myIsEqual(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
{
  const PPMonoidRingExpElem * const expv1 = myExpv(rawpp1);
  const PPMonoidRingExpElem * const expv2 = myExpv(rawpp2);

  for (long i = 0; i < myNumIndets; ++i)
    if (expv1->exponents[i] != expv2->exponents[i]) return false;
  return true;
}


bool PPMonoidRingExpImpl::myIsDivisible(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
{
  const PPMonoidRingExpElem * const expv1 = myExpv(rawpp1);
  const PPMonoidRingExpElem * const expv2 = myExpv(rawpp2);

  for (long i = 0; i < myNumIndets; ++i)
    if (expv1->exponents[i] < expv2->exponents[i]) return false;

  return true;
}


bool PPMonoidRingExpImpl::myIsRadical(ConstRawPtr rawpp) const
{
  const PPMonoidRingExpElem * const expv = myExpv(rawpp);

  for (long i = 0; i < myNumIndets; ++i)
    if (expv->exponents[i] > one(ExponentRing)) return false;

  return true;
}


/* XXX myCmp currently uses lex ordering with no other option */

int PPMonoidRingExpImpl::myCmp(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
{
  const PPMonoidRingExpElem * const expv1 = myExpv(rawpp1);
  const PPMonoidRingExpElem * const expv2 = myExpv(rawpp2);

  for (long i=0; i<myNumIndets; ++i) {
    if (expv1->exponents[i] != expv2->exponents[i]) {
      if (expv1->exponents[i] > expv2->exponents[i]) return 1; else return -1;
    }
  }
  return 0;
}


long PPMonoidRingExpImpl::myStdDeg(ConstRawPtr) const
{
  // Returning zero here makes this act like a ungraded ring
  return 0;
}


void PPMonoidRingExpImpl::myWDeg(degree& d, ConstRawPtr) const
{
  // This is an ungraded ring.
  CoCoA_ASSERT(GradingDim(d) == 0);
}


int PPMonoidRingExpImpl::myCmpWDeg(ConstRawPtr, ConstRawPtr) const
{
  CoCoA_ERROR(ERR::NYI, "myCmpWDeg in PPMonoidRingExp");
}


int PPMonoidRingExpImpl::myCmpWDegPartial(ConstRawPtr, ConstRawPtr, long) const
{
  CoCoA_ERROR(ERR::NYI, "myCmpWDegPartial comparison in PPMonoidRingExp");
}


long PPMonoidRingExpImpl::myExponent(ConstRawPtr rawpp, long indet) const
{
  const PPMonoidRingExpElem * const expv = myExpv(rawpp);
  BigInt N;
  long n;

  CoCoA_ASSERT(indet < myNumIndets);

  if (! IsInteger(N, expv->exponents[indet])) {
    CoCoA_ERROR(ERR::BadConvert, "Exponent extraction in PPMonoidRingExp");
  }

  if (! IsConvertible(n, N)) {
    CoCoA_ERROR(ERR::ExpTooBig, "Exponent extraction in PPMonoidRingExp");
  }

  return n;
}

void PPMonoidRingExpImpl::myBigExponent(BigInt& EXP, ConstRawPtr rawpp, long indet) const
{
  const PPMonoidRingExpElem * const expv = myExpv(rawpp);

  CoCoA_ASSERT(indet < myNumIndets);

  if (! IsInteger(EXP, expv->exponents[indet])) {
    CoCoA_ERROR(ERR::BadConvert, "Exponent extraction in PPMonoidRingExp");
  }
}


void PPMonoidRingExpImpl::myRingElemExponent(RingElem& EXP, ConstRawPtr rawpp, long indet) const
{
  const PPMonoidRingExpElem * const expv = myExpv(rawpp);

  CoCoA_ASSERT(indet < myNumIndets);

  EXP = expv->exponents[indet];
}


void PPMonoidRingExpImpl::myExponents(std::vector<long>& v, ConstRawPtr rawpp) const
{
  const PPMonoidRingExpElem * const expv = myExpv(rawpp);
  BigInt N;
  long n;

  CoCoA_ASSERT(len(v) == myNumIndets);

  // Run this loop twice so we don't modify the vector if there's an exception

  for (long i=0; i < myNumIndets; ++i) {
    if (! IsInteger(N, expv->exponents[i])) {
      CoCoA_ERROR(ERR::BadConvert, "Exponent extraction in PPMonoidRingExp");
    }

    if (! IsConvertible(n, N)) {
      CoCoA_ERROR(ERR::ExpTooBig, "Exponent extraction in PPMonoidRingExp");
    }
  }

  for (long i=0; i < myNumIndets; ++i) {
    if (! IsInteger(N, expv->exponents[i])) {
      CoCoA_ERROR(ERR::BadConvert, "Exponent extraction in PPMonoidRingExp");
    }

    if (! IsConvertible(n, N)) {
      CoCoA_ERROR(ERR::ExpTooBig, "Exponent extraction in PPMonoidRingExp");
    }

    v[i] = n;
  }
}


void PPMonoidRingExpImpl::myBigExponents(std::vector<BigInt>& expvector, ConstRawPtr rawpp) const
{
  const PPMonoidRingExpElem * const expv = myExpv(rawpp);
  BigInt N;

  CoCoA_ASSERT(len(expv->exponents) == myNumIndets);

  // Run this loop twice so we don't modify the vector if there's an exception

  for (long i=0; i < myNumIndets; ++i) {
    if (! IsInteger(N, expv->exponents[i])) {
      CoCoA_ERROR(ERR::BadConvert, "Exponent extraction in PPMonoidRingExp");
    }
  }

  for (long i=0; i < myNumIndets; ++i) {
    if (! IsInteger(expvector[i], expv->exponents[i])) {
      CoCoA_ERROR(ERR::BadConvert, "Exponent extraction in PPMonoidRingExp");
    }
  }
}


void PPMonoidRingExpImpl::myComputeDivMask(DivMask&, const DivMaskRule&, ConstRawPtr) const
{
  // DivMasks are used for "fast" division checks.  We can do nothing here,
  // which forces all division checks to be "slow", i.e. use myIsDivisible()
}


void PPMonoidRingExpImpl::myOutputSelf(std::ostream& out) const
{
  out << "PPMonoidRingExp(" << myNumIndets << ", " << myOrd <<")";
}


void PPMonoidRingExpImpl::myDebugPrint(std::ostream& out, ConstRawPtr rawpp) const
{
  out << "DEBUG PP: myNumIndets=" << myNumIndets << ", exps=[";
  for (long i=0; i < myNumIndets; ++i)
    out << myExponent(rawpp, i) << " ";
  out << "]" << std::endl;
}




PPMonoid NewPPMonoidRing(const std::vector<symbol>& IndetNames, const PPOrdering& ord, const ring& ExponentRing = RingZZ())
{
  // Sanity check on the indet names given.
  const long nvars = NumIndets(ord);

  if (len(IndetNames) != nvars)
    CoCoA_ERROR(ERR::BadNumIndets, "NewPPMonoidRing(IndetNames,ord)");
  if (!AreDistinct(IndetNames))
    CoCoA_ERROR(ERR::BadIndetNames, "NewPPMonoidRing(IndetNames,ord)");
  if (!AreArityConsistent(IndetNames))
    CoCoA_ERROR(ERR::BadIndetNames, "NewPPMonoidRing(IndetNames,ord)");

  return PPMonoid(new PPMonoidRingExpImpl(IndetNames, ord, ExponentRing));
}

PPMonoid NewPPMonoidRing(const std::vector<symbol>& IndetNames, const PPOrderingCtor& ord, const ring& ExponentRing = RingZZ())
{
  //return NewPPMonoidEv(IndetNames, ord.myCtor(len(IndetNames)));
  return NewPPMonoidRing(IndetNames, ord.myCtor(len(IndetNames)), ExponentRing);
}

symbol mksymbol(string name)
{
  return symbol(name);
}

PPMonoid NewPPMonoidRing(const std::vector<string>& IndetNames, const PPOrderingCtor& ord = lex, const ring& ExponentRing = RingZZ())
{
  std::vector<symbol> symbols;
  transform(IndetNames.begin(), IndetNames.end(), back_inserter(symbols), mksymbol);
  return NewPPMonoidRing(symbols, ord.myCtor(len(symbols)), ExponentRing);
}




/* Differential - differential operator on a ring
 *
 * NOT a homomorphism (Leibniz rule doesn't preserve multiplication)
 *
 * Constructor takes a ring and a vector of substitution homomorphisms
 * that tell us how the indeterminates map.  It's done that way to
 * make usage syntax easy; we don't actually use the homomorphisms as
 * homomorphisms.
 */

class Differential {

public:

  const ring R;

  vector<pair<RingElem, RingElem>> difmap;

private:

  void insert(RingElem indet, RingElem result)
  {
    if (IsFractionField(owner(indet)) && IsOne(den(indet))) {
      indet = num(indet);
    }

    if (IsIndet(indet)) {
      //difmap[indet] = result;
      difmap.push_back(pair<RingElem,RingElem>(indet,result));
    } else {
      CoCoA_ERROR(ERR::NotIndet, "difmap insert");
    }
  }

public:

  Differential(ring R) : R(R) { }

  Differential(ring R, std::vector<RingHom> ringhoms) : R(R) {
    vector<RingElem> targets;
    RingHom ItoT = IdentityHom(R);

    if (IsFractionField(R)) {
      targets = indets(BaseRing(R));
      ItoT = CanonicalHom(BaseRing(R), R);
      for (auto &v: targets) v = ItoT(v);
    } else {
      targets = indets(R);
    }

    /* Keep track of which indeterminates have been replaced with
     * targets, otherwise we might replace "x" with "y" and then a
     * later homomorphism maps the "y" to a "z".  This is also why we
     * can't just multiply the vector of homomorphisms together into a
     * single homomorphism.
     */

    vector<bool> changed(targets.size(), false);

    for (RingHom hom: ringhoms) {
      if ((domain(hom) != R) || (codomain(hom) != R)) {
	CoCoA_ERROR(ERR::MixedRings, "creating differential");
      }
      for (unsigned int i=0; i < targets.size(); i++) {
	if (! changed[i] && hom(targets[i]) != targets[i]) {
	  insert(targets[i], hom(targets[i]));
	  targets[i] = hom(targets[i]);
	  changed[i] = true;
	}
      }
    }
  }

  RingElem operator() (const RingElem elem) const
  {
    BigRat q;

    if (IsRational(q,elem)) return zero(R);
    else if (IsFractionField(owner(elem))) {
      RingElem n = num(elem);
      RingElem d = den(elem);
      RingHom RtoF = CanonicalHom(owner(n), R);
      return ((*this)(n)*RtoF(d) - RtoF(n)*(*this)(d))/RtoF(d*d);
    } else {
      RingElem result = zero(R);
      for (auto it=difmap.begin(); it!=difmap.end(); ++it) {
	long index;
	RingHom PtoF = CanonicalHom(owner(elem), R);
	if (IsIndet(index, it->first)) {
	  result += PtoF(deriv(elem, index)) * it->second;
	}
      }
      return result;
    }
  }

  //  const ring& myOwner const { return R; }
};

ring owner(const Differential& D) {
  return D.R;
}

/* OrderedPolyRingElem - a polynomial ring constructed with an ordered coefficient ring
 *
 * Ordering is based on the coefficient of the most significant
 * monomial, according to the PPMonoid ordering.
 */

class OrderedPolyRingBase : public RingDistrMPolyCleanImpl {

  using RingDistrMPolyCleanImpl::RingDistrMPolyCleanImpl;

  virtual bool IamOrderedDomain(void) const {
    return true;
  }

  virtual int mySign(ConstRawPtr rawx) const {
    for (auto it=myBeginIter(rawx); it != myEndIter(rawx); it++) {
      return sign(coeff(it));
    }
    return 0;
  }

  virtual int myCmp(ConstRawPtr rawx, ConstRawPtr rawy) const {
    RawPtr rawdiff = myNew();
    mySub(rawdiff, rawx, rawy);
    int result = mySign(rawdiff);
    myDelete(rawdiff);
    return result;
  }
#if 0
  virtual int myCmpAbs(ConstRawPtr rawx, ConstRawPtr rawy) const;                  ///< equiv to myCmp(abs(x),abs(y))
  virtual bool myFloor(BigInt& N, ConstRawPtr rawx) const;                         ///< true iff x is integer; put floor(x) in N};
#endif
};

SparsePolyRing NewOrderedPolyRing(const ring& CoeffRing, const std::vector<symbol>& IndetNames, const PPOrderingCtor& ord = lex)
{
  return SparsePolyRing(new OrderedPolyRingBase(CoeffRing, NewPPMonoidEv(IndetNames, ord)));
}

/* PowerPolyRingBase - a polynomial ring whose exponents are polynomials
 *
 * It's not a Noetherian ring, so the Buchberger algorithm can run
 * indefinitely.  I patch up the GCD calculation to get something that
 * works by replacing f^p with a new variable g, computing the GCD
 * using the standard algorithm, then substituting back.
 */

class PowerPolyRingBase : public RingDistrMPolyCleanImpl {

private:

  /* Track indets in our exponents.  Something like f^(p-1) would map
   * into g, with f as the lower indet, p as the upper indet, g as the
   * tmpring indet, and -1 as the constant term.
   *
   * When we're inserting elements, or mapping into the temporary
   * ring, we want to map by the lower indet and the upper indet
   * together.
   *
   * We want to count how many total inserts we're made to figure
   * how many new indets to add to our temporary ring.
   *
   * When we're mapping back from the temporary ring, we want to
   * map by the lower indet only, and find out all temporary
   * indets correspond to that lower indet.
   *
   * exponentmap[lower_indet_number][upper_indet_PPMonoid]
   */

  struct indet_in_exponent {
    long tmpring_indet_number;
    RingElem upper_indet;
    long constant_term;
  };

  typedef std::map<long, std::map<PPMonoidElem, indet_in_exponent>> ExponentMap;

  void myGcd_find_RingElem_exponent(ConstRawPtr raw, ExponentMap& exponentMap) const {

    for (SparsePolyIter it=myBeginIter(raw); !IsEnded(it); ++it) {
      for (long indet=0; indet < NumIndets(myPPM()); indet ++) {
	RingElem exp = RingElemExponent(PP(it), indet);
	BigInt N;
	if (! IsInteger(N, exp)) {
	  /* exp is a RingElem in PPM's exponent ring.  We expect it
	   * to be a linear polynomial.  One of the linear terms (the
	   * last one processed) is adjusted to ensure that negative
	   * constant terms can be properly handled.
	   */
	  long a=1,b=0;
	  PPMonoidElem ppme(SparsePolyRingPtr(owner(exp))->myPPM());  // throws error if not a SparsePolyRing
	  for (auto expit = BeginIter(exp); !IsEnded(expit); ++expit) {
	    long i;
	    if (IsOne(PP(expit))) {
	      /* constant term - throws conversion exception if it isn't a long */
	      b = long(coeff(expit));
	    } else if (!IsIndet(i, PP(expit))) {
	      CoCoA_ERROR(ERR::NYI, "exponent is not a linear polynomial in PowerPolyRing GCD");
	    } else {
	      /* linear term - throws conversion exception if it isn't a long */
	      a = long(coeff(expit));
	      ppme = PP(expit);
	      exponentMap[indet][ppme].upper_indet = monomial(owner(exp), 1, ppme);
	    }
	  }

	  /* we track the minimum constant term, so if both f^p and f^(p-1)
	   * appear, we'll use g=f^(p-1) and write f^p as fg.  What about
	   * f^(2p-2)?  Then we want to use g=f^(p-1).  f^(2p-3) requires
	   * g=f^(p-2).  f^(2p+3) requires f^(p+1).
	   */

	  if (a > 1) {
	    if ((b < 0) && (b/a)*a != b) {
	      b = (b/a) - 1;
	    } else {
	      b = (b/a);
	    }
	  }

	  if (b < exponentMap[indet][ppme].constant_term) {
	    exponentMap[indet][ppme].constant_term = b;
	  }
	}
      }
    }
  }

  RingElem myGcd_rewrite_polynomial(ConstRawPtr raw, const SparsePolyRing& newRing, ExponentMap& exponentMap) const {

    RingElem newpoly(newRing);

    for (SparsePolyIter it=myBeginIter(raw); !IsEnded(it); ++it) {
      PPMonoidElem newPP(PPM(newRing));
      for (long indet=0; indet < NumIndets(myPPM()); indet ++) {
	if (exponentMap.count(indet) == 0) {
	  newPP *= IndetPower(PPM(newRing), indet, exponent(PP(it), indet));
	} else {
	  RingElem exp = RingElemExponent(PP(it), indet);
	  long a=0, b=0;
	  for (auto expit = BeginIter(exp); !IsEnded(expit); ++expit) {
	    if (IsOne(PP(expit))) {
	      b += long(coeff(expit));
	    } else {
	      a = long(coeff(expit));
	      newPP *= IndetPower(PPM(newRing), exponentMap[indet][PP(expit)].tmpring_indet_number, a);
	      b -= a * exponentMap[indet][PP(expit)].constant_term;
	    }
	  }
	  // if everything was done right in this function and the last, b should non-negative now
	  newPP *= IndetPower(PPM(newRing), indet, b);
	}
      }
      newpoly += monomial(newRing, coeff(it), newPP);
    }

    return newpoly;
  }

public:

  using RingDistrMPolyCleanImpl::RingDistrMPolyCleanImpl;

  void myGcd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override {

    ExponentMap exponentMap;

    myGcd_find_RingElem_exponent(rawx, exponentMap);
    myGcd_find_RingElem_exponent(rawy, exponentMap);

    /* If we didn't find an exponent that has to be modified, use our
     * underlying GCD implementation.
     */

    if (exponentMap.size() == 0) {
      RingDistrMPolyCleanImpl::myGcd(rawlhs, rawx, rawy);
      return;
    }

    /* Assign indet numbers in the newring to the exponentMap */
    long tmpring_indet_number = NumIndets(myPPM());
    for (auto it=exponentMap.begin(); it != exponentMap.end(); it++) {
      for (auto it2=it->second.begin(); it2 != it->second.end(); it2++) {
	it2->second.tmpring_indet_number = tmpring_indet_number ++;
      }
    }

    /* Otherwise, construct a new ring with a newly appended
     * indeterminate, map our two arguments into the new ring, compute
     * their GCD, then map the result back.
     */

    const std::vector<symbol> IndetNames = NewSymbols(tmpring_indet_number);
    PPMonoid NewPPM = NewPPMonoidNested(myPPM(), IndetNames, 0, WDegPosTO);
    SparsePolyRing NewPR(NewPolyRing(myCoeffRing(), NewPPM));

    RingElem newx = myGcd_rewrite_polynomial(rawx, NewPR, exponentMap);
    RingElem newy = myGcd_rewrite_polynomial(rawy, NewPR, exponentMap);

    RingElem GCD = gcd(newx, newy);

    SparsePolyRing P(this);

    myAssignZero(rawlhs);

    for (SparsePolyIter it=BeginIter(GCD); !IsEnded(it); ++it) {
      PPMonoidElem newPP(myPPM());
      for (long i=0; i < NumIndets(myPPM()); i ++) {
	if (exponentMap.count(i) == 0) {
	  newPP *= IndetPower(myPPM(), i, exponent(PP(it), i));
	} else {
	  RingElem exp(owner(exponentMap[i].begin()->second.upper_indet));
	  exp += exponent(PP(it), i);
	  for (auto it2=exponentMap[i].begin(); it2 != exponentMap[i].end(); it2++) {
	    long pow = exponent(PP(it), it2->second.tmpring_indet_number);
	    exp += pow * (it2->second.upper_indet + it2->second.constant_term);
	  }
	  newPP *= power(indet(myPPM(), i), exp);
	}
      }
      myAdd(rawlhs, rawlhs, raw(monomial(P, coeff(it), newPP)));
    }

  }
};

SparsePolyRing NewPowerPolyRing(const ring& CoeffRing, const PPMonoid& PPM) {
  return SparsePolyRing(new PowerPolyRingBase(CoeffRing, PPM));
}


/* Smith Normal Form
 *
 * Algorithm to compute a matrix's Smith normal form ported from
 * A. Bigatti <bigatti@dima.unige.it> S. DeFrancisci's CoCoA 4.6
 * version.
 *
 * Current version only works on integer matrices.
 *
 * For input matrix M, return a SmithRecord L with L.M in SNF and
 *
 *   L.M = L.U * M * L.V
 */

class SmithRecord {
public:
  matrix U;
  matrix M;
  matrix V;

  const long NR;
  const long NC;

  SmithRecord(const ConstMatrixView& U, const ConstMatrixView& M, const ConstMatrixView& V)
    : U(NewDenseMat(U)), M(NewDenseMat(M)), V(NewDenseMat(V)),
      NR(NumRows(M)), NC(NumCols(M)) {}
};

// LMinInII find the smallest element in the submatrix with row and
// col indices greater than or equal to i, and modifies the matrix
// to exchange the minimum element into entry (i,i)

RingElem FirstNonZero(matrix M, int I)
{
  for (int IP = I; IP < NumRows(M); ++ IP) {
    for (int JP = I; JP < NumCols(M); ++ JP) {
      if (! IsZero(M(IP, JP))) return M(IP, JP);
    }
  }
  return RingElem(RingOf(M));
}

pair<int,int> PosMinAbs(matrix M, int I)
{
  RingElem MinM = abs(FirstNonZero(M, I))+1;
  if (IsOne(MinM)) return pair<int,int>(I,I);

  int MinI = -1;
  int MinJ = -1;
  for (int IP = I; IP < NumRows(M); ++ IP) {
    for (int JP = I; JP < NumCols(M); ++ JP) {
      if (! IsZero(M(IP,JP)) && abs(M(IP,JP))<MinM) {
	MinM = abs(M(IP,JP));
	MinI = IP;
	MinJ = JP;
      }
    }
  }
  return pair<int,int>(MinI, MinJ);
}

pair<int,int> PosMinDeg(matrix M, int I)
{
  long MinDeg = deg(FirstNonZero(M, I))+1;
  if (MinDeg == 1) return pair<int,int>(I,I);

  int MinI = -1;
  int MinJ = -1;
  for (int IP = I; IP < NumRows(M); ++ IP) {
    for (int JP = I; JP < NumCols(M); ++ JP) {
      if (! IsZero(M(IP,JP)) && deg(M(IP,JP))<MinDeg) {
	MinDeg = deg(M(IP,JP));
	MinI = IP;
	MinJ = JP;
      }
    }
  }
  return pair<int,int>(MinI, MinJ);
}

void LMinInII(SmithRecord& L, int I, bool Pol)
{
  pair<int,int> p;

  if (Pol) {
    p = PosMinDeg(L.M, I);
  } else {
    p = PosMinAbs(L.M, I);
  }

  int IL = p.first;
  int JL = p.second;

  if (IL != I) {
    SwapRows(L.M, I, IL);
    SwapRows(L.U, I, IL);
  }

  if (JL != I) {
    SwapCols(L.M, I, JL);
    SwapCols(L.V, I, JL);
  }
}


// PDiv and PMod currently only work for integer RingElems

RingElem PDiv(ConstRefRingElem x, ConstRefRingElem y)
{
  BigInt X, Y;

  if (! IsInteger(X, x)) {
    CoCoA_ERROR(ERR::BadConvert, "Non-integer matrix in SmithNormalForm");
  }
  if (! IsInteger(Y, y)) {
    CoCoA_ERROR(ERR::BadConvert, "Non-integer matrix in SmithNormalForm");
  }

  return RingElem(owner(x), X/Y);
}

RingElem PMod(ConstRefRingElem x, ConstRefRingElem y)
{
  BigInt X, Y;

  if (! IsInteger(X, x)) {
    CoCoA_ERROR(ERR::BadConvert, "Non-integer matrix in SmithNormalForm");
  }
  if (! IsInteger(Y, y)) {
    CoCoA_ERROR(ERR::BadConvert, "Non-integer matrix in SmithNormalForm");
  }

  return RingElem(owner(x), X%Y);
}

// RecDiag transforms the matrix into diagonal form

void RecDiag(SmithRecord& L, int I, bool Pol)
{
  // assumiamo che la matrice abbia NR <= NC
  if (I >= L.NR) return;

  LMinInII(L,I,Pol);

  if (IsZero(L.M(I,I))) return;

  // lavora sulla i-ma riga
  for (int J = I+1; J < L.NC; J++) {
    while (! IsZero(L.M(I,J))) {
      RingElem Q = PDiv(L.M(I,J),L.M(I,I));
      L.M->myAddColMul(J,I,-Q);
      L.V->myAddColMul(J,I,-Q);
      // cerca il nuovo minimo
      LMinInII(L,I,Pol);
    }
  }

  // lavora sulla i-ma colonna
  int IR = I+1;
  bool OK = true;
  while ((IR < L.NR) && OK) {
    if (! IsZero(L.M(IR,I))) {
      RingElem Q1 = PDiv(L.M(IR,I),L.M(I,I));
      L.M->myAddRowMul(IR,I,-Q1);
      L.U->myAddRowMul(IR,I,-Q1);
    }
    if (! IsZero(L.M(IR,I))) {
      OK = false;
    }
    IR ++;
  }

  // chiamata ricorsiva
  if (OK) {
    RecDiag(L, I+1, Pol);
  } else {
    RecDiag(L, I, Pol);
  }
}

// Matrix is in diagonal form; adjust the diagonal so it satisfies the
// SNF divisibility condition.

void LSuperDiag(SmithRecord& L, int I, bool Pol)
{
  RingElem X;

  if (I+1 >= L.NR) return;

  L.M->myAddRowMul(I,I+1,RingElem(RingOf(L.M),1));
  L.U->myAddRowMul(I,I+1,RingElem(RingOf(L.M),1));
  LMinInII(L,I,Pol);

  if (IsZero(L.M(I,I))) return;

  //rende monico M[I,I] o cambia segno a colonna I
  if (Pol) {
    //X = 1/LC(M(I,I));
  } else {
    X = RingElem(RingOf(L.M), sign(L.M(I,I)));
  }
  if (! IsOne(X)) {
    L.M->myColMul(I,X);
    L.V->myColMul(I,X);
  }

  // lavora sulla i-ma riga
  while (! IsZero(L.M(I,I+1))) {
    //rende monico M[I,I+1] o cambia segno a colonna I+1
    if (Pol) {
      //X = 1/LC(M(I,I+1));
    } else {
      X = RingElem(RingOf(L.M), sign(L.M(I,I+1)));
    }
    if (! IsOne(X)) {
      L.M->myAddColMul(I,I+1,X);
      L.V->myAddColMul(I,I+1,X);
    }

    RingElem Q = PDiv(L.M(I,I+1),L.M(I,I));
    L.M->myAddColMul(I+1, I, -Q);
    L.V->myAddColMul(I+1, I, -Q);
    LMinInII(L,I,Pol); //cerca il nuovo minimo

    //rende monico M[I,I] o cambia segno a colonna I
    if (Pol) {
      //X = 1/LC(M(I,I));
    } else {
      X = RingElem(RingOf(L.M), sign(L.M(I,I)));
    }
    if (! IsOne(X)) {
      L.M->myColMul(I,X);
      L.V->myColMul(I,X);
    }
  }

  // lavora sull'elemento diagonale utilizzando l'i-ma colonna
  int IR = I+1;
  bool OK = true;
  while ((IR<L.NR) && OK) {
    if (! IsZero(L.M(IR,I))) {
      //rende monico M[IR,I] oppure cambia segno a M[IR,I]
      if (Pol) {
	//X = 1/LC(M(IR,I));
      } else {
	X = RingElem(RingOf(L.M), sign(L.M(IR,I)));
      }
      if (! IsOne(X)) {
	L.M->myRowMul(IR, X);
	L.U->myRowMul(IR, X);
      }

      RingElem Q1 = PDiv(L.M(IR,I),L.M(I,I));
      L.M->myAddRowMul(IR, I, -Q1);
      L.U->myAddRowMul(IR, I, -Q1);
    }
    if (! IsZero(L.M(IR,I))) {
      OK = false;
    }
    if (Pol) {
      //X = 1/LC(M(IR,IR));
    } else {
      X = RingElem(RingOf(L.M), sign(L.M(IR,IR)));
    }
    if (! IsOne(X)) {
      L.M->myColMul(IR,X);
      L.V->myColMul(IR,X);
    }
    IR = IR+1;
  }

  // chiamata ricorsiva
  if (OK) {
    LSuperDiag(L,I+1,Pol);
    if (! IsZero(PMod(L.M(I+1,I+1),L.M(I,I)))) {
      LSuperDiag(L,I,Pol);
    }
  } else {
    LSuperDiag(L,I,Pol);
  }
}

SmithRecord SmithFactor(matrix A)
{
  const long NR = NumRows(A);
  const long NC = NumCols(A);
  const ConstMatrixView IdNR = IdentityMat(RingOf(A), NR);
  const ConstMatrixView IdNC = IdentityMat(RingOf(A), NC);

  // An ordering is required.  Polynomial rings order by degree,
  // integer rings order by comparison.

  //bool Pol = IsPolyRing(RingOf(A));

  // always use integer algorithm for right now, as we might
  // be in a polynomial ring with only integer elements

  bool Pol = false;

  if (NR <= NC) {
    SmithRecord L(IdNR, A, IdNC);
    RecDiag(L, 0, Pol);
    LSuperDiag(L, 0, Pol);
    return L;
  } else {
    SmithRecord LT(IdNC, transpose(A), IdNR);
    RecDiag(LT, 0, Pol);
    LSuperDiag(LT, 0, Pol);
    SmithRecord L(transpose(L.V), transpose(L.M), transpose(L.U));
    return L;
  }
}

// convention: a function containing a "new" should be named "New.."
template<size_t n> matrix NewMatrixFromC(ring R, int (&cmat)[n][n])
{
  matrix M(NewDenseMat(R,n,n));
  
  for (int i=0; i < n; ++i)
    for (int j=0; j < n; ++j)
      SetEntry(M, i, j, cmat[i][j]);
  return M;
}

void testSmithFactor(void)
{
  GlobalManager CoCoAFoundations;

  cout << boolalpha; // so that bools print out as true/false
  cout << TeX;

#if 0
  // invariant factors 2,6,12

  int C_matrix[3][3] = {{ 2, 4, 4},
                        {-6, 6, 12},
                        {10,-4,-16}};
#else

  // invariant factors 1,3,21,0

  int C_matrix[4][4] = {{-6, 111, -36, 6},
                        { 5,-672, 210, 74},
                        { 0,-255,  81, 24},
                        {-7, 255, -81,-10}};
#endif

  matrix M_Z(NewMatrixFromC(RingZZ(), C_matrix));

  cout << M_Z << endl;

  SmithRecord L = SmithFactor(M_Z);

  cout << L.M << endl;
  cout << L.U << endl;
  cout << L.V << endl;
  cout << L.U * M_Z * L.V << endl;
}

/* Given a element in a polynomial ring, determine if it has any
 * integer solutions.
 */

bool DiophantineSolvable(ConstRefRingElem x)
{
  const factorization<RingElem> FacInfo = factor(x);
  const vector<RingElem>& IrredFacs = FacInfo.myFactors();
  const int NumIrredFacs = len(IrredFacs);

  CoCoA_ASSERT(IsSparsePolyRing(owner(x)));

  // Factor the polynomial and consider each factor individually.

  for (int i = 0; i != NumIrredFacs; ++i)
  {
    //cout << IrredFacs[i] << endl;

    matrix M(NewDenseMat(CoeffRing(owner(x)), 1, NumIndets(owner(x))));
    matrix C(NewDenseMat(CoeffRing(owner(x)), 1, 1));


    // Is this a linear factor?

    // If so, it is in the form ax + by + cz + ... 

    for (auto it=BeginIter(IrredFacs[i]); !IsEnded(it); ++it) {
      long idx;

      if (IsOne(PP(it))) {
	// constant term
	SetEntry(C, 0, 0, -coeff(it));
      } else if (! IsIndet(idx, PP(it))) {
	CoCoA_ERROR(ERR::NYI, "DiophantineSolvable: non-linear factor");
      } else {
	// linear term of form ax
	SetEntry(M, 0, idx, coeff(it));
      }
    }

    SmithRecord L = SmithFactor(M);
    matrix D = L.U * C;

    if (IsDivisible(D(0,0), L.M(0,0))) {
      // this factor can be solved with integers, so the whole equation can be, too
      return true;
    }
  }

  // none of the factors could be solved with integers
  return false;
}

void testDiophantineSolvable(void)
{
  GlobalManager CoCoAFoundations;

  cout << boolalpha; // so that bools print out as true/false
  cout << TeX;

  ring R = NewOrderedPolyRing(RingZZ(), vector<symbol> {symbol("a"), symbol("b"), symbol("c"), symbol("m"), symbol("n")});
  RingElem a(R, "a");
  RingElem b(R, "b");
  RingElem c(R, "c");
  RingElem m(R, "m");
  RingElem n(R, "n");

  RingElem t1 = (m+2)*(m+1);
  RingElem t2 = 2*m + 2*a - 1;

  CoCoA_ASSERT(DiophantineSolvable(t1));
  CoCoA_ASSERT(! DiophantineSolvable(t2));

  //cout << DiophantineSolvable(t1) << endl;
  //cout << DiophantineSolvable(t2) << endl;
}

/* A Weyl ring promoted to an operator algebra
 *
 * A vector of differentials is provided.  All must operate on the
 * same ring.  The symbols provided to the Weyl algebra must also
 * exist in that ring.  We assume that CanonicalHom can lift
 * coefficients from the Weyl algebra's coefficient ring
 * into the target ring of the differentials.
 */

class WeylOperatorAlgebra : public RingWeylImpl
{
private:
  const std::vector<Differential> differentials;
  const std::vector<symbol> SymList;
public:
  WeylOperatorAlgebra(const ring& CoeffRing, const std::vector<symbol>& names, const std::vector<long>& ElimIndets, const std::vector<Differential>& differentials):
    RingWeylImpl(CoeffRing, names, ElimIndets),
    differentials(differentials),
    SymList(names)
  {
    const ring R = owner(differentials[0]);
    for (auto it=differentials.begin(); it != differentials.end(); ++it) {
      if (owner(*it) != R) {
	CoCoA_ERROR(ERR::MixedRings, "WeylOperatorAlgebra: differentials must map the same ring");
      }
    }
    auto Rsymbols = symbols(R);
    for (auto it=names.begin(); it != names.end(); ++it) {
      // if (!any_of(Rsymbols.begin(), Rsymbols.end(), [] a));
    }
  }

  RingElem myLeftMul(ConstRawPtr rawx, ConstRefRingElem y) const
  {
    RingElem ans(owner(y));
    long myNumTrueIndets = myNumIndets()/2;

    if (owner(y) != owner(differentials[0])) {
      CoCoA_ERROR(ERR::MixedRings, "WeylOperatorAlgebra used on wrong ring");
    }
    for (auto it=myBeginIter(rawx); !IsEnded(it); ++it) {
      RingElem term(y);
      // each term in our operator is in canonical form, so operators are applied first
      for (long idx=0; idx < myNumTrueIndets; ++idx) {
	const long Didx = idx + myNumTrueIndets;
	const long d = exponent(PP(it), Didx);

	for (long i=1; i <= d; ++i) {
	  term = (differentials[idx])(term);
	  if (IsZero(term)) break;
	}
      }
      // then multiplication by indets
      for (long idx=0; idx < myNumTrueIndets; ++idx) {
	const long d = exponent(PP(it), idx);

	for (long i=1; i <= d; ++i) {
	  term = term * RingElem(owner(y), SymList[idx]);
	}
      }
      // inject the Weyl algebra's coefficient ring into the target ring
      ans += CanonicalHom(myCoeffRing(),owner(y))(coeff(it))*term;
    }
    return ans;
  }

  // Attempt to factor an operator out of a RingElem.  Pass in a
  // RingElem in the target ring along with an indeterminate in the
  // target ring.  If successful, returns a RingElem in the Weyl
  // algebra that operates on the indeterminate to produce the
  // RingElem.  Otherwise, returns 0.
  //
  // It's all heuristics; I don't have a general algorithm.
  //
  // Typical use is to factor a polynomial using a Weyl algebra
  // acting on the polynomial's fraction field.

  RingElem factor(ConstRefRingElem x, ConstRefRingElem arg) const
  {
    long myNumTrueIndets = myNumIndets()/2;
    ring WA = SparsePolyRing(this);
    RingElem solution(WA);

    ring K = owner(differentials[0]);
    ring R = owner(x);

    RingHom RtoK = CanonicalHom(R, K);

    if (!(K == R) && !(IsFractionField(K) && (BaseRing(K) == R))) {
      CoCoA_ERROR(ERR::MixedRings, "WeylOperatorAlgebra factor used on wrong ring");
    }

    if (K == R) {
      CoCoA_ERROR(ERR::NYI, "WeylOperatorAlgebra factor used on fraction field");
    }

    if (owner(arg) != K) {
      CoCoA_ERROR(ERR::MixedRings, "WeylOperatorAlgebra factor used on wrong ring");
    }

    if (owner(arg) != K) {
      CoCoA_ERROR(ERR::MixedRings, "WeylOperatorAlgebra factor used on wrong ring");
    }

    if (!IsSparsePolyRing(R)) {
      CoCoA_ERROR(ERR::NYI, "WeylOperatorAlgebra factor used on non-SparsePolyRing");
    }

    // compute a monoid consisting of just those indets we can inject from
    // the WeylAlgebra into K/R, and a RingHom that pulls their images
    // from R back into WA

    PPMonoidElem injectable_indets(PPM(R));
    vector<RingElem> pullbacks(NumIndets(R), RingElem(WA));

    for (long idx=0; idx < myNumTrueIndets; ++idx) {
      long indet_num;
      injectable_indets *= LPP(RingElem(R, SymList[idx]));
      CoCoA_ASSERT(IsIndet(indet_num, RingElem(R, SymList[idx])));
      pullbacks[indet_num] = indet(WA, idx);
    }

    RingHom pullback = PolyAlgebraHom(R, WA, pullbacks);

    // now we assume that R is a SparsePolyRing and iterate on x's terms

    for (auto it=BeginIter(x); !IsEnded(it); ++it) {
      RingElem term = myOne();
      PPMonoidElem pp = PP(it);
      long idx;

      //cerr << "trying " << coeff(it) << " " << PP(it) <<endl;

      // try each differential at orders 0, 1, and 2
      for (idx=0; idx < myNumTrueIndets; ++idx) {
	long order;
	for (order=0; order <= 2; ++order) {

	  const long Didx = idx + myNumTrueIndets;

	  // compute the differential applied to the indet

	  RingElem target(arg);
	  for (long i=1; i <= order; ++i) {
	    target = (differentials[idx])(target);
	  }

	  if (IsZero(target)) continue;

	  // can we exactly divide the target out of this term?

	  RingElem rem = RtoK(monomial(R, coeff(it), PP(it))) / target;
	  if (IsOne(den(rem))) {
	    // if so, then is the remainder composed of only injectable indets?
	    rem = num(rem);
	    if (IsMonomial(rem)) {
	      if (lcm(radical(LPP(rem)), injectable_indets) == injectable_indets) {
		// yes?  then we've solved this term
		solution += pullback(rem) * IndetPower(WA, Didx, order);
		//cerr << "solution is " << solution << endl;
		break;
	      }
	    }
	  }
	}

	if (order <= 2) break;

#if 0
	if (IsFractionField(owner(x))) {
	  // if we're in a fraction field, see if division yields a polynomial
	  if (IsOne(den(monomial(owner(x), coeff(it), PP(it)) / target))) break;
	} else {
	  // if we're in a polynomial ring, see if division is exact
	  if (IsDivisible(monomial(owner(x), coeff(it), PP(it)), target)) break;
	}
#endif
      }

      if (idx == myNumTrueIndets) {
	//CoCoA_ERROR(ERR::NYI, "WeylOperatorAlgebra factor failed");
	solution = 0;
	return solution;
      }

    }

    return solution;
  }

  /* Given a differential operator, find a basis for all solutions in C[x,t]
   *
   * Returns a list of RingElems in the Weyl algebra's target ring.
   *
   * XXX this isn't a good enough return value, as the basis set is often infinite
   *
   * Returns a recursion relationship that must be satisfied by the solution polynomials - how?
   */

  long long_gcd(long a, long b) const {
    return b == 0 ? a : long_gcd(b, a % b);
  }

  // I want a class to hold the coordinates of polynomial coefficients.
  //
  // Ex: a x^3y^2 term in the solution polynomial (roughly)
  // corresponds to a coordinate (3,2), except that the coordinates
  // are actually relative to a base coordinate, (0,0) in this case
  //
  // vector<long> doesn't allow me to subtract coordinates
  // component-wise
  //
  // valarray<long> doesn't allow me to compare coordinates
  // non-component-wise, required for using them as a key in a
  // std::map

  class coordinate_t : public vector<long>
  {
  public:

    using vector::vector;

    coordinate_t const operator-(const coordinate_t &b) const
    {
      coordinate_t c(size());
      std::transform(this->begin(),this->end(),b.begin(),c.begin(),std::minus<long>());
      return c;
    }

    // convert coordinate to a comma separated string

    operator std::string() const
    {
      std::string str;
      for (auto it=begin(); it != end(); ++it) {
	if (it != begin()) str += ",";
	str += std::to_string(*it);
      }
      return str;
    }
  };

  RingElem poly_solve(ConstRefRingElem diffop) const
  {
    // diffop is a RingElem in the WeylOperatorAlgebra
    ring WA = owner(diffop);
    long myNumTrueIndets = NumIndets(WA)/2;
    ring CoeffRing = NewPolyRing(RingZZ(), myNumTrueIndets);

    // Start by transforming diffop into a recursion relationship on
    // the solution polynomial's coefficients.

    map<coordinate_t, RingElem> coefficients;

    for (auto it=BeginIter(diffop); !IsEnded(it); ++it) {
      coordinate_t coefficient_coordinate(myNumTrueIndets);
      BigInt N;

      if (! IsInteger(N, coeff(it))) {
	CoCoA_ERROR(ERR::BadConvert, "Non-integer coefficient in Weyl algebra poly_solve");
      }

      RingElem coefficient(CoeffRing, N);

      // x^a dx^b transforms to (m+b-a)_(b) c[m+b-a]
      // where _(b) is Pochhammer falling factorial

      for (long idx=0; idx < myNumTrueIndets; ++idx) {
	long exp = exponent(PP(it), idx);
	long Dexp = exponent(PP(it), idx + myNumTrueIndets);

	coefficient_coordinate[idx] = Dexp - exp;

	for (long m=0; m < Dexp; m++) {
	  coefficient *= indet(CoeffRing, idx) + (Dexp - exp - m);
	}
      }

      if (coefficients.count(coefficient_coordinate) == 0) {
	coefficients[coefficient_coordinate] = coefficient;
      } else {
	coefficients[coefficient_coordinate] = coefficients[coefficient_coordinate] + coefficient;
      }
    }

    if (coefficients.size() == 1) {
      CoCoA_ERROR(ERR::NYI, "zero dimensional operator in poly_solve");
    }

    // can't be a zero-dimensional operator; is it one-dimensional?

    auto it=coefficients.begin();
    coordinate_t starting_coordinate = it->first;

    // use the first two coordinates to compute a slope, and remove
    // its GCD so we can easily tell if other coordinates are on the
    // same line

    it ++;
    coordinate_t slope = it->first - starting_coordinate;
    long first_nonzero_coord_index;

    for (long idx=0; idx < myNumTrueIndets; ++idx) {
      if (slope[idx] != 0) {
	first_nonzero_coord_index = idx;
	break;
      }
    }

    long slope_gcd = labs(slope[first_nonzero_coord_index]);

    for (long idx=0; idx < myNumTrueIndets; ++idx) {
      if (slope[idx] != 0) {
	slope_gcd = long_gcd(slope_gcd, labs(slope[idx]));
      }
    }
    for (long idx=0; idx < myNumTrueIndets; ++idx) {
      slope[idx] = slope[idx] / slope_gcd;
    }

    // if the coefficients all lie on the same line, associate with
    // each one a multiple of the slope, so that
    //
    // coordinate[i] = starting_coordinate + slope * multiple[i]

    //coordinate_t multiples(coefficients.size());
    //multiples[0] = 0;
    //multiples[1] = slope_gcd;

    int lowest_multiple = 0;
    coordinate_t lowest_multiple_coordinate = starting_coordinate;
    int highest_multiple = slope_gcd;
    coordinate_t highest_multiple_coordinate = it->first;

    // check any remaining coefficients to see if they're on the same line,
    // and compute the multiple if it is
    //
    // all we need to remember are the highest and lowest multiples

    for (++it; it != coefficients.end(); ++it) {

      coordinate_t difference = it->first - starting_coordinate;

      long multiple = difference[first_nonzero_coord_index] / slope[first_nonzero_coord_index];
      for (long idx=0; idx < myNumTrueIndets; ++idx) {
	if (multiple * slope[idx] != difference[idx]) {
	  CoCoA_ERROR(ERR::NYI, "operator dimension greater than one in poly_solve");
	}
      }
      //multiples[i] = multiple;
      if (multiple < lowest_multiple) {
	lowest_multiple = multiple;
	lowest_multiple_coordinate = it->first;
      }
      if (multiple > highest_multiple) {
	highest_multiple = multiple;
	highest_multiple_coordinate = it->first;
      }
    }

    // Now we know for sure that we've got a one-dimensional operator

    // We want to change variables into a system that separates a
    // coordinate along the line with coordinate orthogonal to it,
    // but I'll skip this step as unnecessary.

    // The coefficients of the lowest and highest multiple give
    // Diophantine equations that must be solvable for the operator
    // to have a solution.

    if (! DiophantineSolvable(coefficients[lowest_multiple_coordinate])
	|| ! DiophantineSolvable(coefficients[highest_multiple_coordinate])) {
      return zero(RingQQ());
    }

    // now construct a ring with symbols for each coefficient, use
    // it to construct the recursion relation, and return it

    vector<symbol> coefficient_symbols;

    for (it=coefficients.begin(); it != coefficients.end(); ++it) {
      coefficient_symbols.push_back(symbol("c_{" + (std::string)(it->first) + "}"));
    }

    ring RecursionRing = NewPolyRing(CoeffRing, coefficient_symbols);
    RingElem recursion(RecursionRing);

    int idx=0;
    for (it=coefficients.begin(); it != coefficients.end(); ++it, ++idx) {
      recursion += CoeffEmbeddingHom(RecursionRing)(it->second) * indet(RecursionRing, idx);
    }

    return recursion;
  }

};

const WeylOperatorAlgebra* WeylOperatorAlgebraPtr(const ring& R)
{
  return dynamic_cast<const WeylOperatorAlgebra*>(R.myRawPtr());
}

SparsePolyRing NewWeylOperatorAlgebra(const ring& CoeffRing, const std::vector<symbol>& names, const std::vector<Differential>& differentials)
{
  std::vector<long> ElimIndets;   // empty set
  return SparsePolyRing(new WeylOperatorAlgebra(CoeffRing, WANAMES(names), ElimIndets, differentials));
}

/* Return the smallest exponent with which an indeterminate appears in a polynomial */

RingElem minExponent(RingElem in, RingElem indet)
{
  RingElem poly;
  long index;
  bool valid = false;
  RingElem result;

  if (IsFractionField(owner(indet))) {
    CoCoA_ASSERT(IsOne(den(indet)));
    CoCoA_ASSERT(IsIndet(index, num(indet)));
  } else {
    CoCoA_ASSERT(IsIndet(index, indet));
  }

  if (IsFractionField(owner(in))) {
    CoCoA_ASSERT(IsOne(den(in)));
    poly = num(in);
  } else {
    poly = in;
  }

  for (SparsePolyIter it=BeginIter(poly); !IsEnded(it); ++it) {
    RingElem myexp = RingElemExponent(PP(it), index);
    if (!valid || (myexp < result)) {
      result = myexp;
      valid = true;
    }
  }

  return result;
}

/* Returns the factor by which an indeterminate appears to its minimal power */

RingElem minCoeff(RingElem in, RingElem myindet)
{
  RingElem poly;
  long index;

  if (IsFractionField(owner(myindet))) {
    CoCoA_ASSERT(IsOne(den(myindet)));
    CoCoA_ASSERT(IsIndet(index, num(myindet)));
  } else {
    CoCoA_ASSERT(IsIndet(index, myindet));
  }

  if (IsFractionField(owner(in))) {
    CoCoA_ASSERT(IsOne(den(in)));
    poly = num(in);
  } else {
    poly = in;
  }

  RingElem minexp = minExponent(in, myindet);
  RingElem result(owner(poly));

  for (SparsePolyIter it=BeginIter(poly); !IsEnded(it); ++it) {
    RingElem myexp = RingElemExponent(PP(it), index);
    if (myexp == minexp) {
      result += monomial(owner(poly), coeff(it), PP(it)/power(indet(owner(PP(it)), index), minexp));
    }
  }

  return result;
}

void program()
{
  GlobalManager CoCoAFoundations;

  cout << boolalpha; // so that bools print out as true/false
  cout << TeX;

  ring ZZ = RingZZ();
  ring QQ = RingQQ();

  // ExponentRing - these are the indeterminates that can appear in powers

  ring ExponentRing = NewOrderedPolyRing(ZZ, vector<symbol> {symbol("p"), symbol("a"), symbol("i"), symbol("b"), symbol("c")});
  RingElem p(ExponentRing, "p");
  RingElem a(ExponentRing, "a");
  RingElem b(ExponentRing, "b");
  RingElem c(ExponentRing, "c");
  RingElem i(ExponentRing, "i");

  // We now create a K[Z[p]] ring whose coefficient and exponent rings are ExponentRing,
  // along with its fraction field.

  PPMonoid PPM = NewPPMonoidRing(vector<string> {"x", "t", "z", "r", "T", "T_t", "(t+1)",
	"f", "f_x", "f_{xx}", "f_t", "q", "q_x", "q_{xx}", "q_t",
	"n", "n_{x}", "n_{xx}", "n_{t}",
	"n_i", "n_{ix}", "n_{ixx}", "n_{it}",
	"n_r", "n_{rx}", "n_{rxx}", "n_{rt}",
	"N", "N_x", "N_{xx}", "N_t", "D", "D_x", "D_{xx}", "D_t"}, lex, ExponentRing);
  ring R = NewPowerPolyRing(ExponentRing, PPM);
  ring K = NewFractionField(R);

  // x,t are in our field of definition
  // z = exp(-x^2/(4(t+1)))
  // r = sqrt(t)

  RingElem x(K, "x");
  RingElem t(K, "t");
  RingElem z(K, "z");
  RingElem r(K, "r");

  RingElem tpo(K, "(t+1)");

  // T is a polynomial in C[t] (doesn't involve x or z)
  RingElem T(K, "T");
  RingElem Tt(K, "T_t");

  // N is the numerator
  RingElem N(K, "N");
  RingElem Nx(K, "N_x");
  RingElem Nxx(K, "N_{xx}");
  RingElem Nt(K, "N_t");

  // D is the denominator
  RingElem D(K, "D");
  RingElem Dx(K, "D_x");
  RingElem Dxx(K, "D_{xx}");
  RingElem Dt(K, "D_t");

  // f and q are factors of something.  Typically D=f^p q,
  // where f is irreducible and q is coprime to f.
  RingElem f(K, "f");
  RingElem fx(K, "f_x");
  RingElem fxx(K, "f_{xx}");
  RingElem ft(K, "f_t");

  RingElem q(K, "q");
  RingElem qx(K, "q_x");
  RingElem qxx(K, "q_{xx}");
  RingElem qt(K, "q_t");

  RingElem n(K, "n");
  RingElem n_x(K, "n_{x}");
  RingElem n_xx(K, "n_{xx}");
  RingElem n_t(K, "n_{t}");

  RingElem n_r(K, "n_r");
  RingElem n_rx(K, "n_{rx}");
  RingElem n_rxx(K, "n_{rxx}");
  RingElem n_rt(K, "n_{rt}");

  // n_i is one coefficient in a numerator sum
  RingElem n_i(K, "n_i");
  RingElem n_ix(K, "n_{ix}");
  RingElem n_ixx(K, "n_{ixx}");
  RingElem n_it(K, "n_{it}");

  // setup our differentials (acting on K)

  Differential dx(K, vector<RingHom> {x >> 1, t >> 0, z >> -x/(2*t)*z, r >> 0, tpo >> 0,
	N >> Nx, Nx >> Nxx, D >> Dx, Dx >> Dxx, T >> 0,
	f >> fx, fx >> fxx, q >> qx, qx >> qxx,
	n >> n_x, n_x >> n_xx,
	n_r >> n_rx, n_rx >> n_rxx,
	n_i >> n_ix, n_ix >> n_ixx});

  Differential dt(K, vector<RingHom> {x >> 0, t >> 1, z >> power(x,2)/(4*power(t,2))*z, r >> r/(2*t), tpo >> 1,
	N >> Nt, D >> Dt, T >> Tt, f >> ft, q >> qt, n >> n_t, n_r >> n_rt, n_i >> n_it});

  // Create a Weyl algebra, with ExponentRing as the coefficient ring.
  // I don't actually use operators with coefficients not in QQ, but WA.factor() currently won't work
  // unless the operator algebra and the target ring share the same coefficient ring.

  ring WA = NewWeylOperatorAlgebra(ExponentRing, vector<symbol> {symbol("x"), symbol("t")}, vector<Differential> {dx, dt});

  RingElem WA_x(WA, "x");
  RingElem WA_t(WA, "t");
  RingElem WA_dx(WA, "dx");
  RingElem WA_dt(WA, "dt");

  // our operator

  RingElem O = WA_dx*WA_dx - WA_dt;

  cout << WeylOperatorAlgebraPtr(WA)->poly_solve(O) << endl;
  cout << WeylOperatorAlgebraPtr(WA)->poly_solve(-2*WA_x*WA_dx + 2*WA_t*WA_dx*WA_dx -2*WA_t*WA_dt -1) << endl;

  RingElem e = N/D;

  //RingHom rh = N >> 1;
  //cout << rh << endl;
  //cout << rh(e) << endl;

  //cout << power(x,p)/x << endl;
  //cout << dx(power(x,p)*D) << endl;
  //cout << deriv(power(x,p)*D,x) << endl;
  //cout << deriv(power(N,p)*D,N) << endl;

  cout << num(dx(dx(e)) - dt(e)) << endl;
  cout << num(2*t*dx(dx(e)) - 2*t*dt(e) - e) << endl;
  cout << num(-2*x*dx(e) + 2*t*dx(dx(e)) -2*t*dt(e) -e) << endl;
  cout << num(-2*x*dx(e) + 2*t*dx(dx(e)) -2*t*dt(e) -2*e) << endl;

  cout << endl;

  cout << O << endl;
  cout << num(O*e) << endl;

  cout << endl;
  cout << "try an irreducible factor f^p (p >= 2) in denominator = f^p q" << endl;
  cout << endl;

  RingElem d = (power(f,p) * q);

  //cout << d << endl;
  //cout << dx(d) << endl;
  //cout << power(dx(d),2) << endl;
  //cout << dx(dx(d)) << endl;

  //cout << (2*power(dx(d),2) - d*dx(dx(d)))/power(f,2*p-2) << endl;

  //cout << (power(f,p)) /power(f,2*p) << endl;
  //cout << (power(f,p-1))/power(f,2*p) << endl;
  //cout << (power(f,p) + power(f,p-1))/power(f,2*p) << endl;
  //cout << (power(f,2*p) - 1)/(power(f,p)-1) << endl;

  // This next problem's solved by making ExponentRing Z[p] instead of Q[p]

  //cout << "GCD: " << gcd(num(32*t*t + 64*t + 32), num(8*q)) << endl;
  //cout << "GCD: " << content(num(32*t*t + 64*t + 32)) << content(num(8*q)) << endl;
  //cout << "GCD: " << gcd(RingElem(ExponentRing, 32), RingElem(ExponentRing, 64)) << endl;
  //cout << "GCD: " << gcd(RingElem(ExponentRing, 32), RingElem(ExponentRing, 64)) << endl;

  RingElem eq = num(O*(N/d));

  cout << eq << endl;
  //cout << factor(eq) << endl;
  //cout << minExponent(eq, f) << endl;
  cout << "minCoeff(eq,f) = " << minCoeff(eq, f) << endl;

  cout << endl;
  cout << "try a square-free factor f in the denominator = f q" << endl;
  cout << endl;

  d = f * q;
  eq = num(O*(N/d));

  cout << eq << endl;
  //cout << minExponent(eq, f) << endl;
  cout << "minCoeff(eq, f) = " <<minCoeff(eq, f) << endl;




  cout << endl;
  cout << "try denominator z^p f" << endl;
  cout << endl;

  d = power(z,p) * f;
  eq = num(O*(N/d));
  //cout << eq << endl;
  //eq = (t >> (tpo - 1)) (CanonicalHom(R,K)(eq));
  cout << eq << endl;
  cout << "minCoeff(eq, f) = " << minCoeff(eq, f) << endl;


  cout << endl;
  cout << "try denominator z^p t^a f" << endl;
  cout << endl;

  d = power(z,p) * power(t,a) * f;
  eq = num(O*(N/d));
  //cout << eq << endl;
  //eq = (t >> (tpo - 1)) (CanonicalHom(R,K)(eq));
  cout << eq << endl;
  cout << "minCoeff(eq, f) = " << minCoeff(eq, f) << endl;
  //cout << den(dx(dx(N/d)) - dt(N/d)) << endl;


  cout << endl;
  cout << "try denominator z^p t^a T where T_x = 0" << endl;
  cout << endl;

  d = power(z,p) * power(t, a) * T;
  eq = num(O*(N/d));
  //cout << eq << endl;
  //eq = (t >> (tpo - 1)) (CanonicalHom(R,K)(eq));
  cout << eq << endl;
  cout << "minCoeff(eq, T) = " << minCoeff(eq, T) << endl;
  //cout << "minCoeff(eq, t+1) = " << minCoeff(eq, tpo) << endl;
  //cout << den(dx(dx(N/d)) - dt(N/d)) << endl;

  cout << endl;
  cout << "try denominator z^p" << endl;
  cout << endl;

  d = power(z,p);
  eq = num(O*(N/d));
  //cout << eq << endl;
  //eq = (t >> (tpo - 1)) (CanonicalHom(R,K)(eq));
  cout << eq << endl;
  //cout << minExponent(eq, tpo) << endl;
  cout << "minCoeff(eq, t) = " << minCoeff(eq, t) << endl;
  //cout << "minCoeff(eq, T) = " << minCoeff(eq, T) << endl;
  //cout << den(dx(dx(N/d)) - dt(N/d)) << endl;

  cout << endl;
  cout << "try denominator z^a with numerator t^b N where N has no t factor" << endl;
  cout << endl;

  d = power(z,a);
  RingElem NN = power(t,b) * N;
  eq = num(O*(NN/d));
  //cout << eq << endl;
  //eq = (t >> (tpo - 1)) (CanonicalHom(R,K)(eq));
  cout << eq << endl;
  cout << "minCoeff(eq, t) = " << minCoeff(eq, t) << endl;
  //cout << "minCoeff(eq, T) = " << minCoeff(eq, T) << endl;
  //cout << den(dx(dx(N/d)) - dt(N/d)) << endl;

  cout << endl;
  cout << "try denominator z^a with numerator n t^b + n_r t^c r where n and n_r have no t factor" << endl;
  cout << endl;

  d = power(z,a);
  NN = power(t,b) * n + power(t,c) * n_r * r;
  eq = num(O*(NN/d));
  //cout << eq << endl;
  //eq = (t >> (tpo - 1)) (CanonicalHom(R,K)(eq));
  cout << eq << endl;
  cout << "minCoeff(eq, r) = " << minCoeff(eq, r) << endl;
  cout << "eq - minCoeff(eq, r) = " << (eq - minCoeff(eq, r))/num(r) << endl;
  cout << "minCoeff(eq, t) = " << minCoeff(eq, t) << endl;
  //cout << "minCoeff(eq, T) = " << minCoeff(eq, T) << endl;
  //cout << den(dx(dx(N/d)) - dt(N/d)) << endl;

  cout << endl;
  cout << "try numerator n_i z^i" << endl;
  cout << endl;

  NN = n_i * power(z,i);
  eq = num(O*NN);
  //eq = (t >> (tpo - 1)) (CanonicalHom(R,K)(eq));
  cout << eq << endl;
  cout << (eq=minCoeff(eq,z)) << endl;
  cout << "minCoeff(eq, t) = " << minCoeff(eq,t) << endl;

  cout << endl;
  cout << "try numerator n_i t^a z^i" << endl;
  cout << endl;

  NN = n_i * power(z,i) * power(t,a);
  eq = num(O*NN);
  //eq = (t >> (tpo - 1)) (CanonicalHom(R,K)(eq));
  cout << eq << endl;
  cout << (eq=minCoeff(eq,z)) << endl;
  cout << minCoeff(eq,t) << endl;

  cout << endl;
  cout << "assume i = 1; numerator is n_i t^a z" << endl;
  cout << endl;

  NN = n_i * power(z,1) * power(t,a);
  eq = num(O*NN);
  cout << eq << endl;
  cout << (eq=minCoeff(eq,z)) << endl;
  cout << minCoeff(eq,t) << endl;

  cout << endl;
  cout << "assume i = 1; numerator is n_i z" << endl;
  cout << endl;

  NN = n_i * power(z,1);
  eq = num(O*NN);
  cout << eq << endl;
  cout << (eq=minCoeff(eq,z)) << endl;
  RingElem O1 = WeylOperatorAlgebraPtr(WA)->factor(eq, n_i);
  cout << O1 << endl;

  cout << endl;
  cout << "numerator is n + n_r r;    n_xx = n_t (i = 0)" << endl;
  NN = n + n_r * r;
  eq = num(O*NN);
  cout << eq << endl;
  cout << "minCoeff(eq, r) = " << minCoeff(eq, r) << endl;
  cout << "eq - minCoeff(eq, r) = " << (eq - minCoeff(eq, r))/num(r) << endl;

  RingElem O0a = WeylOperatorAlgebraPtr(WA)->factor(minCoeff(eq, r), n);
  RingElem O0b = WeylOperatorAlgebraPtr(WA)->factor((eq - minCoeff(eq, r))/num(r), n_r);

  cout << "O0a = " << O0a << endl;
  cout << "O0b = " << O0b << endl;

  // eq = num(-2*dx(NN)*x - NN + 2*t*dx(dx(NN)) - 2*t*dt(NN));
  //O = -2*WA_x*WA_dx - 1 + 2*WA_t*WA_dx*WA_dx - 2*WA_t*WA_dt;

  cout << endl;
  cout << "numerator is n + n_r r;    -2n_x x - n + 2 t n_{xx} = 2 t n_t  (i = 1)" << endl;
  NN = n + n_r * r;
  eq = num(O1*NN);
  cout << eq << endl;
  cout << "minCoeff(eq, r) = " << minCoeff(eq, r) << endl;
  cout << "eq - minCoeff(eq, r) = " << (eq - minCoeff(eq, r))/num(r) << endl;

  RingElem O1a = WeylOperatorAlgebraPtr(WA)->factor(minCoeff(eq, r), n);
  RingElem O1b = WeylOperatorAlgebraPtr(WA)->factor((eq - minCoeff(eq, r))/num(r), n_r);

  cout << "O1a = " << O1a << endl;
  cout << "O1b = " << O1b << endl;

  cout << endl;
  cout << "operator 0b" << endl;
  cout << "try a square-free factor f in the denominator = f q" << endl;
  cout << endl;

  //eq = num(2*t*dx(dx(e)) - 2*t*dt(e) - e);
  //O = 2*WA_t*WA_dx*WA_dx - 2*WA_t*WA_dt - 1;

  d = f * q;
  e = N/d;
  eq = num(O0b*e);
  cout << eq << endl;
  cout << "minCoeff(eq, f) = " <<minCoeff(eq, f) << endl;

  cout << endl;
  cout << "operator 0b" << endl;
  cout << "try a denominator = f t^a q" << endl;
  cout << endl;

  d = f * power(t,a) * q;
  e = N/d;
  eq = num(O0b*e);
  cout << eq << endl;
  cout << "minCoeff(eq, f) = " <<minCoeff(eq, f) << endl;

  cout << endl;
  cout << "operator 0b" << endl;
  cout << "try a denominator = f^a t^b q" << endl;
  cout << endl;

  d = power(f,a) * power(t,b) * q;
  e = N/d;
  eq = num(O0b*e);
  cout << eq << endl;
  cout << "minCoeff(eq, f) = " <<minCoeff(eq, f) << endl;


  //eq = num(-2*x*dx(e) + 2*t*dx(dx(e)) -2*t*dt(e) -e);
  //O = -2*WA_x*WA_dx + 2*WA_t*WA_dx*WA_dx - 2*WA_t*WA_dt - 1;

  cout << endl;
  cout << "operator 1a" << endl;
  cout << "try a denominator = f t^b q" << endl;
  cout << endl;

  d = f * power(t,b) * q;
  e = N/d;
  eq = num(O1a*e);
  cout << eq << endl;
  cout << "minCoeff(eq, f) = " <<minCoeff(eq, f) << endl;

  cout << endl;
  cout << "operator 1a" << endl;
  cout << "try a denominator = f^a t^b q" << endl;
  cout << endl;

  d = power(f,a) * power(t,b) * q;
  e = N/d;
  eq = num(O1a*e);
  cout << eq << endl;
  cout << "minCoeff(eq, f) = " <<minCoeff(eq, f) << endl;

  //eq = num(-2*x*dx(e) + 2*t*dx(dx(e)) -2*t*dt(e) -2*e);
  //O = -2*WA_x*WA_dx + 2*WA_t*WA_dx*WA_dx - 2*WA_t*WA_dt - 2;


  cout << endl;
  cout << "operator 1b" << endl;
  cout << "try a denominator = f t^b q" << endl;
  cout << endl;

  d = f * power(t,b) * q;
  e = N/d;
  eq = num(O1b*e);
  cout << eq << endl;
  cout << "minCoeff(eq, f) = " <<minCoeff(eq, f) << endl;

  cout << endl;
  cout << "operator 1b" << endl;
  cout << "try a denominator = f^a t^b q" << endl;
  cout << endl;

  d = power(f,a) * power(t,b) * q;
  e = N/d;
  eq = num(O1b*e);
  cout << eq << endl;
  cout << "minCoeff(eq, f) = " <<minCoeff(eq, f) << endl;



}

//----------------------------------------------------------------------
// Use main() to handle any uncaught exceptions and warn the user about them.
int main()
{
  try
  {
    program();
    //testSmithFactor();
    //testDiophantineSolvable();
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
