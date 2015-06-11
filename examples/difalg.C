
#include "CoCoA/library.H"
// #include <unordered_map>
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

private:

  ring R;

  vector<pair<RingElem, RingElem>> difmap;

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

  RingElem operator() (RingElem elem)
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
};

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

  using RingDistrMPolyCleanImpl::RingDistrMPolyCleanImpl;

  void myGcd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override {
    long modified_indet;
    RingElem exponent_RI_indet;
    long exponent_indet;
    long min_constant_term;
    bool have_modified_indet = false;

    /* XXX get both x and y */
    for (SparsePolyIter it=myBeginIter(rawx); !IsEnded(it); ++it) {
      for (long indet=0; indet < NumIndets(myPPM()); indet ++) {
	RingElem exp = RingElemExponent(PP(it), indet);
	BigInt N;
	if (! IsInteger(N, exp)) {
	  /* exp is a RingElem in PPM's exponent ring.  We expect it
	   * to be a polynomial in the form at+b.  Begin by extracing
	   * a and b, and ensuring that t is the same indeterminate
	   * seen in previous exponents (if any).
	   */
	  long a=1,b=0;
	  for (auto expit = BeginIter(exp); !IsEnded(expit); ++expit) {
	    long i;
	    if (IsOne(PP(expit))) {
	      /* b term - throws conversion exception if it isn't a long */
	      b = long(coeff(expit));
	    } else if (!IsIndet(i, PP(expit))) {
	      CoCoA_ERROR(ERR::NYI, "exponent not of form at+b in PowerPolyRing GCD");
	    } else if (have_modified_indet && (exponent_indet != i)) {
	      CoCoA_ERROR(ERR::NYI, "multiple indeterminates in exponent in PowerPolyRing GCD");
	    } else {
	      /* at term - throws conversion exception if it isn't a long */
	      a = long(coeff(expit));
	      exponent_indet = i;
	      //exponent_RI_indet = monomial(owner(exp), coeff(expit), PP(expit));
	      exponent_RI_indet = monomial(owner(exp), 1, PP(expit));
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

	  if (have_modified_indet) {
	    CoCoA_ASSERT(modified_indet == indet);
	    if (b < min_constant_term) {
	      min_constant_term = b;
	    }
	  } else {
	    have_modified_indet = true;
	    modified_indet = indet;
	    min_constant_term = b;
	  }
	}
      }
    }

    /* If we didn't find an exponent that has to be modified, use our
     * underlying GCD implementation.
     */

    if (! have_modified_indet) {
      RingDistrMPolyCleanImpl::myGcd(rawlhs, rawx, rawy);
      return;
    }

    /* Otherwise, construct a new ring with a newly appended indeterminate,
     * map our two arguments into the new ring, and compute their GCD.
     */

    const std::vector<symbol> IndetNames = NewSymbols(NumIndets(myPPM()) + 1);
    PPMonoid NewPPM = NewPPMonoidNested(myPPM(), IndetNames, 0, WDegPosTO);
    SparsePolyRing NewPR(NewPolyRing(myCoeffRing(), NewPPM));
    RingElem newx(NewPR), newy(NewPR);

    for (SparsePolyIter it=myBeginIter(rawx); !IsEnded(it); ++it) {
      PPMonoidElem newPP(NewPPM);
      for (long indet=0; indet < NumIndets(myPPM()); indet ++) {
	if (indet != modified_indet) {
	  newPP *= IndetPower(NewPPM, indet, exponent(PP(it), indet));
	} else {
	  RingElem exp = RingElemExponent(PP(it), indet);
	  long a=0, b=0;
	  for (auto expit = BeginIter(exp); !IsEnded(expit); ++expit) {
	    if (IsOne(PP(expit))) {
	      b = long(coeff(expit));
	    } else {
	      a = long(coeff(expit));
	    }
	  }
	  newPP *= IndetPower(NewPPM, indet, b - a*min_constant_term);
	  newPP *= IndetPower(NewPPM, NumIndets(myPPM()), a);
	}
      }
      newx += monomial(NewPR, coeff(it), newPP);
    }

    for (SparsePolyIter it=myBeginIter(rawy); !IsEnded(it); ++it) {
      PPMonoidElem newPP(NewPPM);
      for (long indet=0; indet < NumIndets(myPPM()); indet ++) {
	if (indet != modified_indet) {
	  newPP *= IndetPower(NewPPM, indet, exponent(PP(it), indet));
	} else {
	  RingElem exp = RingElemExponent(PP(it), indet);
	  long a=0, b=0;
	  for (auto expit = BeginIter(exp); !IsEnded(expit); ++expit) {
	    if (IsOne(PP(expit))) {
	      b = long(coeff(expit));
	    } else {
	      a = long(coeff(expit));
	    }
	  }
	  newPP *= IndetPower(NewPPM, indet, b - a*min_constant_term);
	  newPP *= IndetPower(NewPPM, NumIndets(myPPM()), a);
	}
      }
      newy += monomial(NewPR, coeff(it), newPP);
    }

    RingElem GCD = gcd(newx, newy);
    SparsePolyRing P(this);

    myAssignZero(rawlhs);

    for (SparsePolyIter it=BeginIter(GCD); !IsEnded(it); ++it) {
      PPMonoidElem newPP(myPPM());
      for (long i=0; i < NumIndets(myPPM()); i ++) {
	if (i != modified_indet) {
	  newPP *= IndetPower(myPPM(), i, exponent(PP(it), i));
	} else {
	  //RingElem exp(*exponent_ring, exponent(PP(it), i));
	  //exp += monomial(*exponent_ring, exponent(PP(it), NumIndets(myPPM())), indet(*exponent_PPM, exponent_indet));
	  RingElem exp = exponent(PP(it), NumIndets(myPPM())) * exponent_RI_indet + exponent(PP(it), i) + min_constant_term * exponent(PP(it), NumIndets(myPPM()));
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
    }
  }

  return result;
}

void program()
{
  GlobalManager CoCoAFoundations;

  cout << boolalpha; // so that bools print out as true/false
  cout << TeX;

  ring QQ = RingQQ();

  ring ExponentRing = NewOrderedPolyRing(QQ, vector<symbol> {symbol("p")});
  RingElem p(ExponentRing, "p");

  PPMonoid PPM = NewPPMonoidRing(vector<string> {"x", "t", "z",
	"f", "f_x", "f_{xx}", "f_t", "q", "q_x", "q_{xx}", "q_t",
	"N", "N_x", "N_{xx}", "N_t", "D", "D_x", "D_{xx}", "D_t"}, lex, ExponentRing);
  //ring R = NewPolyRing(ExponentRing, PPM);
  ring R = NewPowerPolyRing(ExponentRing, PPM);
  ring K = NewFractionField(R);

  RingElem x(K, "x");
  RingElem t(K, "t");
  RingElem z(K, "z");

  RingElem N(K, "N");
  RingElem Nx(K, "N_x");
  RingElem Nxx(K, "N_{xx}");
  RingElem Nt(K, "N_t");

  RingElem D(K, "D");
  RingElem Dx(K, "D_x");
  RingElem Dxx(K, "D_{xx}");
  RingElem Dt(K, "D_t");

  RingElem f(K, "f");
  RingElem fx(K, "f_x");
  RingElem fxx(K, "f_{xx}");
  RingElem ft(K, "f_t");

  RingElem q(K, "q");
  RingElem qx(K, "q_x");
  RingElem qxx(K, "q_{xx}");
  RingElem qt(K, "q_t");

  Differential dx(K, vector<RingHom> {x >> 1, t >> 0,
	N >> Nx, Nx >> Nxx, D >> Dx, Dx >> Dxx,
	f >> fx, fx >> fxx, q >> qx, qx >> qxx});

  Differential dt(K, vector<RingHom> {x >> 0, t >> 1,
	N >> Nt, D >> Dt, f >> ft, q >> qt});

  RingElem e = N/D;

  //RingHom rh = N >> 1;
  //cout << rh << endl;
  //cout << rh(e) << endl;

  //cout << power(x,p)/x << endl;
  //cout << dx(power(x,p)*D) << endl;
  //cout << deriv(power(x,p)*D,x) << endl;
  //cout << deriv(power(N,p)*D,N) << endl;

  cout << num(dx(dx(e)) - dt(e)) << endl;

  RingElem d = (power(f,p) * q);

  cout << d << endl;
  cout << dx(d) << endl;
  cout << power(dx(d),2) << endl;
  cout << dx(dx(d)) << endl;

  cout << (2*power(dx(d),2) - d*dx(dx(d)))/power(f,2*p-2) << endl;

  cout << (power(f,p)) /power(f,2*p) << endl;
  cout << (power(f,p-1))/power(f,2*p) << endl;
  cout << (power(f,p) + power(f,p-1))/power(f,2*p) << endl;

  RingElem eq = num(dx(dx(N/d)) - dt(N/d));

  cout << eq << endl;
  cout << minExponent(eq, f) << endl;

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
