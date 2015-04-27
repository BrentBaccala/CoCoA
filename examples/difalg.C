
#include "CoCoA/library.H"
// #include <unordered_map>
using namespace CoCoA;
using namespace std;

/* PPMonoidRingExpImpl
 *
 * A PPMonoid whose exponents are elements in a ring that is also an ordered domain.
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
      CoCoA_ASSERT("Exponent Overflow" && expv1[i] <= std::numeric_limits<SmallExponent_t>::max()-expv2[i]);
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
      CoCoA_ASSERT("Exponent Underflow" && expv1[i] >= expv2[i]);
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
    if (! IsZero(expv1->exponents[i]) != 0 && ! IsZero(expv2->exponents[i])) return false;

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


long PPMonoidRingExpImpl::myStdDeg(ConstRawPtr rawpp) const
{
  CoCoA_ERROR(ERR::NYI, "PPMonoid comparison in PPMonoidRingExp");
}


void PPMonoidRingExpImpl::myWDeg(degree& d, ConstRawPtr rawpp) const
{
  CoCoA_ERROR(ERR::NYI, "PPMonoid comparison in PPMonoidRingExp");
}


int PPMonoidRingExpImpl::myCmpWDeg(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
{
  CoCoA_ERROR(ERR::NYI, "PPMonoid comparison in PPMonoidRingExp");
}


int PPMonoidRingExpImpl::myCmpWDegPartial(ConstRawPtr rawpp1, ConstRawPtr rawpp2, long i) const
{
  CoCoA_ERROR(ERR::NYI, "PPMonoid comparison in PPMonoidRingExp");
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

  CoCoA_ASSERT(len(expv) == myNumIndets);

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


void PPMonoidRingExpImpl::myComputeDivMask(DivMask& dm, const DivMaskRule& DivMaskImpl, ConstRawPtr rawpp) const
{
  CoCoA_ERROR(ERR::NYI, "computeDivMask in PPMonoidRingExp");
  /* DivMaskImpl->myAssignFromExpv(dm, myExpv(rawpp), myNumIndets); */
}


void PPMonoidRingExpImpl::myOutputSelf(std::ostream& out) const
{
  out << "PPMonoidEv(" << myNumIndets << ", " << myOrd <<")";
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




/* Differential - differential operator on a ring
 *
 * NOT a homomorphism (Leibniz rule doesn't preserve multiplication)
 *
 * constructor takes a ring and a vector of substitution homomorphisms
 * that tell us how the indeterminates map
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

/* MyRingElem - extended version of RingElem
 *
 * 1. constructor MyRingElem(ring, string) as shorthand for RingElem(ring, symbol("string"))
 *
 * 2. Substitution homomorphism created with "x >> y" syntax ("x" must be an indeterminate)
 *    (only polynomial rings and their fraction fields are supported)
 *
 */

class MyRingElem : public RingElem {
public:
  using RingElem::RingElem;

  MyRingElem(ring R, std::string str) : RingElem(R, symbol(str)) {}

  RingHom operator>>(RingElem target) {
    long index;
    RingElem indet;
    ring R;

    if (owner(*this) != owner(target)) {
      CoCoA_ERROR(ERR::MixedRings, "creating substitution homomorphism");
    }

    if (IsFractionField(owner(*this)) && IsOne(den(*this))) {
      indet = num(*this);
    } else {
      indet = *this;
    }

    R = owner(indet);
    if (! IsPolyRing(R)) {
      CoCoA_ERROR(ERR::NotPolyRing, "creating substitution homomorphism");
    }
    if (! IsIndet(index, indet)) {
      CoCoA_ERROR(ERR::NotIndet, "creating substitution homomorphism");
    }

    vector<RingElem> targets = indets(R);

    if (R != owner(target)) {
      RingHom RtoT = CanonicalHom(R, owner(target));
      for (auto &v: targets) v = RtoT(v);
    }

    targets[index] = target;

    if (IsFractionField(owner(*this))) {
      return InducedHom((FractionField)(owner(*this)), PolyAlgebraHom(R, owner(target), targets));
    } else {
      return PolyAlgebraHom(R, owner(target), targets);
    }
  }

  RingHom operator>>(long target) {
    return (*this) >> ZZEmbeddingHom(owner(*this))(target);
  }
};

void program()
{
  GlobalManager CoCoAFoundations;

  cout << boolalpha; // so that bools print out as true/false

  ring QQ = RingQQ();
  ring ExponentRing = NewPolyRing(QQ, vector<symbol> {symbol("p")});
  //ring R = NewPolyRing(QQ, vector<symbol> {symbol("x"), symbol("t"), symbol("z"), symbol("N"), symbol("Nx"), symbol("Nxx"), symbol("Nt"), symbol("D"), symbol("Dx"), symbol("Dxx"), symbol("Dt")});
  PPMonoid PPM = NewPPMonoidRing(vector<symbol> {symbol("x"), symbol("t"), symbol("z"), symbol("N"), symbol("Nx"), symbol("Nxx"), symbol("Nt"), symbol("D"), symbol("Dx"), symbol("Dxx"), symbol("Dt")}, lex, ExponentRing);
  ring R = NewPolyRing(QQ, PPM);
  ring K = NewFractionField(R);

  MyRingElem x(K, "x");
  MyRingElem t(K, "t");
  MyRingElem z(K, "z");
  MyRingElem N(K, "N");
  MyRingElem Nx(K, "Nx");
  MyRingElem Nxx(K, "Nxx");
  MyRingElem Nt(K, "Nt");
  MyRingElem D(K, "D");
  MyRingElem Dx(K, "Dx");
  MyRingElem Dxx(K, "Dxx");
  MyRingElem Dt(K, "Dt");

  Differential dx(K, vector<RingHom> {x >> 1, t >> 0, N >> Nx, Nx >> Nxx, D >> Dx, Dx >> Dxx});

  Differential dt(K, vector<RingHom> {x >> 0, t >> 1, N >> Nt, D >> Dt});

  MyRingElem e = N/D;

  RingHom rh = N >> 1;

  cout << rh << endl;

  cout << rh(e) << endl;

  cout << num(dx(dx(e)) - dt(e)) << endl;

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
