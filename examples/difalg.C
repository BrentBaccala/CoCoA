
#include "CoCoA/library.H"
#include <unordered_set>
#include <algorithm>
#include <functional>
#include <string>
#include <regex>

#define USE_BLAD 1

#ifdef USE_BLAD
#include <blad.h>
#endif

using namespace CoCoA;
using namespace std;

/* C++ constructs */

/* Class 'bimap' is like 'map', but lookup is bi-directional.  The two
 * types have to be incompatible to allow overload resolution to
 * determine which lookup direction to use.
 *
 * XXX could be more efficient
 */

template<typename A, typename B, class Compare = less<A> >
class bimap : public std::map<A, B, Compare> {
public:
  using std::map<A, B, Compare>::map;    // inherit the constructors

  using std::map<A, B, Compare>::at;     // inherit the original A->B at()
  const A& at (const B& value)
  {
    for (auto it = std::map<A,B,Compare>::begin(); it != std::map<A,B,Compare>::end(); it++) {
      if (it->second == value) return it->first;
    }
    throw std::out_of_range("bimap");
  }

  using std::map<A, B, Compare>::operator[];   // inherit the original A->B operator[]
  const A& operator[] (const B& value)
  {
    return at(value);
  }

  using typename std::map<A, B, Compare>::size_type;
  using std::map<A, B, Compare>::count;

  size_type count (const B& value)
  {
    size_type result = 0;
    for (auto it = std::map<A,B,Compare>::begin(); it != std::map<A,B,Compare>::end(); it++) {
      if (it->second == value) result ++;
    }
    return result;
  }
};



/* PPMonoidRingExpImpl
 *
 * A PPMonoid whose exponents are elements in a ring that is also an
 * ordered domain.  Only "positive" exponents are allowed (positive in
 * the sense that they compare greater than zero).
 *
 * We also allow new symbols to be inserted into the PPMonoid, to
 * support differential rings with an infinite number of (algebraic)
 * generators.
 *
 * The only monomial ordering used is lexicographic, based on a total
 * degree ranking of the derivative terms.  The derivatives are
 * attached as tails by PowerPolyRing, separated from the heads by
 * underscores.  Identical heads are grouped together, along with
 * tails of the same length.
 *
 * For example, if the constructor was passed f, t, and q as initial
 * symbols, default lexicographic ranking is f > t > q, and the
 * derivative ranking might be like:
 *
 *       q < q_t < t < f < f_t < f_{tt} < f_{qt}
 *
 * Based on this ranking, lex is then used to compare monomials.
 *
 * Other orderings are not used because our implementation of the HDT
 * (highest derivative term) function assumes that the HDT appears in
 * the leading monomial, which can only be guaranteed by lex ordering.
 */

class PPMonoidRingExpImpl;

/* PPMonoidRingExpElem : an element of a PPMonoidRingExp
 *
 * internally implemented as a vector<RingElem> in the ExponentRing
 *
 * Because of the ability to dynamically add new symbols, different
 * elements might have different length vectors, so we overload
 * operator[] in a const (i.e, non-assignment) context and return zero
 * if we attempt to access past the end of the vector.  So you still
 * have to resize() when you assign new elements, but accessing
 * elements looks like an infinite zero-padded vector.
 */

class PPMonoidRingExpElem : protected std::vector<RingElem> {
  friend PPMonoidRingExpImpl;

  PPMonoidRingExpElem(long numIndets) {
    resize(numIndets);
  }

  using std::vector<RingElem>::operator[];

  const RingElem operator[](size_type pos) const {
    if (pos >= size()) return RingElem(owner(front()));
    else return std::vector<RingElem>::operator[](pos);
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
  const std::vector<PPMonoidElem>& myIndets() const;               ///< std::vector whose n-th entry is n-th indet as PPMonoidElem

  // The functions below are operations on power products owned by PPMonoidRingExpImpl
  ConstRefPPMonoidElem myNewSymbolValue(const symbol& s, ConstRefPPMonoidElem) const override;
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
  PPMonoidRingExpElem & myExpv(RawPtr) const;
  const PPMonoidRingExpElem & myExpv(ConstRawPtr) const;

  void myComputeDivMask(DivMask& dm, const DivMaskRule& DivMaskImpl, ConstRawPtr rawpp) const; ///< used by PPWithMask
  bool myCheckExponents(const std::vector<long>& expv) const;

private: // data members
  ring ExponentRing;
  vector<PPMonoidElem> myIndetVector; ///< the indets as PPMonoidElems
  unique_ptr<PPMonoidElem> myOnePtr;
  vector<long> ranking;  ///< the indices of the indets, in ranking order, from least to greatest
};


// File local inline functions

inline PPMonoidRingExpElem & PPMonoidRingExpImpl::myExpv(RawPtr rawpp) const
{
  return * (static_cast<PPMonoidRingExpElem *>(rawpp.myRawPtr()));
}


inline const PPMonoidRingExpElem & PPMonoidRingExpImpl::myExpv(ConstRawPtr rawpp) const
{
  return * (static_cast<const PPMonoidRingExpElem *>(rawpp.myRawPtr()));
}


bool PPMonoidRingExpImpl::myCheckExponents(const std::vector<long>& expv) const
{
  // CoCoA_ERROR(ERR::NYI, "myCheckExponents in PPMonoidRingExp");
  // Check vector is no bigger than our current number of indets
  if ((long) expv.size() > myNumIndets) return false;
  // Check exps are non-neg
  for (unsigned long i=0; i < expv.size(); ++i)
    if (expv[i] < 0) return false;
  return true;
}


//----   Constructors & destructor   ----//

PPMonoidRingExpImpl::PPMonoidRingExpImpl(const std::vector<symbol>& IndetNames, const PPOrdering& ord, const ring& ExponentRing):
  PPMonoidBase(ord, IndetNames),
  ExponentRing(ExponentRing),
  myIndetVector(),
  ranking(myNumIndets)
{
  if (! IsLex(ord)) {
    CoCoA_ERROR(ERR::NYI, "non-lexicographic ordering in PPMonoidRingExp");
  }
  if (! IsOrderedDomain(ExponentRing)) {
    CoCoA_ERROR(ERR::BadRing, "unordered ring used for ExponentRing in PPMonoidRingExp");
  }

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

  // initialize the ranking to the order passed in
  for (long i = 0; i < myNumIndets; ++i) {
    ranking[i] = myNumIndets - i - 1;
  }
}


PPMonoidRingExpImpl::~PPMonoidRingExpImpl()
{}

/////////////////////////////////////////////////////////////////////////////



const std::vector<PPMonoidElem>& PPMonoidRingExpImpl::myIndets() const
{
  return myIndetVector;
}

void BWBprint(ConstRefPPMonoidElem x)
{
  std::cerr << x << endl;
}

void BWBprint(ConstRefRingElem x)
{
  std::cerr << x << endl;
}

void BWBprint(ring x)
{
  std::cerr << x << endl;
}

void BWBprint(std::vector<RingElem> x)
{
  std::cerr << x << endl;
}

/* myNewSymbolValue() - add a new symbol to the PPMonoid, inserting it
 * into the ranking order before PPMonoidElem 'next'.  If the given
 * PPMonoidElem is one, insert at the end of the ranking order, so the
 * new element will become the least significant element.
 *
 * Yuck!  I now want to modify a const object by adding a new symbol
 * to the PPMonoid.  The rest of the library is so heavily invested in
 * PPMonoid's being const that changing it elsewhere is too painful to
 * consider.  Instead we play games with const_cast, to create local
 * variables that override the instance variables and makes them
 * writable.  My greatest concern is that myNumIndets can change
 * underfoot.
 *
 * There's a problem here related to CoCoA's memory management, which
 * is done with reference counting.  myIndetVector, in particular,
 * contains PPMonoidElem's that themselves contain pointers to the
 * parent PPMonoid, which creates a circular chain.  The standard
 * PPMonoid constructor resets the reference count to zero after
 * constructing myIndetVector (which never changes).  Yet here we add
 * a new element, so we should probably decrement our reference count
 * with myRefCountDec().  Yuck.  How brittle reference counting is.
 */

ConstRefPPMonoidElem PPMonoidRingExpImpl::myNewSymbolValue(const symbol& s, ConstRefPPMonoidElem next) const
{
  auto non_const_this = const_cast<PPMonoidRingExpImpl *>(this);
  int after_i = -1;

  if (! IsOne(next) && ! IsIndet(next)) {
    CoCoA_ERROR(ERR::NotIndet, "PPMonoidRingExpImpl::myNewSymbolValue");
  }

  for (long i=0; i < myNumIndets; ++i) {
    if ( s == myIndetSymbols[ranking[i]] ) {
      return myIndetVector[ranking[i]];
    }
    if (next == myIndetVector[ranking[i]]) {
      after_i = i;
    }
  }

  non_const_this->myIndetSymbols.push_back(s);
  non_const_this->myNumIndets ++;

  PPMonoidElem pp(PPMonoid(this));
  myMulIndetPower(raw(pp), myNumIndets-1, 1);

  non_const_this->myIndetVector.push_back(pp);

  // insert the new indet into the ranking, based on how its symbol compares

  // this local variable overrides the instance variable and makes it writable
  std::vector<long> & ranking = non_const_this->ranking;

  if (after_i == -1) {
    ranking.push_back(myNumIndets-1);
  } else {
    ranking.insert(ranking.begin() + after_i, myNumIndets-1);
  }

  // Now update the PPOrdering to be a lex ordering on the new ranking

  matrix OrderMatrix = NewDenseMat(RingZZ(), myNumIndets, myNumIndets);

  for (long i=0; i < myNumIndets; ++i) {
    SetEntry(OrderMatrix, i, ranking[myNumIndets - i - 1], 1);
  }

  non_const_this->myOrd = NewMatrixOrdering(myNumIndets, 0, OrderMatrix);
  //myOrd = NewMatrixOrdering(myNumIndets, 0, OrderMatrix);

  return myIndetVector.back();
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
  PPMonoidRingExpElem & expv = myExpv(rawpp);
  for (unsigned long i = 0; i < expv.size(); ++i)
    expv[i] = zero(ExponentRing);
}


void PPMonoidRingExpImpl::myAssign(RawPtr rawpp, ConstRawPtr rawpp1) const
{
  if (rawpp == rawpp1) return;

  myExpv(rawpp) = myExpv(rawpp1);
}

void PPMonoidRingExpImpl::myAssign(RawPtr rawpp, const vector<long>& expv1) const
{
  CoCoA_ASSERT(myCheckExponents(expv1));

  PPMonoidRingExpElem & expv = myExpv(rawpp);
  expv.resize(expv1.size());
  for (unsigned long i = 0; i < expv1.size(); ++i)
    expv[i] = RingElem(ExponentRing, expv1[i]);
}


void PPMonoidRingExpImpl::myDelete(RawPtr rawpp) const
{
  delete &myExpv(rawpp);
}


void PPMonoidRingExpImpl::mySwap(RawPtr rawpp1, RawPtr rawpp2) const
{
  if (rawpp1 == rawpp2) return;

  PPMonoidRingExpElem & expv1 = myExpv(rawpp1);
  PPMonoidRingExpElem & expv2 = myExpv(rawpp2);

  expv1.swap(expv2);
}


void PPMonoidRingExpImpl::myMul(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
{
  // No worries about aliasing.
  PPMonoidRingExpElem & expv = myExpv(rawpp);
  const PPMonoidRingExpElem & expv1 = myExpv(rawpp1);
  const PPMonoidRingExpElem & expv2 = myExpv(rawpp2);

  const auto nexps = max(expv1.size(), expv2.size());
  expv.resize(nexps);

  for (unsigned long i=0; i < nexps; ++i)
    {
      expv[i] = expv1[i] + expv2[i];
    }
}


void PPMonoidRingExpImpl::myMulIndetPower(RawPtr rawpp, long indet, long exp) const  // assumes exp >= 0
{
  CoCoA_ASSERT(exp >= 0);
  CoCoA_ASSERT(0 <= indet && indet < myNumIndets);
  PPMonoidRingExpElem & expv = myExpv(rawpp);
  if ((long) expv.size() <= indet) expv.resize(indet+1);
  expv[indet] += exp;
}


void PPMonoidRingExpImpl::myDiv(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
{
  // No worries about aliasing.
  PPMonoidRingExpElem & expv = myExpv(rawpp);
  const PPMonoidRingExpElem & expv1 = myExpv(rawpp1);
  const PPMonoidRingExpElem & expv2 = myExpv(rawpp2);

  const auto nexps = max(expv1.size(), expv2.size());
  expv.resize(nexps);

  for (unsigned long i=0; i < nexps; ++i)
    {
      CoCoA_ASSERT("Exponent Underflow" && expv1[i] >= expv2[i]);
      expv[i] = expv1[i] - expv2[i];
    }
}


void PPMonoidRingExpImpl::myColon(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
{
  // No worries about aliasing.
  PPMonoidRingExpElem & expv = myExpv(rawpp);
  const PPMonoidRingExpElem & expv1 = myExpv(rawpp1);
  const PPMonoidRingExpElem & expv2 = myExpv(rawpp2);

  const auto nexps = max(expv1.size(), expv2.size());
  expv.resize(nexps);

  for (unsigned long i=0; i < nexps; ++i)
    if (expv1[i] > expv2[i])
      expv[i] = expv1[i] - expv2[i];
    else
      expv[i] = zero(ExponentRing);
}


void PPMonoidRingExpImpl::myGcd(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
{
  // No worries about aliasing.
  PPMonoidRingExpElem & expv = myExpv(rawpp);
  const PPMonoidRingExpElem & expv1 = myExpv(rawpp1);
  const PPMonoidRingExpElem & expv2 = myExpv(rawpp2);

  const auto nexps = max(expv1.size(), expv2.size());
  expv.resize(nexps);

  for (unsigned long i = 0; i < nexps; ++i)
    expv[i] = min(expv1[i], expv2[i]);
}


void PPMonoidRingExpImpl::myLcm(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
{
  // No worries about aliasing.
  PPMonoidRingExpElem & expv = myExpv(rawpp);
  const PPMonoidRingExpElem & expv1 = myExpv(rawpp1);
  const PPMonoidRingExpElem & expv2 = myExpv(rawpp2);

  const auto nexps = max(expv1.size(), expv2.size());
  expv.resize(nexps);

  for (unsigned long i = 0; i < nexps; ++i)
    expv[i] = max(expv1[i], expv2[i]);
}


void PPMonoidRingExpImpl::myRadical(RawPtr rawpp, ConstRawPtr rawpp1) const
{
  PPMonoidRingExpElem & expv = myExpv(rawpp);
  const PPMonoidRingExpElem & expv1 = myExpv(rawpp1);

  const auto nexps = expv1.size();
  expv.resize(nexps);

  for (unsigned long i = 0; i < nexps; ++i)
    expv[i] = (expv1[i] > 0) ? one(ExponentRing) : zero(ExponentRing);
}


void PPMonoidRingExpImpl::myPowerSmallExp(RawPtr rawpp, ConstRawPtr rawpp1, long LongExp) const  // assumes exp >= 0
{
  CoCoA_ASSERT(LongExp >= 0);

  PPMonoidRingExpElem & expv = myExpv(rawpp);
  const PPMonoidRingExpElem & expv1 = myExpv(rawpp1);

  const auto nexps = expv1.size();
  expv.resize(nexps);

  for (unsigned long i = 0; i < nexps; ++i)
    expv[i] = LongExp * expv1[i];
}


void PPMonoidRingExpImpl::myPowerRingElem(RawPtr rawpp, ConstRawPtr rawpp1, ConstRefRingElem exp) const  // assumes exp >= 0
{
  PPMonoidRingExpElem & expv = myExpv(rawpp);
  const PPMonoidRingExpElem & expv1 = myExpv(rawpp1);

  const auto nexps = expv1.size();
  expv.resize(nexps);

  for (unsigned long i = 0; i < nexps; ++i)
    expv[i] = exp * expv1[i];
}


bool PPMonoidRingExpImpl::myIsOne(ConstRawPtr rawpp) const
{
  const PPMonoidRingExpElem & expv = myExpv(rawpp);
  const auto nexps = expv.size();

  for (unsigned long i = 0; i < nexps; ++i)
    if (! IsZero(expv[i])) return false;

  return true;
}


bool PPMonoidRingExpImpl::myIsIndet(long& index, ConstRawPtr rawpp) const
{
  const PPMonoidRingExpElem & expv = myExpv(rawpp);
  const auto nexps = expv.size();

  unsigned long j = nexps;
  for (unsigned long i = 0; i < nexps; ++i)
    {
      if (IsZero(expv[i])) continue;
      if (j != nexps || ! IsOne(expv[i])) return false;
      j = i;
    }
  if (j == nexps) return false;
  index = j;
  return true;
}


bool PPMonoidRingExpImpl::myIsCoprime(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
{
  const PPMonoidRingExpElem & expv1 = myExpv(rawpp1);
  const PPMonoidRingExpElem & expv2 = myExpv(rawpp2);

  const auto nexps = min(expv1.size(), expv2.size());

  for (unsigned long i = 0; i < nexps; ++i)
    if (! IsZero(expv1[i]) && ! IsZero(expv2[i])) return false;

  return true;
}


bool PPMonoidRingExpImpl::myIsEqual(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
{
  const PPMonoidRingExpElem & expv1 = myExpv(rawpp1);
  const PPMonoidRingExpElem & expv2 = myExpv(rawpp2);

  const auto nexps = max(expv1.size(), expv2.size());

  for (unsigned long i = 0; i < nexps; ++i)
    if (expv1[i] != expv2[i]) return false;

  return true;
}


bool PPMonoidRingExpImpl::myIsDivisible(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
{
  const PPMonoidRingExpElem & expv1 = myExpv(rawpp1);
  const PPMonoidRingExpElem & expv2 = myExpv(rawpp2);

  const auto nexps = max(expv1.size(), expv2.size());

  for (unsigned long i = 0; i < nexps; ++i)
    if (expv1[i] < expv2[i]) return false;

  return true;
}


bool PPMonoidRingExpImpl::myIsRadical(ConstRawPtr rawpp) const
{
  const PPMonoidRingExpElem & expv = myExpv(rawpp);
  const auto nexps = expv.size();

  for (unsigned long i = 0; i < nexps; ++i)
    if (expv[i] > one(ExponentRing)) return false;

  return true;
}


int PPMonoidRingExpImpl::myCmp(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
{
  const PPMonoidRingExpElem & expv1 = myExpv(rawpp1);
  const PPMonoidRingExpElem & expv2 = myExpv(rawpp2);

  // lex comparison

  for (long i=0; i < myNumIndets; ++i) {
    long r = ranking[myNumIndets - i - 1];
    if (expv1[r] != expv2[r]) {
      if (expv1[r] > expv2[r]) return 1; else return -1;
    }
  }
  return 0;
}


long PPMonoidRingExpImpl::myStdDeg(ConstRawPtr rawpp) const
{
  // Returning zero here makes this act like a ungraded ring
  // return 0;

  // ...but there's a problem.  SparsePolyRing's CoeffVecWRT assumes
  // that StdDeg returns the power of an indeterminate.

  // const PPMonoidRingExpElem & expv = myExpv(rawpp);
  long degree = 0;

  for (long i=0; i < myNumIndets; ++i) {
    // this doesn't work because we're adding a RingElem to a long
    // degree += expv[i];

    // ...but this does (and might throw an exception)
    degree += myExponent(rawpp, i);
  }

  return degree;
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
  const PPMonoidRingExpElem & expv = myExpv(rawpp);
  BigInt N;
  long n;

  CoCoA_ASSERT(indet < myNumIndets);

  if (indet < (long) expv.size()) {

    if (! IsInteger(N, expv[indet])) {
      CoCoA_ERROR(ERR::BadConvert, "Exponent extraction in PPMonoidRingExp");
    }

    if (! IsConvertible(n, N)) {
      CoCoA_ERROR(ERR::ExpTooBig, "Exponent extraction in PPMonoidRingExp");
    }

    return n;

  } else {

    return 0;
  }
}

void PPMonoidRingExpImpl::myBigExponent(BigInt& EXP, ConstRawPtr rawpp, long indet) const
{
  const PPMonoidRingExpElem & expv = myExpv(rawpp);

  CoCoA_ASSERT(indet < myNumIndets);

  if (indet < (long) expv.size()) {
    if (! IsInteger(EXP, expv[indet])) {
      CoCoA_ERROR(ERR::BadConvert, "Exponent extraction in PPMonoidRingExp");
    }
  } else {
    EXP = 0;
  }
}


void PPMonoidRingExpImpl::myRingElemExponent(RingElem& EXP, ConstRawPtr rawpp, long indet) const
{
  const PPMonoidRingExpElem & expv = myExpv(rawpp);

  CoCoA_ASSERT(indet < myNumIndets);

  if (indet < (long) expv.size()) {
    EXP = expv[indet];
  } else {
    EXP = RingElem(ExponentRing);
  }
}


void PPMonoidRingExpImpl::myExponents(std::vector<long>& v, ConstRawPtr rawpp) const
{
  const PPMonoidRingExpElem & expv = myExpv(rawpp);
  BigInt N;
  long n;

  CoCoA_ASSERT(len(v) == myNumIndets);

  // Run this loop twice so we don't modify the vector if there's an exception

  for (long i=0; i < myNumIndets; ++i) {
    if (! IsInteger(N, expv[i])) {
      CoCoA_ERROR(ERR::BadConvert, "Exponent extraction in PPMonoidRingExp");
    }

    if (! IsConvertible(n, N)) {
      CoCoA_ERROR(ERR::ExpTooBig, "Exponent extraction in PPMonoidRingExp");
    }
  }

  for (long i=0; i < myNumIndets; ++i) {
    if (! IsInteger(N, expv[i])) {
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
  const PPMonoidRingExpElem & expv = myExpv(rawpp);
  BigInt N;

  CoCoA_ASSERT((long) expv.size() == myNumIndets);

  // Run this loop twice so we don't modify the vector if there's an exception

  for (long i=0; i < myNumIndets; ++i) {
    if (! IsInteger(N, expv[i])) {
      CoCoA_ERROR(ERR::BadConvert, "Exponent extraction in PPMonoidRingExp");
    }
  }

  for (long i=0; i < myNumIndets; ++i) {
    if (! IsInteger(expvector[i], expv[i])) {
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
 *
 * A Differential can be modified, because it's sometimes useful to
 * use some of the indeterminates in the differential to define other
 * indeterminates.  For example, if z = exp(n/d), then dx(z) = z
 * dx(n/d), which would have to expanded by hand if we couldn't define
 * dx's action on n and d before defining it's action on z.
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

  void update(RingElem source, ConstRefRingElem target)
  {
    RingHom ItoT = IdentityHom(R);

    if ((owner(source) != R) || (owner(target) != R)) {
      CoCoA_ERROR(ERR::MixedRings, "Differential update source or target in wrong ring");
    }

    if (IsFractionField(R)) {
      if (! IsOne(den(source))) {
	CoCoA_ERROR(ERR::NotIndet, "Differential update");
      }
      source = num(source);
    }
    if (! IsIndet(source)) {
      CoCoA_ERROR(ERR::NotIndet, "Differential update");
    }

    // if a mapping already exists for this indet, replace it
    for (auto it = difmap.begin(); it != difmap.end(); ++ it) {
      if (it->first == source) {
	it->second = target;
	return;
      }
    }

    // otherwise, insert a new mapping
    insert(source, target);
  }

  virtual RingElem operator() (const RingElem elem) const
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

class BasicDifferential : public Differential {
private:
  ConstRefRingElem d;
public:

  BasicDifferential(ConstRefRingElem d) : Differential(owner(d)), d(d) { }

  RingElem operator() (const RingElem elem) const override
  {
    return deriv(elem, d);
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

	if (! IsInteger(exp)) {
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

void testPowerPolyRing(void)
{
  ring ZZ = RingZZ();
  ring QQ = RingQQ();

  // ExponentRing - these are the indeterminates that can appear in powers

  ring ExponentRing = NewOrderedPolyRing(ZZ, vector<symbol> {symbol("p")});
  RingElem p(ExponentRing, "p");

  // We now create a K[Z[p]] ring whose coefficient and exponent rings are ExponentRing,
  // along with its fraction field.

  PPMonoid PPM = NewPPMonoidRing(vector<string> {"f", "t", "q"}, lex, ExponentRing);
  ring R = NewPowerPolyRing(ExponentRing, PPM);
  ring K = NewFractionField(R);

  RingElem f(K, "f");
  RingElem t(K, "t");
  RingElem q(K, "q");

  // check term ordering

  CoCoA_ASSERT(LPP(num(f)) > LPP(num(t)));
  CoCoA_ASSERT(LPP(num(f)) > LPP(num(q)));
  CoCoA_ASSERT(LPP(num(t)) > LPP(num(q)));

  // CoCoA forbids mixed ring operations. 'p' in ExponentRing is
  // different from 'p' in K, which we now create.

  RingElem pK = EmbeddingHom(K)(CoeffEmbeddingHom(R)(p));

  cout << gcd(num(power(f,p)),num(f)) << endl;
  cout << power(f,p)/f << endl;
  //cout << (power(f,p)) /power(f,2*p) << endl;
  //cout << (power(f,p-1))/power(f,2*p) << endl;
  //cout << (power(f,p) + power(f,p-1))/power(f,2*p) << endl;
  //cout << (power(f,2*p) - 1)/(power(f,p)-1) << endl;

  CoCoA_ASSERT(power(f,p)/f == power(f,p-1));
  CoCoA_ASSERT(power(f,p) / power(f,2*p) == 1 / power(f,p));
  CoCoA_ASSERT(power(f,p-1) / power(f,2*p) == 1 / power(f,p+1));
  CoCoA_ASSERT((power(f,p) + power(f,p-1)) / power(f,2*p) == (f+1) / power(f,p+1));
  CoCoA_ASSERT((power(f,2*p) - 1) / (power(f,p) - 1) == power(f,p) + 1);

  // This next problem's solved by making ExponentRing Z[p] instead of Q[p]

  //cout << "GCD: " << gcd(num(32*t*t + 64*t + 32), num(8*q)) << endl;
  //cout << "GCD: " << content(num(32*t*t + 64*t + 32)) << " " << content(num(8*q)) << endl;
  //cout << "GCD: " << gcd(RingElem(ExponentRing, 32), RingElem(ExponentRing, 64)) << endl;
  //cout << "GCD: " << gcd(RingElem(ExponentRing, 32), RingElem(ExponentRing, 64)) << endl;

  CoCoA_ASSERT(gcd(num(32*t*t + 64*t + 32), num(8*q)) == 8);
  CoCoA_ASSERT(gcd(RingElem(ExponentRing, 32), RingElem(ExponentRing, 64)) == 32);

  //cout << deriv(power(f,p),f) << endl;

  CoCoA_ASSERT(deriv(power(f,p), f) == pK*power(f,p-1));
}

/* PowerPolyDifferentialRing - a polynomial ring whose exponents are
 * polynomials, and whose differentiation operator generates an
 * infinite number of indeterminates.
 *
 * The symbols in the ring have the form head_{tail} (LaTeX jet
 * notation).
 *
 * We assume that everything in the coefficient ring is a constant
 * (i.e, differentiates to zero).
 */

static std::string symbol_head(const std::string & symbol)
{
  const std::size_t underscore = symbol.find("_");

  if (underscore == std::string::npos) {
    return symbol;
  } else {
    return symbol.substr(0, underscore);
  }
}

static std::string symbol_tail(const std::string & symbol)
{
  const std::size_t underscore = symbol.find("_");
  const std::size_t openbracket = symbol.find("{");
  const std::size_t closebracket = symbol.find("}");

  if (underscore == std::string::npos) {
    return "";
  } else if (openbracket == std::string::npos) {
    return symbol.substr(underscore+1, std::string::npos);
  } else {
    return symbol.substr(openbracket+1, closebracket-openbracket-1);
  }
}

/* Splits a differential indeterminate into its base and derivative parts.
 *
 * Ex:  f_xt splits into f and x*t
 */

static void split_differential_indet (const PPMonoidElem &e, PPMonoidElem &base, PPMonoidElem &deriv)
{
  long i;

  if (! IsIndet(i, e)) {
    CoCoA_ERROR(ERR::NotIndet, "split_differential_indet");
  }

  const std::string sym = head(IndetSymbol(owner(e), i));

  base = owner(e)->mySymbolValue(symbol(symbol_head(sym)));

  AssignOne(deriv);
  for (char &c: symbol_tail(sym)) {
    deriv *= owner(e)->mySymbolValue(symbol(std::string(1, c)));
  }
}

const symbol& Symbol(ConstRefPPMonoidElem indet)
{
  long idx;
  CoCoA_ASSERT(IsIndet(idx, indet));
  return IndetSymbol(owner(indet), idx);
}

// scoped enum to avoid conflicting with existing 'lex'
enum class DifferentialRanking { lex, grlexA };

class PowerPolyDifferentialRingBase : public PowerPolyRingBase {

private:

  /* Given two strings (from two symbols), return the symbol string
   * corresponding to the derivative of the first by the second.
   *
   * (f,t) -> f_t      (f_t,t) -> f_{tt}     (f_t,x) -> f_{tx}
   *
   * Differentiation with respect to different variables is
   * commutative, so subscript symbols are always sorted into
   * accending alphabetical order.
   */

  std::string append_symbol(const std::string head, const std::string tail) const
  {
    const std::size_t underscore = head.find("_");

    // XXX wrong if tail contains extended TeX like \lambda
    CoCoA_ASSERT(tail.length() == 1);

    if (underscore == std::string::npos) {
      return head + "_" + tail;
    }

    std::string result = head;

    if (result.find("{", underscore) == std::string::npos) {
      result.insert(underscore+1, "{");
      result.append("}");
    }

    // XXX assumes that any previous subscripts are already sorted
    for (auto i=result.find("{")+1; i < result.length(); i++) {
      if ((result[i] > tail[0]) || (result[i] == '}')) {
	result.insert(i, tail);
	break;
      }
    }

    return result;
  }

  /* To speed things up, we keep a map from pairs of indet indices to
   * their derivative.  I.e, if indet 0 is 'x' and indet 1 is 't',
   * then <0,1> will map to 'x_t'.  This mapping is theoretically
   * infinite, so we only compute the entries when we need them, then
   * stash them here for future reference.  It's 'mutable' because
   * 'myDeriv' is declared 'const' way up in PolyRing.H and I don't
   * want to change that just to speed things up.
   *
   * This map creates circular reference that cause the ring to never
   * be destroyed.
   */

  mutable map<pair<long, long>, RingElem> differentiation_map;

public:

  // XXX don't really use this right now; only passed to the blad library
  DifferentialRanking ranking;

  PowerPolyDifferentialRingBase(const ring& R, const PPMonoid& PPM, DifferentialRanking ranking)
    : PowerPolyRingBase(R, PPM), ranking(ranking)
  { }

  /* We're inheriting from SparsePolyRing, whose implementation of
   * mySymbolValue and mySymbols obtains the symbols directly from the
   * underlying PPM, so the PPM's myNewSymbolValue takes care of
   * inserting new symbols.  We do have to insert a new indet into
   * myIndetVector, however.  Also, I want myNewSymbolValue() to be
   * const, so I can call it from myDeriv().  Thus we use const_cast.
   */

  RingElem myNewSymbolValue(const symbol& s)
  {
    vector<symbol> syms = symbols(ring(this));
    syms.push_back(s);

    // take care not to insert a new symbol if it already exists
    if ( ! AreDistinct(syms) ) {
      return mySymbolValue(s);
    }

    // construct a list of indets, sorted into ranking order
    std::vector<PPMonoidElem> indets = myPPM()->myIndets();
    std::sort(indets.begin(), indets.end());

    // comparison logic: first lexicographic ordering on the base symbol
    // (the part before an underscore), this gives us an elimination
    // ordering on the base symbols, then total degree on the
    // derivatives, then lexiographic on the derivatives

    // C++ regex doesn't work prior to gcc 4.9

    // const std::regex re("([^_]+)(_\\{?(.+)\\}?)?");
    // std::smatch newsymbol;
    // CoCoA_ASSERT(std::regex_match (head(s), newsymbol, re));

    // symbol[0] is entire match
    // symbol[1] is head
    // symbol[2] is tail with _ and braces
    // symbol[3] is tail stripped

    const std::string newsymbol_head = symbol_head(head(s));
    const std::string newsymbol_tail = symbol_tail(head(s));

    bool found_matching_head = false;

    // If we fall through the loop below, 'after' will be the identity
    // element in the monoid.

    PPMonoidElem after(myPPM());

    long i;
    for (i=0; i < myNumIndets() - 1; ++i) {

      //std::smatch oldsymbol;
      //CoCoA_ASSERT(std::regex_match (head(myIndetSymbols[ranking[i]]), oldsymbol, re));

      const symbol& oldsymbol = Symbol(indets[i]);

      const std::string oldsymbol_head = symbol_head(head(oldsymbol));
      const std::string oldsymbol_tail = symbol_tail(head(oldsymbol));

      // [Ma91]'s total degree ordering

      if (oldsymbol_head == newsymbol_head) {
	found_matching_head = true;
	if (oldsymbol_tail.length() > newsymbol_tail.length()) {
	  after = indets[i];
	  break;
	}
	if (oldsymbol_tail.length() == newsymbol_tail.length()) {
	  if (oldsymbol_tail < newsymbol_tail) {
	    // If they're the same length, then oldsymbol_tail
	    // compares less when the first character that differs is
	    // less in oldsymbol_tail.  Tail letters are sorted into
	    // ascending order, so "xxyz" < "xyyz".  The highest
	    // ranking derivative that differs ("xx") in this case,
	    // dominates.

	    // XXX using built-in string comparision here assumes that
	    // alphabetical order is the same as ranking order!
	    after = indets[i];
	    break;
	  }
	}
      } else if (found_matching_head) {
	// We've already found the matching head, but oldsymbol_head !=
	// newsymbol_head, so that means we've just passed the matching
	// head and moved onto the next head.  Insert new symbol here.
	after = indets[i];
	break;
      }
    }

    // Everything compared less than the new symbol, so we fell through
    // here.  New symbol goes at the end.

    //RingElem r = monomial(ring(this), 1, myPPM().PPMonoidRingExpImpl::myNewSymbolValue(s, after));
    //RingElem r = monomial(ring(this), 1, const_cast<ConstRefPPMonoidElem *(const symbol&, ConstRefPPMonoidElem)>(&myPPM()->myNewSymbolValue)->(s, after));
    //RingElem r = monomial(ring(this), 1, const_cast<PPMonoid&>(myPPM())->myNewSymbolValue(s, after));
    RingElem r = monomial(ring(this), 1, myPPM()->myNewSymbolValue(s, after));

    const_cast<PowerPolyDifferentialRingBase *>(this)->myIndetVector.push_back(r);
    //myIndetVector.push_back(r);

    return r;
  }


  void myDeriv(RawPtr rawlhs, ConstRawPtr rawf, ConstRawPtr rawx) const override
  {
    auto non_const_this = const_cast<PowerPolyDifferentialRingBase *>(this);

    if (myIsOne(rawx)) { myAssign(rawlhs, rawf); return; }

    const SparsePolyRing P(this);
    long lower_indet;

    CoCoA_ASSERT(IsIndet(lower_indet, myLPP(rawx)));

    RingElem ans(P);

    for (SparsePolyIter itf=myBeginIter(rawf); !IsEnded(itf); ++itf)
    {
      for (long indetn = 0; indetn < myNumIndets(); indetn ++) {

	/* First we compute x^n -> n*x^(n-1) */

	const RingElem d = RingElemExponent(PP(itf), indetn);

	if (IsZero(d)) continue;

	/* Differentiating x^n -> n*x^(n-1) requires multiplying
	 * monomials by exponents.  CanonicalHom is used to inject
	 * exponents into the polynomial ring.
	 */

	const RingElem scale = CanonicalHom(owner(d), P)(d);

	RingElem m(scale * monomial(P, coeff(itf), PP(itf)/indet(myPPM(), indetn)));

	/* Next, if upper and lower indets differ, we multiply by a
	 * partial derivative.  If upper and lower indets are the
	 * same, we leave it alone, since dx/dx = 1.
	 */

	if (indetn != lower_indet) {
	  const auto key = make_pair(indetn, lower_indet);
	  if (differentiation_map.count(key) == 0) {
	    const symbol & indet_symbol = myPPM()->myIndetSymbol(indetn);
	    const symbol & lower_indet_symbol = myPPM()->myIndetSymbol(lower_indet);

	    differentiation_map[key] = non_const_this->myNewSymbolValue(symbol(append_symbol(head(indet_symbol), head(lower_indet_symbol))));
	  }
	  m *= differentiation_map[key];
	}

	ans += m;
      }
    }
    mySwap(raw(ans), rawlhs); // really an assignment
  }

};

SparsePolyRing NewPowerPolyDifferentialRing(const ring& CoeffRing, const PPMonoid& PPM, DifferentialRanking ranking = DifferentialRanking::lex) {
  return SparsePolyRing(new PowerPolyDifferentialRingBase(CoeffRing, PPM, ranking));
}

void testPowerPolyDifferentialRing(void)
{
  ring ZZ = RingZZ();
  ring QQ = RingQQ();

  // ExponentRing - these are the indeterminates that can appear in powers

  ring ExponentRing = NewOrderedPolyRing(ZZ, vector<symbol> {symbol("p")});
  RingElem p(ExponentRing, "p");

  // We now create a K[Z[p]] ring whose coefficient and exponent rings are ExponentRing,
  // along with its fraction field.

  PPMonoid PPM = NewPPMonoidRing(vector<string> {"f", "t", "q"}, lex, ExponentRing);
  ring R = NewPowerPolyDifferentialRing(ExponentRing, PPM);
  ring K = NewFractionField(R);

  RingElem f(K, "f");
  RingElem t(K, "t");
  RingElem q(K, "q");

  // CoCoA forbids mixed ring operations. 'p' in ExponentRing is
  // different from 'p' in K, which we now create.

  RingElem pK = EmbeddingHom(K)(CoeffEmbeddingHom(R)(p));

  // Funny syntax - the derivatives don't exist in the ring until we
  // call 'deriv', so we do that first and then check to see if
  // they're right, instead of the other way around.

  RingElem ft = deriv(f,t);
  RingElem ftt = deriv(deriv(f,t),t);
  RingElem fqt = deriv(deriv(f,q),t);

  CoCoA_ASSERT(ft == RingElem(K, "f_t"));
  CoCoA_ASSERT(ftt == RingElem(K, "f_{tt}"));
  CoCoA_ASSERT(fqt == RingElem(K, "f_{qt}"));
  CoCoA_ASSERT(deriv(deriv(f,t),q) == fqt);

  // Test the ranking on the derivative terms.

  CoCoA_ASSERT(LPP(num(f)) < LPP(num(ft)));
  CoCoA_ASSERT(LPP(num(f)) < LPP(num(ftt)));
  CoCoA_ASSERT(LPP(num(ft)) < LPP(num(ftt)));
  CoCoA_ASSERT(LPP(num(f)) < LPP(num(fqt)));
  CoCoA_ASSERT(LPP(num(ft)) < LPP(num(fqt)));

  //cout << deriv(power(f,p),f) << endl;
  //cout << deriv(power(f,p),t) << endl;
  //cout << deriv(deriv(power(f,p),t),t) << endl;
  //cout << deriv(power(f,p)*q,t) << endl;
  //cout << deriv(power(f,p)/q,t) << endl;
  //cout << deriv(f/q,t) << endl;

  RingElem qt = deriv(q,t);

  CoCoA_ASSERT(deriv(power(f,p),f) == pK*power(f,p-1));
  CoCoA_ASSERT(deriv(power(f,p),t) == pK*power(f,p-1)*ft);
  CoCoA_ASSERT(deriv(deriv(power(f,p),t),t) == pK*power(f,p-1)*ftt + pK*(pK-1)*power(f,p-2)*ft*ft);
  CoCoA_ASSERT(deriv(power(f,p)*q,t) == power(f,p)*qt + pK*power(f,p-1)*q*ft);
  CoCoA_ASSERT(deriv(power(f,p)/q,t) == -(power(f,p)*qt - pK*power(f,p-1)*q*ft)/power(q,2));
  CoCoA_ASSERT(deriv(f/q,t) == -(f*qt - q*ft)/power(q,2));

  // This next problem's solved by making ExponentRing Z[p] instead of Q[p]

  //cout << "GCD: " << gcd(num(32*t*t + 64*t + 32), num(8*q)) << endl;
  //cout << "GCD: " << content(num(32*t*t + 64*t + 32)) << " " << content(num(8*q)) << endl;
  //cout << "GCD: " << gcd(RingElem(ExponentRing, 32), RingElem(ExponentRing, 64)) << endl;
  //cout << "GCD: " << gcd(RingElem(ExponentRing, 32), RingElem(ExponentRing, 64)) << endl;

  CoCoA_ASSERT(gcd(num(32*t*t + 64*t + 32), num(8*q)) == 8);
  CoCoA_ASSERT(gcd(RingElem(ExponentRing, 32), RingElem(ExponentRing, 64)) == 32);

  //cout << deriv(power(f,p),f) << endl;

  CoCoA_ASSERT(deriv(power(f,p), f) == EmbeddingHom(K)(CoeffEmbeddingHom(R)(p))*power(f,p-1));
}

/* A RegularDifferentialIdeal represents a differential ideal
 * (generated by addition, multiplication, and differentiation).
 *
 * A differential ideal J is said to be regular if there exists a
 * regular system Ω such that J = [A] : H Ω.
 *
 * An algebraic ideal J is said to be regular if there exists a
 * regular system Ω such that J = (A) : H Ω
 *
 * A system of differential equations and inequations is said to be
 * regular if the set of its equations is autoreduced and coherent,
 * the initial and separant of each equation appear among the
 * inequations and if its other inequations are partially reduced with
 * respect to the equations.  [Bo95] Definition 2.
 *
 * Rosenfeld-Groebner represents the least radical differential ideal
 * containing Σ as a finite intersection of radical differential
 * ideals, each described by a regular differential system.
 */

/* I want to handle sets of RingElem's, but C++'s std::set required
 * ordered comparison, so instead I used std::unordered_set, which
 * requires a hash function!  This defines a nil hash function
 * for RingElem's, to allow std::unordered_set<RingElem>.
 */

namespace std {
  template <>
  struct hash<RingElem> {
    std::size_t operator()(const RingElem&) const {
      return 0;
    }
  };
}

// The goal of the Rosenfeld-Groebner algorithm is to characterise a
// radical differential ideal as an intersection of regular
// differential ideals, which are in turn characterised by regular
// systems.
//
// [Bo95] lemma 4 tells us how to test membership in a regular
// differential ideal: we partial reduce modulo the equations, then
// see if the result lies in the algebraic ideal generated by the
// equation's ideal quotient'ed by the inequation's multiplicative
// family (i.e, the smallest multiplicatively closed set generated by
// the inequations).
//
// N.B. This ISN'T the ideal quotient, both since we can't use ideals
// with inequations (a ≠ 0 and b ≠ 0 doesn't imply a+b ≠ 0), and also
// since the ideal quotient I:T is ∀ t ∈ T, tq ∈ I, while this
// quotient is ∃ t ∈ T, tq ∈ I, i.e, the universal quantifier is
// replaced with an existential quantifier.

class RegularSystem {
  friend std::ostream & operator<<(std::ostream &out, const RegularSystem);
public:
  //std::unordered_set<RingElem> equations;
  //std::unordered_set<RingElem> inequations;
  std::vector<RingElem> equations;
  std::vector<RingElem> inequations;

  // (equations) colon-quotient (inequations)
  ideal I;

  RegularSystem(std::vector<RingElem> eq, std::vector<RingElem> ineq)
    : equations(eq.begin(), eq.end()), inequations(ineq.begin(), ineq.end()),
      // I(ineq.size() == 0 ? ideal(owner(eq[0]), eq) : colon(ideal(owner(eq[0]), eq), ideal(owner(eq[0]), ineq)))
      I(RingZZ(), std::vector<RingElem> {})
  {
    // compute I
  }
};

std::ostream & operator<<(std::ostream &out, const RegularSystem omega)
{
  out << "RegularSystem(";

#if 0
  bool first_eq_printed = false;
  bool first_ineq_printed = false;

  for (auto eq: omega.equations) {
    if (first_eq_printed) {
      out << ", ";
    }
    out << eq;
    first_eq_printed = true;
  }
  out << " \\ ";
  for (auto ineq: omega.inequations) {
    if (first_ineq_printed) {
    out << ", ";
    }
    out << ineq;
    first_ineq_printed = true;
  }
  out << " ; ";
#endif

  //out << gens(omega.I);
  out << omega.equations;

  out << ")";
    
  return out;
}

class RegularDifferentialIdeal {

  friend std::ostream & operator<<(std::ostream &out, const RegularDifferentialIdeal);

protected:

  const ring R;
  std::vector<RingElem> gens;

public:

  RegularDifferentialIdeal(const ring R, std::vector<RingElem> gens) :
    R(R), gens(gens) { }

  RegularDifferentialIdeal(const RingElem& r1) : R(owner(r1)) {
    gens.push_back(r1);
  }

  template<typename... Args>
  RegularDifferentialIdeal(const RingElem& r1, Args... rest) : RegularDifferentialIdeal(rest...) {
    gens.push_back(r1);
  }
};

namespace diffalg {

  // HDT - highest derivative term
  // Hp - highest power of HDT
  //
  // When passed a element from the coefficient ring, will return (1, 0)

  std::pair<PPMonoidElem, int> HDT_Hp(ConstRefRingElem f1)
  {
    if (IsConstant(f1)) {
      return make_pair(PPMonoidElem(PPM(owner(f1))), 0);
    }

    // XXX Assumes we're using a monomial ordering that places a
    // polynomial's HDT in its LPP.

    const PPMonoidElem& lpp = LPP(f1);
    PPMonoidElem hdt(owner(lpp));
    long hp = 0;

    // I count on 'one' comparing less than any indet, so hdt, as
    // initialized, will always compare less than any actual indet.

    for (long i = 0; i < NumIndets(owner(lpp)); ++ i) {
      if (exponent(lpp, i) > 0) {
	if (hdt < indet(owner(lpp), i)) {
	  hdt = indet(owner(lpp), i);
	  hp = exponent(lpp, i);
	}
	if ((hdt == indet(owner(lpp), i)) && hp < exponent(lpp, i)) {
	  hp = exponent(lpp, i);
	}
      }
    }

    return make_pair(hdt, hp);
  }

  PPMonoidElem HDT(ConstRefRingElem f1)
  {
    return HDT_Hp(f1).first;
  }

  int Hp(ConstRefRingElem f1)
  {
    return HDT_Hp(f1).second;
  }

  // The rank of a polynomial is its highest derivative term, raised
  // to its highest power.

  PPMonoidElem rank(ConstRefRingElem f)
  {
    const auto p = HDT_Hp(f);
    return power(p.first, p.second);
  }

  // The rank of a set of polynomial is the list of their ranks,
  // sorted into ascending order.  Kolchin §I.10

  std::vector<PPMonoidElem> rank(std::vector<RingElem> v)
  {
    std::vector<PPMonoidElem> result;
    for (auto e: v) {
      result.push_back(rank(e));
    }
    std::sort(result.begin(), result.end());
    return result;
  }

  // The total rank of a polynomial is the least common multiple of
  // its constituent monomials.

  PPMonoidElem total_rank(ConstRefRingElem f)
  {
    PPMonoidElem result(PPM(owner(f)));

    for (auto it=BeginIter(f); !IsEnded(it); ++it) {
      result = lcm(result, PP(it));
    }

    return result;
  }

  PPMonoidElem total_rank(std::vector<RingElem> v)
  {
    PPMonoidElem result(PPM(owner(v[0])));
    for (auto e: v) {
      result = lcm(result, total_rank(e));
    }
    return result;
  }

  RingElem Hcoeff(ConstRefRingElem f)
  {
    return CoeffVecWRT(f, monomial(owner(f), 1, HDT(f))).back();
  }

  // Hu - highest unknown, the unknown function appearing in the
  // highest derivative term.  Ex: if HDT is u_xx, then Hu is u.

  // If the RingElem is a constant, then Hu() returns 1.

  PPMonoidElem Hu(ConstRefRingElem f1)
  {
    if (! IsConstant(f1)) {
      PPMonoidElem e = HDT(f1);
      PPMonoidElem base(owner(e));
      PPMonoidElem deriv(owner(e));
      split_differential_indet(e, base, deriv);
      return base;
    } else {
      return PPMonoidElem(PPM(owner(f1)));
    }
  }

  // Elizabeth Mansfield's notation - alpha is the derivative

  PPMonoidElem alpha(ConstRefRingElem f1)
  {
    if (! IsConstant(f1)) {
      PPMonoidElem e = HDT(f1);
      PPMonoidElem base(owner(e));
      PPMonoidElem deriv(owner(e));
      split_differential_indet(e, base, deriv);
      return deriv;
    } else {
      // otherwise, return 1 (identity in the monoid)
      return PPMonoidElem(PPM(owner(f1)));
    }
  }

  // Differentiate a RingElem with respect to *all* of the indets in a
  // PPMonoidElem.  For example, differentiating u w.r.t to x^2*y
  // produces u_xxy, but RingElem can be more complicated than that.

  RingElem multideriv(RingElem r, ConstRefPPMonoidElem p)
  {
    for (long i = 0; i < NumIndets(owner(p)); ++ i) {
      for (long j = exponent(p, i); j > 0; -- j) {
	r = deriv(r, monomial(owner(r), 1, indet(owner(p), i)));
      }
    }
    return r;
  }

  /* Compute the Buchberger S-polynomial of two RingElems. */

  RingElem Spoly(ConstRefRingElem f1, ConstRefRingElem f2)
  {
    ConstRefPPMonoidElem lpp1 = LPP(f1);
    ConstRefPPMonoidElem lpp2 = LPP(f2);

    ConstRefRingElem c1 = LC(f1);
    ConstRefRingElem c2 = LC(f2);

    const PPMonoidElem LCM_pp = lcm(lpp1, lpp2);

    RingElem LCM_c;
    if (IsField(owner(c1))) {
      LCM_c = c1 * c2;
    } else {
      LCM_c = lcm(c1, c2);
    }

    RingElem a1 = monomial(owner(f1), LCM_c/c1, LCM_pp / lpp1);
    RingElem a2 = monomial(owner(f2), LCM_c/c2, LCM_pp / lpp2);

    return a1*f1 - a2*f2;
  }

  /* Compute the D-polynomial of two RingElems */

  RingElem Dpoly(ConstRefRingElem f1, ConstRefRingElem f2)
  {
    // XXX Assume we're using a monomial ordering so that a
    // polynomial's highest ranking unknown will appear in its LPP.

    // We want to find the highest ranking unknown in each polynomial.
    // If it's the same unknown, we can construct differential
    // operators that elevate each to a common derivative.

    // We check to see if f1 is constant to avoid the case where both
    // f1 and f2 are constant and Hu() therefore returns 1 for both.

    if (! IsConstant(f1) && (Hu(f1) == Hu(f2))) {
      PPMonoidElem hdt1 = HDT(f1);
      PPMonoidElem hdt2 = HDT(f2);
      PPMonoidElem base1(owner(hdt1));
      PPMonoidElem base2(owner(hdt2));
      PPMonoidElem deriv1(owner(hdt1));
      PPMonoidElem deriv2(owner(hdt2));

      // If Hu(f1) == Hu(f2), then the bases should be the same.
      split_differential_indet(hdt1, base1, deriv1);
      split_differential_indet(hdt2, base2, deriv2);

      CoCoA_ASSERT(base1 == base2);

      // compute a1 and a2, the differentials that raise HDT(f1) and
      // HDT(f2) to the same derivative

      PPMonoidElem lcm12 = lcm(deriv1, deriv2);
      PPMonoidElem a1 = lcm12 / deriv1;
      PPMonoidElem a2 = lcm12 / deriv2;

      // compute the differentials and cancel the results

      return Spoly(multideriv(f1, a1), multideriv(f2, a2));
    } else {
      //return zero(R);
      return zero(owner(f1));
    }
  }

  std::vector<RingElem> Dpoly(const std::vector<RingElem> set)
  {
    std::vector<RingElem> result;

    for (unsigned int i=0; i < set.size(); i++) {
      for (unsigned int j=i+1; j < set.size(); j++) {
	ConstRefRingElem dp = Dpoly(set[i], set[j]);
	if (! IsZero(dp)) {
	  result.push_back(dp);
	}
      }
    }

    // no simplification
    return result;
  }

  // This is partial differential reduction.

  RingElem partial_rem(ConstRefRingElem f, ConstRefRingElem g, bool full_reduction = false)
  {
    const PPMonoidElem HDT_g = HDT(g);
    const long Hp_g = Hp(g);
    const PPMonoidElem rank_g = rank(g);
    long indet_HDT_g;

    // attempted reduction by a constant will fail this ASSERT

    CoCoA_ASSERT(IsIndet(indet_HDT_g, HDT_g));

    PPMonoidElem base(owner(HDT_g));
    PPMonoidElem deriv(owner(HDT_g));

    if (IsConstant(f)) {
      return f;
    }

    split_differential_indet(HDT_g, base, deriv);

    // Construct a PPMonoidElem that is the product of all higher
    // derivative terms of HDT(g).

    PPMonoidElem higher_terms(owner(HDT_g));

    for (long i = 0; i < NumIndets(owner(HDT_g)); ++i) {
      PPMonoidElem c = indet(owner(HDT_g), i);

      PPMonoidElem base_c(owner(c));
      PPMonoidElem deriv_c(owner(c));

      split_differential_indet(c, base_c, deriv_c);

      if ((base == base_c) && (deriv_c != deriv) && IsDivisible(deriv_c, deriv)) {
	higher_terms *= c;
      }
    }

    // cout << "higher_terms = " << higher_terms << endl;

    // Now run through all of f's terms, to see if any of them involve
    // higher_terms.  If so, pseudo-reduce at the first term we find.

    for (SparsePolyIter it=BeginIter(f); !IsEnded(it); ++it) {
      PPMonoidElem common_terms = gcd(PP(it), higher_terms);
      if (! IsOne(common_terms)) {
	for (long i = 0; i < NumIndets(owner(HDT_g)); ++i) {
	  if (exponent(common_terms, i) != 0) {
	    const long p = exponent(PP(it), i);
	    const PPMonoidElem DT_PP = indet(owner(HDT_g), i);
	    const RingElem DT = monomial(owner(f), 1, DT_PP);

	    RingElem Dag = multideriv(g, alpha(DT)/deriv);

	    // cout << "Partial Reduction Step" << endl;
	    // cout << "f = " << f << endl;
	    // cout << "g = " << g << endl;
	    // cout << "alpha/deriv = " << alpha(DT)/deriv << endl;
	    // cout << "Dag = " << Dag << endl;
	    // cout << "rank = DT^p = " << DT << "^" << p << endl;
	    // cout << "CoeffVecWRT(f, DT) = " << CoeffVecWRT(f, DT) << endl;
	    // cout << "result = " << Hcoeff(Dag) * f - CoeffVecWRT(f, DT)[p] * power(DT, p-1) * Dag << endl;

	    // If I use ConstRefRingElem's here, I lose my signs!  -1 in CoeffVecWRT(f, DT) becomes 1!

	    RingElem gcoeff = Hcoeff(Dag);
	    RingElem fcoeff = CoeffVecWRT(f, DT)[p];
	    RingElem gcd_coeff = gcd(fcoeff, gcoeff);

	    // cout << "gcoeff = " << gcoeff << endl;
	    // cout << "fcoeff = " << fcoeff << endl;
	    // cout << "gcd = " << gcd_coeff << endl;

	    // return Hcoeff(Dag) * f - CoeffVecWRT(f, DT)[p] * power(DT, p-1) * Dag;
	    return (gcoeff * f - fcoeff * power(DT, p-1) * Dag) / gcd_coeff;
	  }
	}
      }
      if (full_reduction) {
	const long p = exponent(PP(it), indet_HDT_g);
	if (p >= Hp_g) {
	  const RingElem DT = monomial(owner(f), 1, HDT_g);

	  // cout << "Full Reduction Step" << endl;
	  // cout << "f = " << f << endl;
	  // cout << "g = " << g << endl;
	  // cout << "DT^p = " << DT << "^" << p << endl;
	  // cout << "CoeffVecWRT(f, DT) = " << CoeffVecWRT(f, DT) << endl;
	  // cout << "result = " << Hcoeff(g) * f - CoeffVecWRT(f, DT)[p] * power(DT, p - Hp_g) * g << endl;

	  RingElem gcoeff = Hcoeff(g);
	  RingElem fcoeff = CoeffVecWRT(f, DT)[p];
	  RingElem gcd_coeff = gcd(fcoeff, gcoeff);

	  // cout << "gcd = " << gcd_coeff << endl;

	  // return Hcoeff(g) * f - CoeffVecWRT(f, DT)[p] * power(DT, p - Hp_g) * g;
	  return (gcoeff * f - fcoeff * power(DT, p - Hp_g) * g) / gcd_coeff;
	}
      }
    }

    return f;
  }

  RingElem partial_rem(RingElem r, const std::vector<RingElem> A, bool full_reduction = false)
  {
    if (IsZero(r)) return r;

    for (unsigned int i=0; i < A.size(); i++) {
      RingElem rem = partial_rem(r, A[i], full_reduction);
      if (rem != r) {
	// std::cerr << r << " % " << A[i] << " = " << rem << endl;
	r = rem;
	i = -1;
	if (IsZero(r)) return r;
      }
    }

    // return r/CoeffEmbeddingHom(owner(r))(content(r));
    return r;
  }

  std::vector<RingElem> partial_rem(const std::vector<RingElem> v, const std::vector<RingElem> A, bool full_reduction = false)
  {
    std::vector<RingElem> result = v;
    for (auto&& e : result) {
      e = partial_rem(e, A, full_reduction);
    }
    return result;
  }

  RingElem rem(ConstRefRingElem f, ConstRefRingElem g)
  {
    return partial_rem(f, g, true);
  }

  RingElem rem(RingElem r, const std::vector<RingElem> A)
  {
    return partial_rem(r, A, true);
  }

  std::vector<RingElem> rem(const std::vector<RingElem> v, const std::vector<RingElem> A)
  {
    return partial_rem(v, A, true);
  }

  std::vector<RingElem> Union(const std::vector<RingElem> a, const std::vector<RingElem> b)
  {
    std::vector<RingElem> c = a;
    c.insert(c.end(), b.begin(), b.end());
    return c;
  }

  std::vector<RingElem> Union(const std::vector<RingElem> a, const RingElem b)
  {
    std::vector<RingElem> c = a;
    c.push_back(b);
    return c;
  }

  RingElem initial(ConstRefRingElem f)
  {
    std::pair<PPMonoidElem, int> hdt_hp = HDT_Hp(f);
    long idx;
    RingElem result(owner(f));

    CoCoA_ASSERT(IsIndet(idx, hdt_hp.first));

    for (auto it=BeginIter(f); !IsEnded(it); ++it) {
      if (exponent(PP(it), idx) == hdt_hp.second) {
	result += monomial(owner(f), coeff(it), PP(it));
      }
    }

    result /= monomial(owner(f), 1, IndetPower(owner(LPP(f)), idx, hdt_hp.second));
    return result;
  }

  RingElem separant(ConstRefRingElem f)
  {
    PPMonoidElem hdt = HDT(f);
    long idx;
    RingElem result(owner(f));

    CoCoA_ASSERT(IsIndet(idx, hdt));

    for (auto it=BeginIter(f); !IsEnded(it); ++it) {
      long exp = exponent(PP(it), idx);
      if (exp > 0) {
	result += exp * monomial(owner(f), coeff(it), PP(it)) / monomial(owner(f), 1, hdt);
      }
    }

    return result;
  }

  std::vector<RingElem> initials_and_separants(std::vector<RingElem> A)
  {
    std::vector<RingElem> result;

    for (RingElem f : A) {
      const RingElem i = initial(f);
      if (! IsConstant(i)) result.push_back(i);
      const RingElem s = separant(f);
      if (! IsConstant(s)) result.push_back(s);
    }

    return result;
  }

  bool is_autoreduced(const std::vector<RingElem> set)
  {
    for (unsigned int i=0; i < set.size(); i++) {
      for (unsigned int j=i+1; j < set.size(); j++) {
	if (rem(set[i], set[j]) != set[i]) return false;
      }
    }
    return true;
  }

  // An autoreduced subset C of a set E is called a characteristic set
  // of E if E contains no non-zero elements reduced w.r.t. C.  If
  // they're reduced, then taking their remainder doesn't change them.

  bool is_characteristic_set(const std::vector<RingElem> set, const std::vector<RingElem> subset)
  {
    for (auto eq: set) {
      // if (! IsZero(rem(eq, subset))) return false;
      // cerr << "eq: " << eq << endl;
      // cerr << "rem(eq,subset): " << rem(eq,subset) << endl;
      if (rem(eq, subset) == eq) return false;
    }
    return is_autoreduced(subset);
  }

  // This is a native version of Rosenfeld_Groebner that uses an older
  // algorithm ([Bo95] instead of [Bo09]), doesn't work right and is
  // very slow.

#if 0
  void Rosenfeld_Groebner(std::vector<RingElem> equations, std::vector<RingElem> inequations,
			  std::vector<RegularSystem> & results, int nesting_level=0)
  {
    // check for "obvious" inconsistencies

    for (auto eq: equations) {
      if (IsConstant(eq) && !IsZero(eq)) return;
    }
    for (auto ineq: inequations) {
      if (IsZero(ineq)) return;
    }

    std::cerr << "Rosenfeld_Groebner(" << nesting_level << "): equations: " << equations << " inequations: " << inequations << endl;

    // Produce A, a characteristic set of 'equations'.
    //
    // Another definition of characteristic sets is that they are the
    // autoreduced sets of lowest rank ([Ko73] p. 81).  Since rank of
    // a set is dominated by the lowest ranking elements, this leads
    // directly to an algorithm for computing characteristic sets.
    //
    // At each step, select the element of lowest rank and add it to
    // A.  Remove from consideration all elements that can now be
    // reduced by A.  Repeat until no elements remain.
    //
    // XXX Can be sped up by not using full partial reduction; only
    // enough to see if reduction is possible.

    std::vector<RingElem> A;

    std::vector<RingElem> remaining = equations;

    // XXX could we use pointers or references to avoid copying here?

    // XXX I think we can just use 'equations' directly, unless we want it to be const

    std::sort(remaining.begin(), remaining.end(), [] (ConstRefRingElem a, ConstRefRingElem b) -> bool {
	return rank(a) > rank(b);
      });

    while (remaining.size() > 0) {
      A.push_back(remaining.back());
      remaining.erase(std::remove_if(remaining.begin(), remaining.end(), [&A](ConstRefRingElem o) -> bool {
	    return (o != rem(o, A));
	  }), remaining.end());
    }

    // remove A's elements from 'equations'

    for (auto element: A) {
      equations.erase(find(equations.begin(), equations.end(), element));
    }

    // simplify A by dividing out by each polynomial's content

    std::for_each(A.begin(), A.end(), [](RingElem &r) {
	r /= CoeffEmbeddingHom(owner(r))(content(r));
      });

    std::cerr << "A: " << A << endl;
    std::cerr << "rank: " << rank(A) << endl;
    //std::cerr << "equations: " << equations << endl;
    //std::cerr << "Union(A,equations): " << Union(A,equations) << endl;

    CoCoA_ASSERT(is_characteristic_set(Union(A, equations), A));

    std::vector<RingElem> h = initials_and_separants(A);
    std::vector<RingElem> R = rem(Union(equations, Dpoly(A)), A); // (remaining equations union D-polynomials of A) rem A

    // Having zeros in R creates problems later, so remove them now.
    // Specifically, elements in R get passed recursively as
    // 'equations', which lets zeros creep into our characteristic
    // sets, so then we're trying to reduce by zero, and that fails
    // asserts in the reduction code.

    // This strikes me as unnecessarily verbose, but it's too terse to even compile.
    //   std::remove_if(R.begin(), R.end(), IsZero);
    // We need to either specify the template arguments...
    //   std::remove_if<std::vector<RingElem>::iterator, bool(ConstRefRingElem)> (R.begin(), R.end(), IsZero);
    // or use a lambda function
    //   std::remove_if(R.begin(), R.end(), [](const RingElem &o) { return IsZero(o);} );
    // but even that doesn't actually remove the elements, we have to call erase() to do that.
    //   http://stackoverflow.com/questions/7958216

    R.erase(std::remove_if(R.begin(), R.end(), [](const RingElem &o) { return IsZero(o);} ), R.end());

    std::for_each(R.begin(), R.end(), [](RingElem &r) {
	r /= CoeffEmbeddingHom(owner(r))(content(r));
      });

    std::cerr << "R: " << R << endl;
    std::cerr << "h: " << h << endl;

    if (R.empty() || IsZero(R[0])) {
      RegularSystem result(A, Union(partial_rem(inequations, A), h));
      // std::cerr << "Got a result " << A << " \\ " << Union(partial_rem(inequations, A), h) << endl;
      if (! IsElem(one(RegularDifferentialIdeal::R), result.I)) {
	// std::cerr << "Pushed a result " << result << endl;
	results.push_back(result);
      }
    } else {
      Rosenfeld_Groebner(Union(A, R), Union(inequations, h), results, nesting_level+1);
    }

    while (! h.empty()) {
      RingElem hi = h.back();
      h.pop_back();
      Rosenfeld_Groebner(Union(Union(equations, A), hi), h, results, nesting_level+1);
    }
  }

  std::vector<RegularSystem> Rosenfeld_Groebner(std::vector<RingElem> equations, std::vector<RingElem> inequations)
  {
    std::vector<RegularSystem> results;

    Rosenfeld_Groebner(equations, inequations, results);

    return results;
  }
#endif

#if 0
  std::vector<RegularSystem> Rosenfeld_Groebner(std::vector<RingElem> inequations)
  {
    return Rosenfeld_Groebner(gens, inequations);
  }

  std::vector<RegularSystem> Rosenfeld_Groebner(void)
  {
    return Rosenfeld_Groebner(gens, std::vector<RingElem>());
  }
#endif

}

//using diffalg::Rosenfeld_Groebner;

#ifdef USE_BLAD

// This version of Rosenfeld_Groebner uses Francois Boulier's blad
// library.

/* The first issue is that the blad library requires dependent and
 * independent variables to be completely distinct, while we don't.
 * So we need to use a separate naming scheme for blad variables and
 * map them back and forth to our variables.  A single CoCoA
 * indeterminate can correspond to two blad variables, one to
 * differentiate by and one that can be differentiated.
 */

// To maintain maps of elements from different rings (i.e, a
// polynomial ring and its coefficient ring), we need to compare
// PPMonoidElem's from different PPMonoids, which normally isn't
// allowed, but we can fix that by ordering first on the address of
// the PPMonoid's indets vector.

struct globalCmp {
  bool operator()(ConstRefPPMonoidElem a, ConstRefPPMonoidElem b) const {
    if (owner(a) == owner(b)) {
      return (a < b);
    } else {
      return &indets(owner(a)) < &indets(owner(b));
    }
  }
};

class bladDifferentialRing;

class bladDifferentialRingBase {
  friend bladDifferentialRing;

  const ring R;
  const PowerPolyDifferentialRingBase * DR;

  bladDifferentialRingBase(const ring& R) : R(R),
					    DR(dynamic_cast<const PowerPolyDifferentialRingBase *>(R.myRawPtr()))
  { }

  bimap<PPMonoidElem, std::string, globalCmp> dependent_strings;
  bimap<PPMonoidElem, struct bav_variable *, globalCmp> dependent_vars;

  bimap<PPMonoidElem, std::string, globalCmp> independent_strings;
  bimap<PPMonoidElem, struct bav_variable *, globalCmp> independent_vars;

  const std::string next_name(const std::string& orig)
  {
    const std::regex re("^\\w*$");

    if (! std::regex_match(orig, re)) {
      return next_name("");
    }

    if ((orig.length() > 0) && (dependent_strings.count(orig) == 0) && (independent_strings.count(orig) == 0)) {
      return orig;
    }

    char next_name_str[3] = {'a', 'a', '\0'};

    while ((dependent_strings.count(orig + std::string(next_name_str)) != 0)
	   || (independent_strings.count(orig + std::string(next_name_str)) != 0)) {
      if (next_name_str[1] != 'z') {
	next_name_str[1] ++;
      } else {
	next_name_str[0] ++;
	next_name_str[1] = 'a';
      }
    }

    return orig + std::string(next_name_str);
  }

  std::vector<std::string> insert_symbols_into_dependent_strings(const ring& R)
  {
    if (symbols(R).size() == 0) return std::vector<std::string> {};
    if (IsFractionField(R)) {
      return insert_symbols_into_dependent_strings(BaseRing(R));
    } else {
      std::vector<std::string> results = insert_symbols_into_dependent_strings(CoeffRing(R));
      for (auto ind: indets(PPM(R))) {
	if (dependent_strings.count(ind) == 0) {
	  dependent_strings[ind] = next_name(head(Symbol(ind)));
	}
	results.push_back(dependent_strings[ind]);
      }
      return results;
    }
  }

public:

  // Construct a blad ordering using all of the indeterminates in a
  // PPMonoidElem, plus any symbols in our underlying coefficent ring.

  void blad_ordering(bav_Iordering * r)
  {
    // We current don't distinguish between independent and dependent
    // variables (perhaps we should), but blad does.  Thus, we run
    // through our indets, split them apart into base and deriv, and
    // save the derivations.  The independent variables will be the
    // derivations, while the dependent variables (blad's blocks)
    // will be the bases less the derivations.

    const PPMonoid & ppm(PPM(R));

    PPMonoidElem bases(ppm);
    PPMonoidElem derivations(ppm);

    for (auto ind: indets(ppm)) {
      PPMonoidElem base(owner(ind));
      PPMonoidElem deriv(owner(ind));
      split_differential_indet(ind, base, deriv);
      bases = radical(bases * base);
      derivations = radical(derivations * deriv);
    }

    // bases /= derivations;
    // CoCoA_ASSERT(IsOne(gcd(bases, derivations)));

    std::vector<PPMonoidElem> vbases;
    std::vector<PPMonoidElem> vderivations;

    for (auto ind: indets(ppm)) {
      if (IsDivisible(bases, ind)) {
	vbases.push_back(ind);
      }
      if (IsDivisible(derivations, ind)) {
	vderivations.push_back(ind);
      }
    }

    std::sort(vbases.begin(), vbases.end());
    std::sort(vderivations.begin(), vderivations.end());

    // blad wants everything in DECENDING order
    std::reverse(vbases.begin(), vbases.end());
    std::reverse(vderivations.begin(), vderivations.end());

    std::string blad_ordering = "ordering (derivations=[";

    for (auto ind: vderivations) {
      //vars2[ind] = head(Symbol(ind));
      if (independent_strings.count(ind) == 0) {
	independent_strings[ind] = next_name(head(Symbol(ind)));
      }
      blad_ordering += independent_strings[ind] + ",";
    }
    blad_ordering.pop_back();

    blad_ordering += "], blocks=[";

    if (DR->ranking != DifferentialRanking::lex) {
      blad_ordering += "[";
    }

    for (auto ind: vbases) {
      if (dependent_strings.count(ind) == 0) {
	dependent_strings[ind] = next_name(head(Symbol(ind)));
      }
      if (DR->ranking == DifferentialRanking::lex) {
	blad_ordering += "[" + dependent_strings[ind] + "],";
      } else {
	blad_ordering += dependent_strings[ind] + ",";
      }
    }

    if (DR->ranking != DifferentialRanking::lex) {
      blad_ordering.pop_back();
      blad_ordering += "],";
    }

    // If the underlying coefficient ring has any symbols, add
    // them at the end as the lowest ranking block

#if 0
    const auto coeff_syms = symbols(CoeffRing(PolyRing(R)));
    if (coeff_syms.size() != 0) {
      blad_ordering += "[";
      for (auto s: coeff_syms) {
	blad_ordering += head(s) + ",";
      }
      blad_ordering.pop_back();   // remove trailing comma
      blad_ordering += "],";
      insert_symbols_into_dependent_strings(CoeffRing(PolyRing(R)));
    }
#else
    const auto coeff_syms = insert_symbols_into_dependent_strings(CoeffRing(PolyRing(R)));
    if (coeff_syms.size() != 0) {
      blad_ordering += "[";
      for (auto s: coeff_syms) {
	blad_ordering += s + ",";
      }
      blad_ordering.pop_back();   // remove trailing comma
      blad_ordering += "],";
    }
#endif

    blad_ordering.pop_back();   // remove trailing comma
    blad_ordering += "])";

    // std::cerr << blad_ordering << endl;

    // blad doesn't use const qualifiers when declaring ba0_sscanf2()
    ba0_sscanf2 (const_cast<char *>(blad_ordering.c_str()), const_cast<char *>("%ordering"), r);
  }

  struct bav_variable * PPMonoidElem_to_blad_variable(ConstRefPPMonoidElem ind, bool independent = false)
  {
    struct bav_variable * result;

    if (independent) {

      if (independent_vars.count(ind) == 1) {
	return independent_vars[ind];
      }
      ba0_sscanf2 (const_cast<char *>(independent_strings.at(ind).c_str()), const_cast<char *>("%v"), &result);

      independent_vars[ind] = result;

    } else {

      if (dependent_vars.count(ind) == 1) {
	return dependent_vars[ind];
      }

      PPMonoidElem base(owner(ind));
      PPMonoidElem deriv(owner(ind));

      split_differential_indet(ind, base, deriv);

      if (IsOne(deriv)) {
	//ba0_sscanf2 (const_cast<char *>(head(Symbol(base)).c_str()), const_cast<char *>("%v"), &result);
	ba0_sscanf2 (const_cast<char *>(dependent_strings.at(base).c_str()), const_cast<char *>("%v"), &result);
      } else {
	result = PPMonoidElem_to_blad_variable(base);
	for (auto ind2: indets(owner(ind))) {
	  while (IsDivisible(deriv, ind2)) {
	    bav_term t;
	    bav_init_term(&t);
	    bav_set_term_variable(&t, PPMonoidElem_to_blad_variable(ind2, true), 1);
	    result = bav_diff2_variable(result, &t);

	    deriv /= ind2;
	  }
	}
      }

      // Stash the result in a bimap so we convert back the reverse way later.
      dependent_vars[ind] = result;
    }

    return result;
  }

  void RingElem_to_blad_polynomial(ConstRefRingElem x, struct bap_polynom_mpz * const result)
  {
    // struct bap_polynom_mpz result;

    // bap_init_polynom_mpz(&result);

    // Assumes that result has already been initialized.

    if (IsZZ(owner(x))) {
      BigInt N;
      CoCoA_ASSERT(IsInteger(N, x));
      bap_init_polynom_one_mpz(result);
      bap_mul_polynom_numeric_mpz(result, result, mpzref(N));
      return;
    }

    if (IsFractionField(owner(x))) {
      CoCoA_ASSERT(IsOne(den(x)));
      RingElem_to_blad_polynomial(num(x), result);
      return;
    }

    CoCoA_ASSERT(IsSparsePolyRing(owner(x)));

    for (auto it=BeginIter(x); !IsEnded(it); ++it) {
      bav_term t;
      bav_init_term(&t);
      for (long idx=0; idx < NumIndets(owner(PP(it))); ++idx) {
	long exp = exponent(PP(it), idx);
	//bav_mul_term_variable(&t, &t, vars[indet(owner(PP(it)), idx)], exp);
	if (exp != 0) {
	  bav_mul_term_variable(&t, &t, PPMonoidElem_to_blad_variable(indet(owner(PP(it)), idx)), exp);
	}
      }
      struct bap_polynom_mpz c;
      bap_init_polynom_mpz(&c);
      RingElem_to_blad_polynomial(coeff(it), &c);
      bap_mul_polynom_term_mpz(&c, &c, &t);
      bap_add_polynom_mpz(result, result, &c);
    }
  }

  void RingElems_to_blad_polynomial_table(std::vector<RingElem> v, struct bap_tableof_polynom_mpz * const result)
  {
    struct ba0_list * L = nullptr;
    for (auto e: v) {
      struct bap_polynom_mpz * P = bap_new_polynom_mpz();
      RingElem_to_blad_polynomial(e, P);
      L = ba0_endcons_list(P, L);
    }
    ba0_set_table_list((ba0_table *) result, L);
    //ba0_printf(const_cast<char *>("%t[%Az]\n"), result);
  }

  RingElem blad_variable_to_RingElem(struct bav_variable * v, const ring & R)
  {
    // We must have previously constructed this variable from a
    // PPMonoidElem, so retrieve the stashed value.
    if (dependent_vars.count(v) == 1) {
      PPMonoidElem elem = dependent_vars.at(v);
      if (owner(elem) == PPM(R)) {
	return monomial(R, 1, elem);
      } else if (owner(elem) == PPM(CoeffRing(R))) {
	return CoeffEmbeddingHom(R)(monomial(CoeffRing(R), 1, elem));
      } else {
	throw "mistake";
      }
    } else {
      return monomial(R, 1, independent_vars.at(v));
    }
  }

  // This function is DESTRUCTIVE of the bav_term!

  RingElem blad_term_to_RingElem(struct bav_term * t, const ring & R)
  {
    RingElem result(one(R));

    while (! bav_is_one_term(t)) {
      result *= power(blad_variable_to_RingElem(bav_leader_term(t), R), bav_leading_degree_term(t));
      bav_shift_term(t, t);
    }

    return result;
  }

  // Converting polynomials to RingElem's using the technique of
  // printing them to a string and then parsing them is problematic.
  // For one thing, we can have higher derivatives in the polynomial
  // than (currently) appear in the CoCoA ring, which would trigger a
  // parse error.  Also, if the ordering of the derivatives gets mixed
  // up, we might try to parse something like f_{tx} when we really
  // wanted f_{xt}.  Therefore, we iterate through the polynomial's
  // monomials, iterate through each variable in the monomial, and
  // construct each derivative using multideriv().

  RingElem blad_polynomial_to_RingElem(struct bap_polynom_mpz * P)
  {
    struct bap_itermon_mpz iter;
    RingElem result(zero(R));

    // loop over monomials
    bap_begin_itermon_mpz (&iter, P);
    while (! bap_outof_itermon_mpz (&iter)) {
      struct bav_term T;
      mpz_t *c;

      bav_init_term(&T);
      bap_term_itermon_mpz (&T, &iter);
      c = bap_coeff_itermon_mpz (&iter);

      RingElem monomial(R, *c);

      // loop over variables
      while (! bav_is_one_term(&T)) {
	bav_variable * lt = bav_leader_term(&T);
	bav_Idegree ltd = bav_leading_degree_term(&T);

	bav_variable * lt_base = bav_order_zero_variable(lt);
	bav_term diff_op;

	bav_init_term(&diff_op);
	bav_operator_between_derivatives(&diff_op, lt_base, lt);

	RingElem base = blad_variable_to_RingElem(lt_base, R);
	RingElem deriv = blad_term_to_RingElem(&diff_op, R);

	monomial *= power(diffalg::multideriv(base, LPP(deriv)), ltd);

	bav_exquo_term_variable(&T, &T, lt, ltd);
      }

      result += monomial;

      bap_next_itermon_mpz (&iter);
    }

    return result;
  }

};

// static map<const PowerPolyDifferentialRingBase *, bladDifferentialRingBase *> bDRcache;

/* I'd like to keep libblad context around, but that turns out to be
 * problematic.  My library's ability to create new differential
 * variables on the fly wrecks havoc with libblad, because we can't
 * change orderings without creating brand new variables that aren't
 * compatible with the ones we had before.
 *
 * So we create a new libblad context every time we create a
 * bladDifferentialRing.
 */

class bladDifferentialRing {
  bladDifferentialRingBase * bDRbase;

public:
  bladDifferentialRing(const ring& R)
  {
    bDRbase = new bladDifferentialRingBase(R);
    bad_restart(0,0);
    bav_R_init();

    bav_Iordering r;

    bDRbase->blad_ordering(&r);
    bav_R_push_ordering (r);

  }

  ~bladDifferentialRing()
  {
    bad_terminate(ba0_done_level);
  }

  //const bladDifferentialRingBase* operator->() const { return bDRbase; }  ///< Allow const member fns to be called.
  bladDifferentialRingBase* operator->() const { return bDRbase; }  ///< Allow all member fns to be called.
};

// We need to first map our indeterminates to the strings we use to
// name them in the blad ordering, then set the ordering, then
// create the blad variables and map them back and forth.

std::vector<RegularSystem> Rosenfeld_Groebner(std::vector<RingElem> equations, std::vector<RingElem> inequations)
{
  // std::cerr << "Rosenfeld_Groebner equations: " << equations << " inequations: " << inequations << endl;

  struct bap_tableof_polynom_mpz * eqns;
  struct bap_tableof_polynom_mpz * ineqs;

  struct bad_intersectof_regchain * T;
  std::vector<RegularSystem> results;

  // XXX assumes that equations is not an empty vector
  const ring& R = owner(equations[0]);
  bladDifferentialRing bdr(R);

  // memory management (see blad docs §2.2.4.1)
  // This creates problems when new CoCoA objects are created during the session.

  // struct ba0_mark M;
  // ba0_record(&M);

  eqns = (struct bap_tableof_polynom_mpz *) ba0_new_table ();
  ineqs = (struct bap_tableof_polynom_mpz *) ba0_new_table ();

#if 0
  ba0_sscanf2 (const_cast<char *>(RingElems_to_blad_string(equations).c_str()), const_cast<char *>("%t[%Az]"), eqns);
  ba0_sscanf2 (const_cast<char *>(RingElems_to_blad_string(inequations).c_str()), const_cast<char *>("%t[%Az]"), ineqs);
#else
  bdr->RingElems_to_blad_polynomial_table(equations, eqns);
  bdr->RingElems_to_blad_polynomial_table(inequations, ineqs);
#endif

  // ba0_printf(const_cast<char *>("Equations: %t[%Az]\n"), eqns);
  // ba0_printf(const_cast<char *>("Inequations: %t[%Az]\n"), ineqs);

  // We have to initialize this structure with a list of desired properties.

  T = bad_new_intersectof_regchain ();
  ba0_sscanf2
    (const_cast<char *> ("intersectof_regchain ([], [differential, squarefree, coherent, primitive, autoreduced])"),
     const_cast<char *> ("%intersectof_regchain"), T);

  bad_base_field * bf = bad_new_base_field();

  if ((equations.size() == 4) && (IsZero(equations[2]))) {
    struct bav_tableof_variable * generators = (struct bav_tableof_variable *) ba0_new_table();
    struct bav_tableof_parameter * parameters = (struct bav_tableof_parameter *) ba0_new_table();

    ba0_sscanf2 (const_cast<char *>("[a]"), const_cast<char *>("%t[%v]"), generators);
    ba0_sscanf2 (const_cast<char *>("[a]"), const_cast<char *>("%t[%param]"), parameters);

    bad_set_base_field_generators_and_relations(bf, generators, 0, parameters, false, false);

    // ba0_printf (const_cast<char *>("Base field: %base_field\n"), bf);
  }

  bad_Rosenfeld_Groebner(T, eqns, ineqs, bf, (struct bad_splitting_control *)0);

  bad_remove_redundant_components_intersectof_regchain(T, T, bf);

  // ba0_printf (const_cast<char *>("%intersectof_regchain\n"), T);

  for (int i=0; i < T->inter.size; i ++) {
    std::vector<RingElem> v;

    struct bad_regchain * chain = (struct bad_regchain *)(T->inter.tab[i]);
    // ba0_printf (const_cast<char *>("%regchain\n"), chain);

    for (int j=0; j < chain->decision_system.size; j ++) {
      struct bap_polynom_mpz * P = chain->decision_system.tab[j];

      // If e is a derivative of something in our coefficent ring,
      // it will translate to zero, since our coefficients are
      // assumed to be constants, but libblad doesn't know that.

      RingElem e = bdr->blad_polynomial_to_RingElem(P);
      //std::cout << e << endl;
      if (!IsZero(e)) {
	v.push_back(e);
      }
    }
    results.push_back(RegularSystem(v, std::vector<RingElem>()));
  }

  //ba0_restore(&M);

  //bav_R_pull_ordering();

  return results;
}

// Passing RingElem's directly to Rosenfeld_Groebner constructs them
// into a vector of equations, with no inequations.

template<typename... Args>
std::vector<RegularSystem> Rosenfeld_Groebner(Args... rest) {
  std::vector<RingElem> gens{rest...};
  return Rosenfeld_Groebner(gens, std::vector<RingElem>());
}

// Compute reduced form of a RingElem against a RegularSystem
//
// Reduced form of p is a differential polynomial q, equivalent to p
// module the ideal defined by rs.

RingElem operator% (RingElem e, const RegularSystem& rs)
{
  // XXX assumes that equations is not an empty vector
  bladDifferentialRing bdr(owner(rs.equations[0]));

  struct bap_tableof_polynom_mpz * eqns;
  bad_regchain * A = bad_new_regchain();

  eqns = (struct bap_tableof_polynom_mpz *) ba0_new_table ();
  bdr->RingElems_to_blad_polynomial_table(rs.equations, eqns);

#if 0
  ba0_sscanf2(const_cast<char *> ("regchain ([], [differential, squarefree, coherent, primitive, autoreduced])"),
	      const_cast<char *> ("%regchain"), A);
#endif

  bad_set_and_extend_regchain_tableof_polynom_mpz(A, eqns, 0, 0, false, false);

  struct bap_polynom_mpz * F = bap_new_polynom_mpz();

  bdr->RingElem_to_blad_polynomial(e, F);

  struct bap_ratfrac_mpz * R = bap_new_ratfrac_mpz();

  bad_reduced_form_polynom_mod_regchain(R, F, 0, A);

  return bdr->blad_polynomial_to_RingElem(& (R->numer));
}
#endif

#if 0
std::ostream & operator<<(std::ostream &out, const RegularDifferentialIdeal ideal)
{
  out << "RegularDifferentialIdeal(";

  std::copy(ideal.gens.begin(), ideal.gens.end() - 1, std::ostream_iterator<RingElem>(out, ", "));

  out << ideal.gens.back() << ")";

  return out;
}
#endif

void testRegularDifferentialIdeal(void)
{
  ring ZZ = RingZZ();
  ring QQ = RingQQ();

  // ExponentRing - these are the indeterminates that can appear in powers

  //ring ExponentRing = NewOrderedPolyRing(ZZ, vector<symbol> {symbol("p")});
  //RingElem p(ExponentRing, "p");

  // We now create a K[Z[p]] ring whose coefficient and exponent rings are ExponentRing,
  // along with its fraction field.

  // Current implementation of RegularDifferentialIdeal can't handle polynomials in exponents, so this doesn't work
  //PPMonoid PPM = NewPPMonoidRing(vector<string> {"f", "t", "q", "x", "y", "z"}, lex, ExponentRing);
  //ring R = NewPowerPolyDifferentialRing(ExponentRing, PPM);

  // This doesn't work.  Neither for_each nor transform can return a range; it has to be passed in.
  //std::vector<symbol> symbols = std::for_each(names.begin(), names.end(), [](string &n){ new symbol(n); });

  // This SEGVs.  Not sure why.
  //std::vector<string> names {"f", "t", "q", "x", "y", "z"};
  //std::vector<symbol> symbols;
  //std::transform(names.begin(), names.end(), symbols.begin(), [](const string &n){ return symbol(n); });

  //std::vector<symbol> symbols={symbol("f"), symbol("t"), symbol("q"), symbol("u"), symbol("v"), symbol("x"), symbol("y"), symbol("z")};

  //PPMonoid PPM = NewPPMonoid(symbols, lex);

  // NYI: Can't create ideals with coefficents not in a field
  //ring R = NewPolyRing(ExponentRing, PPM);

  //ring R = NewPolyRing(QQ, PPM);
  //ring K = NewFractionField(R);

  // Most significant first: f > t > q > u > v > z > y > x

  PPMonoid PPM = NewPPMonoidRing(vector<string> {"f", "t", "q", "u", "v", "z", "y", "x"}, lex, QQ);
  ring R = NewPowerPolyDifferentialRing(QQ, PPM);

  RingElem f(R, "f");
  RingElem t(R, "t");
  RingElem q(R, "q");
  RingElem u(R, "u");
  RingElem v(R, "v");
  RingElem x(R, "x");
  RingElem y(R, "y");
  RingElem z(R, "z");

  //cerr << colon(ideal(4*power(x,3) -4*power(x,2) -3*x), ideal(R, std::vector<RingElem>{12*power(x,2) -8*x -3, 3595744*power(x,3) -2711796*power(x,2) -1338412*x +3747, -1923560*power(x,4) +4510032*power(x,3) -923526*power(x,2) -1407034*x -109815, 19488*power(x,5) -81280*power(x,4) +92256*power(x,3) -528*power(x,2) -23858*x -1920, -4480*power(x,6) +17760*power(x,5) -23360*power(x,4) +8992*power(x,3) +3864*power(x,2) -2374*x -240, 2048*power(x,7) -3584*power(x,6) -4512*power(x,5) +7840*power(x,4) +2336*power(x,3) -2712*power(x,2) -894*x -36, -2*power(x,2) +4*x -2, -24*power(x,3) +40*power(x,2) -10*x -6, 4*power(x,2) -4*x -3, zero(R), 144*power(x,4) -192*power(x,3) -8*power(x,2) +48*x +9})) << endl;
  //cerr << colon(ideal(4*power(x,3) -4*power(x,2) -3*x), ideal(R, std::vector<RingElem>{4*power(x,2) -4*x -3, 144*power(x,4) -192*power(x,3) -8*power(x,2) +48*x +9})) << endl;
  //cerr << colon(ideal(4*power(x,3) -4*power(x,2) -3*x), ideal(R, std::vector<RingElem>{4*power(x,2) -4*x -3, 48*x +9})) << endl;

  //exit(1);

  //cerr << colon(ideal(4*x*x*x-4*x*x-3*x), ideal(zero(R), 4*x*x-4*x-3)) << endl;
  //exit(1);

  // CoCoA forbids mixed ring operations. 'p' in ExponentRing is
  // different from 'p' in K, which we now create.

  // RingElem pK = EmbeddingHom(K)(CoeffEmbeddingHom(R)(p));

  // Funny syntax - the derivatives don't exist in the ring until we
  // call 'deriv', so we do that first and then check to see if
  // they're right, instead of the other way around.

  RingElem ft = deriv(f,t);
  RingElem ftt = deriv(deriv(f,t),t);
  RingElem fqt = deriv(deriv(f,q),t);

  RingElem fx = deriv(f,x);
  RingElem fy = deriv(f,y);
  RingElem fz = deriv(f,z);

  RingElem ux = deriv(u,x);
  RingElem uxy = deriv(ux, y);
  RingElem vy = deriv(v,y);
  RingElem vx = deriv(v,x);
  RingElem vxx = deriv(vx,x);

  //RegularDifferentialIdeal di(f, t);

  RingElem xt = deriv(x,t);
  RingElem yt = deriv(y,t);

  RingElem xtt = deriv(xt,t);
  RingElem ytt = deriv(yt,t);

  RingElem xttt = deriv(xtt,t);
  RingElem xtttt = deriv(xttt,t);

  // First example from [Bo95]

  RingElem s1 = (2*xtt+1)*yt+y;
  RingElem s2 = xt*xt + x;

  // Test some of our basic differential algebra functions

  CoCoA_ASSERT(monomial(R, 1, diffalg::rank(s1)) == yt);
  CoCoA_ASSERT(diffalg::initial(s1) == (2*xtt+1));
  CoCoA_ASSERT(diffalg::separant(s1) == (2*xtt+1));

  CoCoA_ASSERT(monomial(R, 1, diffalg::rank(s2)) == xt*xt);
  CoCoA_ASSERT(diffalg::initial(s2) == 1);
  CoCoA_ASSERT(diffalg::separant(s2) == 2*xt);

  //RegularDifferentialIdeal Bo95_ex1(s1, s2);

  // To check this better, we should use unordered_set's rather than
  // vectors (because the ordering really doesn't matter), then compare
  // the unordered_sets.  The check below depends on the order that
  // the algorithm produces the equations, which is irrelevant.

  std::vector<RegularSystem> Bo95_ex1_RG = Rosenfeld_Groebner(s1, s2);
  CoCoA_ASSERT(Bo95_ex1_RG.size() == 2);
  CoCoA_ASSERT(Bo95_ex1_RG[0].equations.size() == 2);
  CoCoA_ASSERT(Bo95_ex1_RG[1].equations.size() == 2);
  CoCoA_ASSERT(Bo95_ex1_RG[0].equations[0] == (xt*xt + x));
  CoCoA_ASSERT(Bo95_ex1_RG[0].equations[1] == y);
  CoCoA_ASSERT(Bo95_ex1_RG[1].equations[0] == x);
  CoCoA_ASSERT(Bo95_ex1_RG[1].equations[1] == yt + y);
  //std::cerr << Bo95_ex1_RG << endl;

  // Check reduced forms
  CoCoA_ASSERT((y % Bo95_ex1_RG[0]) == 0);
  CoCoA_ASSERT((yt % Bo95_ex1_RG[1]) == -(y % Bo95_ex1_RG[1]));
  //std::cerr << (yt % Bo95_ex1_RG[1]) << endl;

  // Third example from [Bo95]

  // ((ẍ − ẋ² + 1)x⁽⁴⁾ − x⁽³⁾² + (3 ẋẍ + 2xẋ² − ẋ − 2x)x⁽³⁾ − 2ẍ³ +
  // ((−6x − 2) ẋ − 1)ẍ² + ((2x + 2) ẋ³ + ẋ² + (−2x − 2) ẋ + 2)ẍ − ẋ² + 1)

  // ((ẍ − ẋ² + 1)y − x⁽³⁾ + ( ẋ + 2x)ẍ − ẋ)

  // ((ẍ − ẋ² + 1)z − x⁽³⁾ + ( ẋ + 2x)ẍ − ẋ)

  auto Bo95_ex3_RG = Rosenfeld_Groebner(y*xt+yt+1-xtt, 2*yt*x+2*xt*y*x+y+2*x-ytt, z-y);
  CoCoA_ASSERT(Bo95_ex3_RG.size() == 1);
  CoCoA_ASSERT(Bo95_ex3_RG[0].equations.size() == 3);
  CoCoA_ASSERT(Bo95_ex3_RG[0].equations[0] == (xtt-xt*xt+1)*xtttt - xttt*xttt + (3*xt*xtt + 2*x*xt*xt - xt - 2*x)*xttt - 2*xtt*xtt*xtt + ((-6*x-2)*xt-1)*xtt*xtt + ((2*x+2)*xt*xt*xt + xt*xt + (-2*x-2)*xt+2)*xtt - xt*xt + 1);
  CoCoA_ASSERT(Bo95_ex3_RG[0].equations[1] == (xtt-xt*xt+1)*y - xttt + (xt+2*x)*xtt - xt);
  CoCoA_ASSERT(Bo95_ex3_RG[0].equations[2] == (xtt-xt*xt+1)*z - xttt + (xt+2*x)*xtt - xt);
  //std::cerr << Bo95_ex3_RG << endl;

  //std::cerr << RegularDifferentialIdeal(R, std::vector<RingElem> {19488*power(x,5) -81280*power(x,4) +92256*power(x,3) -528*power(x,2) -23858*x -1920,  48720*power(x,4) -162560*power(x,3) +138384*power(x,2) -528*x -11929,  194880*power(x,3) -487680*power(x,2) +276768*x -528,  194880*power(x,3) -487680*power(x,2) +276768*x -528,  584640*power(x,2) -975360*x +276768,  -3299572*power(x,4) +11636124*power(x,3) -7331817*power(x,2) -4491566*x -243840,  -13198288*power(x,3) +34908372*power(x,2) -14663634*x -4491566,  -13198288*power(x,3) +34908372*power(x,2) -14663634*x -4491566,  -39594864*power(x,2) +69816744*x -14663634,  -105936364564*power(x,4) +287632590324*power(x,3) -102295429833*power(x,2) -83289434024*x -4781890560,  -423745458256*power(x,3) +862897770972*power(x,2) -204590859666*x -83289434024,  -423745458256*power(x,3) +862897770972*power(x,2) -204590859666*x -83289434024,  -1271236374768*power(x,2) +1725795541944*x -204590859666,  -423745458256*power(x,3) +862897770972*power(x,2) -204590859666*x -83289434024,  -1271236374768*power(x,2) +1725795541944*x -204590859666,  -1271236374768*power(x,2) +1725795541944*x -204590859666,  -2542472749536*x +1725795541944}).Rosenfeld_Groebner();

  //std::cerr << RegularDifferentialIdeal(R, std::vector<RingElem> {power(x,5) -power(x,4) +power(x,3) -power(x,2) -x -1920,  power(x,4) -power(x,3) +power(x,2) -x -11,  power(x,3) -power(x,2) +x -5,  power(x,2) -x +27,  power(x,2) +4*x -14, -26*x +14}).Rosenfeld_Groebner();

  // BLAD's test rg0

  // by putting a constant in the equations we should generate an inconsistent system
  CoCoA_ASSERT(Rosenfeld_Groebner(ux*ux - 4*u, uxy*vy - u + 1, RingElem(R, 421), vxx - ux).size() == 0);

  //CoCoA_ASSERT(Rosenfeld_Groebner(ux*ux - 4*u, uxy*vy - u + 1, vxx - ux).size() == 3);
  //std::cerr << Rosenfeld_Groebner(ux*ux - 4*u, uxy*vy - u + 1, vxx - ux) << endl;

  // Sixth example from Mansfield thesis

  // std::cerr << Rosenfeld_Groebner(power(fx,2) - 1, power(fy,2) - 1, (fx+fy)*fz - 1) << endl;
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
  return make_pair(MinI, MinJ);
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
  return make_pair(MinI, MinJ);
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

/* For testing purposes, we want to create some CoCoA matrices from
 * arrays in C.  C++ template tricks let us do this with functions
 * that take references to arbitrary sized arrays.
 */

template<int n> matrix NewMatrixFromC(ring R, int (&cmat)[n][n])
{
  matrix M(NewDenseMat(R,n,n));
  
  for (int i=0; i < n; ++i)
    for (int j=0; j < n; ++j)
      SetEntry(M, i, j, cmat[i][j]);
  return M;
}

template<int n> matrix NewDiagMatrixFromC(ring R, int (&cmat)[n])
{
  matrix M(NewDenseMat(R,n,n));
  
  for (int i=0; i < n; ++i)
    SetEntry(M, i, i, cmat[i]);
  return M;
}

// Give a test matrix and its invariant factors, check Smith Normal Form

template<int n> void testOneSmithFactor(int (&M_int)[n][n], int (&F_int)[n])
{
  matrix M_Z(NewMatrixFromC(RingZZ(), M_int));
  matrix F_Z(NewDiagMatrixFromC(RingZZ(), F_int));

  // cout << M_Z << endl;

  SmithRecord L = SmithFactor(M_Z);

  //cout << L.M << endl;
  //cout << L.U << endl;
  //cout << L.V << endl;
  //cout << L.U * M_Z * L.V << endl;

  CoCoA_ASSERT(L.M == F_Z);
  CoCoA_ASSERT(L.U * M_Z * L.V == L.M);
}

void testSmithFactor(void)
{
  // two test matrices and their invariant factors

  int matrix1[3][3] = {{ 2, 4, 4},
		       {-6, 6, 12},
		       {10,-4,-16}};
  int factors1[3] = {2,6,12};

  int matrix2[4][4] = {{-6, 111, -36, 6},
		       { 5,-672, 210, 74},
		       { 0,-255,  81, 24},
		       {-7, 255, -81,-10}};
  int factors2[4] = {1,3,21,0};

  testOneSmithFactor(matrix1, factors1);
  testOneSmithFactor(matrix2, factors2);
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
 * A vector of pointers to Differentials is provided.  All must
 * operate on the same ring.  The symbols provided to the Weyl algebra
 * must also exist in that ring.  We assume that CanonicalHom can lift
 * coefficients from the Weyl algebra's coefficient ring into the
 * target ring of the differentials.
 *
 * We use POINTERS to Differentials because Differentials can be
 * modified (maybe not a very good idea)...
 */

class WeylOperatorAlgebra : public RingWeylImpl
{
private:
  const std::vector<Differential *> differentials;
  const std::vector<symbol> SymList;
public:
  WeylOperatorAlgebra(const ring& CoeffRing, const std::vector<symbol>& names, const std::vector<long>& ElimIndets, const std::vector<Differential *>& differentials):
    RingWeylImpl(CoeffRing, names, ElimIndets),
    differentials(differentials),
    SymList(names)
  {
    const ring R = owner(*differentials[0]);
    for (auto it=differentials.begin(); it != differentials.end(); ++it) {
      if (owner(**it) != R) {
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

    if (owner(y) != owner(*differentials[0])) {
      CoCoA_ERROR(ERR::MixedRings, "WeylOperatorAlgebra used on wrong ring");
    }
    for (auto it=myBeginIter(rawx); !IsEnded(it); ++it) {
      RingElem term(y);
      // each term in our operator is in canonical form, so operators are applied first
      for (long idx=0; idx < myNumTrueIndets; ++idx) {
	const long Didx = idx + myNumTrueIndets;
	const long d = exponent(PP(it), Didx);

	for (long i=1; i <= d; ++i) {
	  term = (*differentials[idx])(term);
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

    ring K = owner(*differentials[0]);
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
	    target = (*differentials[idx])(target);
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

  // Transform diffop into a recursion relationship on the solution
  // polynomial's coefficients, returning 0 if we can determine that
  // the recursion can't finitely terminate, and would thus be
  // unsuitable for a polynomial.

  RingElem poly_solve(ConstRefRingElem diffop) const
  {
    // diffop is a RingElem in the WeylOperatorAlgebra
    ring WA = owner(diffop);
    long myNumTrueIndets = NumIndets(WA)/2;

    // The coefficients will involve powers of the exponents, i.e, x^m
    // differentiates into m x^(m-1) and thus 'm' is in the
    // coefficient.  Create CoeffRing with an indet (like 'm') for
    // each indet in the Weyl algebra.

    // XXX we're currently restricted to integer coefficients

    ring CoeffRing = NewPolyRing(RingZZ(), myNumTrueIndets);

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

#if 0
  RingElem rat_solve(ConstRefRingElem diffop) const
  {
    // diffop is a RingElem in the WeylOperatorAlgebra
    ring WA = owner(diffop);
    long myNumTrueIndets = NumIndets(WA)/2;

  }
#endif

};

const WeylOperatorAlgebra* WeylOperatorAlgebraPtr(const ring& R)
{
  return dynamic_cast<const WeylOperatorAlgebra*>(R.myRawPtr());
}

SparsePolyRing NewWeylOperatorAlgebra(const ring& CoeffRing, const std::vector<symbol>& names, const std::vector<Differential *>& differentials)
{
  std::vector<long> ElimIndets;   // empty set
  return SparsePolyRing(new WeylOperatorAlgebra(CoeffRing, WANAMES(names), ElimIndets, differentials));
}

/* Return the smallest exponent with which an indeterminate appears in a polynomial.
 *
 * WARNING!  Exponents can be incomparable.  For example, we can't
 * tell if f^p has smaller power than f^2 without knowing something
 * about p.
 */

RingElem minExponent(RingElem in, RingElem indet)
{
  RingElem poly;
  long index;
  bool valid = false;
  RingElem minexp;

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
    if (valid && !IsInteger(myexp - minexp)) {
      CoCoA_ERROR(ERR::NYI, "incompatible exponents in minExponent");
    }
    if (!valid || (myexp < minexp)) {
      minexp = myexp;
      valid = true;
    }
  }

  return minexp;
}

/* Returns the factor by which an indeterminate appears to its minimal power.
 *
 * WARNING!  Exponents can be incomparable.  For example, we can't
 * tell if f^p has smaller power than f^2 without knowing something
 * about p.
 */

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

  bool valid = false;
  RingElem minexp;
  RingElem result(owner(poly));

  for (SparsePolyIter it=BeginIter(poly); !IsEnded(it); ++it) {
    RingElem myexp = RingElemExponent(PP(it), index);
    if (valid && !IsInteger(myexp - minexp)) {
      CoCoA_ERROR(ERR::NYI, "incompatible exponents in minCoeff");
    }
    if (!valid || (myexp < minexp)) {
      minexp = myexp;
      valid = true;
      result = 0;
    }
    if (myexp == minexp) {
      result += monomial(owner(poly), coeff(it), PP(it)/power(indet(owner(PP(it)), index), minexp));
    }
  }

  return result;
}

/* This variant of minCoeff doesn't require the indeterminate to be
 * specified, instead searching all of the indeterminates to find
 * those with monomial minimum coefficients.  It isn't as useful as I
 * thought it would be, because we need to limit our attention to
 * those indeterminates that are irreducible, and we currently have no
 * way to tell which they are.
 */

void minCoeff2(RingElem in)
{
  vector<RingElem> targets;

  if (IsFractionField(owner(in))) {
    targets = indets(BaseRing(owner(in)));
  } else {
    targets = indets(owner(in));
  }

  for (auto it=targets.begin(); it != targets.end(); ++ it) {
    RingElem result = minCoeff(in, *it);
    if (IsMonomial(result)) {
      cout << "minCoeff(eq," << *it << ") = " << result << endl;
    }
  }
}

void program()
{
  cout << boolalpha; // so that bools print out as true/false
  cout << TeX;

  const ring ZZ = RingZZ();
  const ring QQ = RingQQ();

  // ExponentRing - these are the indeterminates that can appear in powers

  const ring ExponentRing = NewOrderedPolyRing(ZZ, vector<symbol> {symbol("p"), symbol("a"), symbol("i"), symbol("b"), symbol("c")});
  const RingElem p(ExponentRing, "p");
  const RingElem a(ExponentRing, "a");
  const RingElem b(ExponentRing, "b");
  const RingElem c(ExponentRing, "c");
  const RingElem i(ExponentRing, "i");

  // We now create a K[Z[p]] ring whose coefficient and exponent rings are ExponentRing,
  // along with its fraction field.

  const PPMonoid PPM = NewPPMonoidRing(vector<string> {"x", "t", "z", "r", "T", "T_t", "(t+1)",
	"f", "f_x", "f_{xx}", "f_t", "q", "q_x", "q_{xx}", "q_t",
	"n", "n_{x}", "n_{xx}", "n_{t}",
	"n_i", "n_{ix}", "n_{ixx}", "n_{it}",
	"n_r", "n_{rx}", "n_{rxx}", "n_{rt}",
	"n_e", "n_{ex}", "n_{exx}", "n_{et}",
	"d_e", "d_{ex}", "d_{exx}", "d_{et}",
	"N", "N_x", "N_{xx}", "N_t", "D", "D_x", "D_{xx}", "D_t"}, lex, ExponentRing);
  const ring R = NewPowerPolyRing(ExponentRing, PPM);
  const ring K = NewFractionField(R);

  // x,t are in our field of definition
  // z = exp(-x^2/(4(t+1)))
  // r = sqrt(t)

  const RingElem x(K, "x");
  const RingElem t(K, "t");
  const RingElem z(K, "z");
  const RingElem r(K, "r");

  const RingElem tpo(K, "(t+1)");

  // T is a polynomial in C[t] (doesn't involve x or z)
  const RingElem T(K, "T");
  const RingElem Tt(K, "T_t");

  // N is the numerator
  const RingElem N(K, "N");
  const RingElem Nx(K, "N_x");
  const RingElem Nxx(K, "N_{xx}");
  const RingElem Nt(K, "N_t");

  // D is the denominator
  const RingElem D(K, "D");
  const RingElem Dx(K, "D_x");
  const RingElem Dxx(K, "D_{xx}");
  const RingElem Dt(K, "D_t");

  // f and q are factors of something.  Typically D=f^p q,
  // where f is irreducible and q is coprime to f.
  const RingElem f(K, "f");
  const RingElem fx(K, "f_x");
  const RingElem fxx(K, "f_{xx}");
  const RingElem ft(K, "f_t");

  const RingElem q(K, "q");
  const RingElem qx(K, "q_x");
  const RingElem qxx(K, "q_{xx}");
  const RingElem qt(K, "q_t");

  const RingElem n(K, "n");
  const RingElem n_x(K, "n_{x}");
  const RingElem n_xx(K, "n_{xx}");
  const RingElem n_t(K, "n_{t}");

  const RingElem n_r(K, "n_r");
  const RingElem n_rx(K, "n_{rx}");
  const RingElem n_rxx(K, "n_{rxx}");
  const RingElem n_rt(K, "n_{rt}");

  // numerator and denominator of exponent in z = e^(n/d)

  const RingElem n_e(K, "n_e");
  const RingElem n_ex(K, "n_{ex}");
  const RingElem n_exx(K, "n_{exx}");
  const RingElem n_et(K, "n_{et}");

  const RingElem d_e(K, "d_e");
  const RingElem d_ex(K, "d_{ex}");
  const RingElem d_exx(K, "d_{exx}");
  const RingElem d_et(K, "d_{et}");

  // n_i is one coefficient in a numerator sum
  const RingElem n_i(K, "n_i");
  const RingElem n_ix(K, "n_{ix}");
  const RingElem n_ixx(K, "n_{ixx}");
  const RingElem n_it(K, "n_{it}");

  // setup our differentials (acting on K)

  Differential dx(K, vector<RingHom> {x >> 1, t >> 0, r >> 0, tpo >> 0,
	N >> Nx, Nx >> Nxx, D >> Dx, Dx >> Dxx, T >> 0,
	f >> fx, fx >> fxx, q >> qx, qx >> qxx,
	n >> n_x, n_x >> n_xx,
	n_e >> n_ex, n_ex >> n_exx,
	n_r >> n_rx, n_rx >> n_rxx,
	n_i >> n_ix, n_ix >> n_ixx,
	d_e >> d_ex, d_ex >> d_exx});

  Differential dt(K, vector<RingHom> {x >> 0, t >> 1, r >> r/(2*t), tpo >> 1,
	N >> Nt, D >> Dt, T >> Tt, f >> ft, q >> qt, n >> n_t, n_e >> n_et, n_r >> n_rt, n_i >> n_it, d_e >> d_et});

  RingElem z_exp = -x*x/(4*t);
  //RingElem z_exp = n_e/d_e;

  dx.update(z, z*dx(z_exp));
  dt.update(z, z*dt(z_exp));

  cout << dx(z) << endl;

  // Create a Weyl algebra, with ExponentRing as the coefficient ring.
  // I don't actually use operators with coefficients not in QQ, but WA.factor() currently won't work
  // unless the operator algebra and the target ring share the same coefficient ring.

  const ring WA = NewWeylOperatorAlgebra(ExponentRing,
					 vector<symbol> {symbol("x"), symbol("t")},
					 vector<Differential *> {&dx, &dt});

  const RingElem WA_x(WA, "x");
  const RingElem WA_t(WA, "t");
  const RingElem WA_dx(WA, "dx");
  const RingElem WA_dt(WA, "dt");

  // our operator

  const RingElem O = WA_dx*WA_dx - WA_dt;

  // cout << WeylOperatorAlgebraPtr(WA)->poly_solve(O) << endl;
  // cout << WeylOperatorAlgebraPtr(WA)->poly_solve(-2*WA_x*WA_dx + 2*WA_t*WA_dx*WA_dx -2*WA_t*WA_dt -1) << endl;

  RingElem e = N/D;

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
  cout << "try N/Tq where T_x = 0" << endl;
  cout << endl;

  eq = num(O*(N/(T*q)));

  cout << eq << endl;
  cout << "minCoeff(eq, T) = " <<minCoeff(eq, T) << endl;



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
  cout << "eq - minCoeff(eq, r) = " << (eq = (eq - minCoeff(eq, r))/num(r)) << endl;
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



  cout << endl;
  cout << "try denominator z^p    with z = exp(n_e / d_e)" << endl;
  cout << endl;

  z_exp = n_e/d_e;

  dx.update(z, z*dx(z_exp));
  dt.update(z, z*dt(z_exp));

  d = power(z,p);
  eq = num(O*(N/d));
  cout << eq << endl;
  minCoeff2(eq);


  cout << endl;
  cout << "try denominator z^p    with z = exp(n_e / f^a d_e)" << endl;
  cout << endl;

  z_exp = n_e/(power(f,a)*d_e);

  dx.update(z, z*dx(z_exp));
  dt.update(z, z*dt(z_exp));

  d = power(z,p);
  eq = num(O*(N/d));
  cout << eq << endl;
  //minCoeff2(eq);


  cout << endl;
  cout << "try denominator z^p    with z = exp(n_e / T^a d_e)" << endl;
  cout << endl;

  z_exp = n_e/(power(T,a)*d_e);

  dx.update(z, z*dx(z_exp));
  dt.update(z, z*dt(z_exp));

  d = power(z,p);
  eq = num(O*(N/d));
  cout << eq << endl;
  minCoeff2(eq);


  cout << endl;
  cout << "try denominator z^p f^a q   with z = exp(n_e / (f^b d_e))" << endl;
  cout << endl;

  z_exp = n_e/(power(f,b)*d_e);

  dx.update(z, z*dx(z_exp));
  dt.update(z, z*dt(z_exp));

  d = power(z,p) * power(f,a) * q;
  //d = power(z,p) * q;
  eq = num(O*(N/d));
  //cout << eq << endl;
  //eq = (t >> (tpo - 1)) (CanonicalHom(R,K)(eq));
  cout << eq << endl;
  cout << "minCoeff(eq, f) = " << minCoeff(eq, f) << endl;
  //cout << "minCoeff(eq, D) = " << minCoeff(eq, D) << endl;


  cout << endl;
  cout << "try denominator z^p T^a q   with z = exp(n_e / (T^b d_e)) where T_x = 0" << endl;
  cout << endl;

  z_exp = n_e/(power(T,b)*d_e);

  dx.update(z, z*dx(z_exp));
  dt.update(z, z*dt(z_exp));

  d = power(z,p) * power(T, a) * q;
  eq = num(O*(N/d));
  //cout << eq << endl;
  //eq = (t >> (tpo - 1)) (CanonicalHom(R,K)(eq));
  cout << eq << endl;
  cout << "minCoeff(eq, T) = " << minCoeff(eq, T) << endl;
  cout << "eq - minCoeff(eq, T) = " << (eq - num(T)*minCoeff(eq, T))/num(T) << endl;
  minCoeff2((eq - num(T)*minCoeff(eq, T))/num(T));
  // XXX this factor() call doesn't work
  //cout << factor(minCoeff(eq,T)).myFactors() << endl;
  //cout << "minCoeff(eq, t+1) = " << minCoeff(eq, tpo) << endl;
  //cout << den(dx(dx(N/d)) - dt(N/d)) << endl;

  cout << endl;
  cout << "try denominator z^p T^a q   with z = exp(n_e / (T d_e)) where T_x = 0" << endl;
  cout << endl;

  z_exp = n_e/(T*d_e);

  dx.update(z, z*dx(z_exp));
  dt.update(z, z*dt(z_exp));

  d = power(z,p) * power(T, a) * q;
  eq = num(O*(N/d));
  //cout << eq << endl;
  //eq = (t >> (tpo - 1)) (CanonicalHom(R,K)(eq));
  cout << eq << endl;
  cout << "minCoeff(eq, T) = " << minCoeff(eq, T) << endl;
  cout << "eq - minCoeff(eq, T) = " << (eq - num(T)*minCoeff(eq, T))/num(T) << endl;
  minCoeff2((eq - num(T)*minCoeff(eq, T))/num(T));
  // XXX this factor() call doesn't work
  //cout << factor(minCoeff(eq,T)).myFactors() << endl;
  //cout << "minCoeff(eq, t+1) = " << minCoeff(eq, tpo) << endl;
  //cout << den(dx(dx(N/d)) - dt(N/d)) << endl;

  cout << endl;
  cout << "try denominator z^p   with z = exp(n_e / (f^a d_e))" << endl;
  cout << endl;

  z_exp = n_e/(power(f,a)*d_e);

  dx.update(z, z*dx(z_exp));
  dt.update(z, z*dt(z_exp));

  //cout << "dx(z) = " << dx(z) << endl;
  //cout << "dt(z) = " << dt(z) << endl;

  d = power(z,p);
  eq = num(O*(N/d));
  //cout << eq << endl;
  //eq = (t >> (tpo - 1)) (CanonicalHom(R,K)(eq));
  //cout << "NUM: " << eq << endl;
  //cout << "DEN: " << den(O*(N/d)) << endl;
  //cout << minExponent(eq, tpo) << endl;
  cout << "minCoeff(eq, f) = " << minCoeff(eq, f) << endl;
  //minCoeff2(eq);
  //cout << "minCoeff(eq, T) = " << minCoeff(eq, T) << endl;
  //cout << den(dx(dx(N/d)) - dt(N/d)) << endl;

#if 0
  cout << endl;
  cout << "try denominator z^p with numerator d_e^b N" << endl;
  cout << endl;

  eq = num(O*((power(d_e,b) * N)/power(z,p)));
  cout << eq << endl;
  cout << "minCoeff(eq, d_e) = " << minCoeff(eq, d_e) << endl;
#endif

}

void program2()
{
  cout << boolalpha; // so that bools print out as true/false
  cout << TeX;

  const ring ZZ = RingZZ();
  const ring QQ = RingQQ();

  // ExponentRing - these are the indeterminates that can appear in powers

  const ring ExponentRing = NewOrderedPolyRing(ZZ, vector<symbol> {symbol("p"), symbol("a"), symbol("i"), symbol("b"), symbol("c")});
  const RingElem p(ExponentRing, "p");
  const RingElem a(ExponentRing, "a");
  const RingElem b(ExponentRing, "b");
  const RingElem c(ExponentRing, "c");
  const RingElem i(ExponentRing, "i");

  // We now create a K[Z[p]] ring whose coefficient and exponent rings are ExponentRing,
  // along with its fraction field.

  const PPMonoid PPM = NewPPMonoidRing(vector<string> {"x", "t", "z", "r", "T", "(t+1)",
	"f", "q", "n", "ni", "nr",
	"ne", "de", "N", "D"}, lex, ExponentRing);

  const ring R = NewPowerPolyDifferentialRing(ExponentRing, PPM, DifferentialRanking::grlexA);
  const ring K = NewFractionField(R);

  // x,t are in our field of definition
  // z = exp(-x^2/(4(t+1)))
  // r = sqrt(t)

  const RingElem x(K, "x");
  const RingElem t(K, "t");

  const RingElem z(K, "z");
  const RingElem z_x = deriv(z,x);
  const RingElem z_t = deriv(z,t);

  const RingElem r(K, "r");

  const RingElem tpo(K, "(t+1)");

  // T is a polynomial in C[t] (doesn't involve x or z)
  const RingElem T(K, "T");
  const RingElem Tt = deriv(T,t);
  const RingElem Tx = deriv(T,x);

  // N is the numerator
  const RingElem N(K, "N");
  const RingElem Nx = deriv(N,x);
  const RingElem Nxx = deriv(Nx,x);
  const RingElem Nt = deriv(N,t);

  // D is the denominator
  const RingElem D(K, "D");
  const RingElem Dx = deriv(D,x);
  const RingElem Dxx = deriv(Dx,x);
  const RingElem Dt = deriv(D,t);

  // f and q are factors of something.  Typically D=f^p q,
  // where f is irreducible and q is coprime to f.
  const RingElem f(K, "f");
  const RingElem fx = deriv(f,x);
  const RingElem fxx = deriv(fx,x);
  const RingElem ft = deriv(f,t);

  const RingElem q(K, "q");
  const RingElem qx = deriv(q,x);
  const RingElem qxx = deriv(qx,x);
  const RingElem qt = deriv(q,t);

  const RingElem n(K, "n");
  const RingElem n_x = deriv(n,x);
  const RingElem n_xx = deriv(n_x,x);
  const RingElem n_t = deriv(n,t);

  const RingElem n_r(K, "nr");
  const RingElem n_rx = deriv(n_r,x);
  const RingElem n_rxx = deriv(n_rx,x);
  const RingElem n_rt = deriv(n_r,t);

  // numerator and denominator of exponent in z = e^(n/d)

  const RingElem n_e(K, "ne");
  const RingElem n_ex = deriv(n_e,x);
  const RingElem n_exx = deriv(n_ex,x);
  const RingElem n_et = deriv(n_e,t);

  const RingElem d_e(K, "de");
  const RingElem d_ex = deriv(d_e,x);
  const RingElem d_exx = deriv(d_ex,x);
  const RingElem d_et = deriv(d_e,t);

  // n_i is one coefficient in a numerator sum
  const RingElem n_i(K, "ni");
  const RingElem n_ix = deriv(n_i,x);
  const RingElem n_ixx = deriv(n_ix,x);
  const RingElem n_it = deriv(n_i,t);

  // setup our differentials (acting on K)

  BasicDifferential dx(x);
  BasicDifferential dt(t);

  // Create a Weyl algebra, with ExponentRing as the coefficient ring.
  // I don't actually use operators with coefficients not in QQ, but WA.factor() currently won't work
  // unless the operator algebra and the target ring share the same coefficient ring.

  const ring WA = NewWeylOperatorAlgebra(ExponentRing,
					 vector<symbol> {symbol("x"), symbol("t")},
					 vector<Differential *> {&dx, &dt});

  const RingElem WA_x(WA, "x");
  const RingElem WA_t(WA, "t");
  const RingElem WA_dx(WA, "dx");
  const RingElem WA_dt(WA, "dt");

  // our operator

  const RingElem O = WA_dx*WA_dx - WA_dt;

  // free exponential

  // q = f^a
  //RingElem z_exp = n_e/(power(f,a) * d_e);
  // q_x = a f^(a-1) f_x = a q f_x / f
  // q_t = a f^(a-1) f_t = a q f_t / f

  RingElem z_exp = n_e/(q * d_e);
  //RingElem e = N/power(z,p);
  RingElem e = N/(D*q);

  cout << O*e << endl;

  //auto RG = Rosenfeld_Groebner(num(z_x - z*deriv(z_exp,x)), num(z_t - z*deriv(z_exp,t)));

  //auto RG = Rosenfeld_Groebner(num(z_x - z*deriv(z_exp,x)), num(z_t - z*deriv(z_exp,t)),
  //			       num(qx - CanonicalHom(ExponentRing, K)(a) * q * fx / f),
  //			       num(qt - CanonicalHom(ExponentRing, K)(a) * q * ft / f));

  //auto RG = Rosenfeld_Groebner(num(qx - CanonicalHom(ExponentRing, K)(a) * q * fx / f),
  //			       num(qt - CanonicalHom(ExponentRing, K)(a) * q * ft / f));

  auto RG = Rosenfeld_Groebner(num(qx - CanonicalHom(ExponentRing, K)(a) * q * fx / f),
			       num(qt - CanonicalHom(ExponentRing, K)(a) * q * ft / f),
			       num(deriv(CanonicalHom(ExponentRing, K)(a), x)),
			       num(deriv(CanonicalHom(ExponentRing, K)(a), t)));

  for (auto s: RG) {
    cout << s << endl;
    //cout << num(O*e) % s << endl;
  }

  //cout << dx(z) << endl;

}

//----------------------------------------------------------------------
// Use main() to handle any uncaught exceptions and warn the user about them.
int main()
{
  GlobalManager CoCoAFoundations;

  try
  {
    testPowerPolyRing();
    testPowerPolyDifferentialRing();
    testRegularDifferentialIdeal();
    testSmithFactor();
    testDiophantineSolvable();
    program2();
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
