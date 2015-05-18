//   Copyright (c)  2015 Brent Baccala

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


// PPMonoidNested
//
// A PPMonoid that expands an underlying PPMonoid by adding indeterminates.
//
// Used to create modules from an underlying ring.

#include "CoCoA/matrix.H"
#include "CoCoA/MatrixView.H"
#include "CoCoA/OpenMath.H"
#include "CoCoA/PPMonoidNested.H"

using namespace std;

namespace CoCoA
{

// File local inline functions

inline PPMonoidNestedElem * PPMonoidNestedImpl::myElem(RawPtr rawpp) const
{
  return static_cast<PPMonoidNestedElem *>(rawpp.myRawPtr());
}


inline const PPMonoidNestedElem * PPMonoidNestedImpl::myElem(ConstRawPtr rawpp) const
{
  return static_cast<const PPMonoidNestedElem *>(rawpp.myRawPtr());
}


bool PPMonoidNestedImpl::myCheckExponents(const std::vector<long>& expv) const
{
  // Check expv.size == myNumIndets.
  // Check exps are non-neg
  if (len(expv) != myNumIndets) return false;
  for (long i=0; i < myNumIndets; ++i)
    if (expv[i] < 0) return false;
  return true;
}


  PPOrdering NewNestedOrdering(PPOrdering nestedOrdering, long GradingDim, long NumIndets);

  namespace PPOrd
  {
    class NestedImpl : public PPOrderingBase
    {
    private:
      friend PPOrdering CoCoA::NewNestedOrdering(PPOrdering nestedOrdering, long GradingDim, long NumIndets);
      NestedImpl(PPOrdering nestedOrdering, long GradingDim, long NumIndets);
      NestedImpl(const NestedImpl&);            ///< NEVER DEFINED -- copy ctor disabled
      NestedImpl& operator=(const NestedImpl&); ///< NEVER DEFINED -- assignment disabled
    public:
      virtual void myOutputSelf(std::ostream& out) const;
      virtual void myOutputSelf(OpenMathOutput& OMOut) const;
      virtual void myOrdMatCopy(MatrixView& M) const;
      virtual bool IamStdGraded() const {return false;}
    };

    NestedImpl::NestedImpl(PPOrdering nestedOrdering, long GradingDim, long NumIndets):
      PPOrderingBase(NumIndets, GradingDim)
    {}


    void NestedImpl::myOutputSelf(std::ostream& out) const
    {
      out << "PPOrderingNested(" << myNumIndets << ")";
    }


    void NestedImpl::myOutputSelf(OpenMathOutput& OMOut) const
    {
      OMOut->mySendApplyStart();
      OMOut << OpenMathSymbol("cocoa", "PPOrderingNested");
      OMOut << myNumIndets;
      OMOut->mySendApplyEnd();
    }


    void NestedImpl::myOrdMatCopy(MatrixView& M) const
    {
      CoCoA_ASSERT(NumRows(M) == myNumIndets);
      CoCoA_ASSERT(NumCols(M) == myNumIndets);
      //??? There should be a better way than this!
      AssignZero(M);
      for (long i=0; i < myNumIndets; ++i)
        SetEntry(M, i, i, 1);
    }

  }

  PPOrdering NewNestedOrdering(PPOrdering nestedOrdering, long GradingDim, long NumIndets)
  {
    if (NumIndets <= 0) CoCoA_ERROR(ERR::BadArg, "NewNestedOrdering(PPOrdering,GradingDim,NumIndets)");
    return PPOrdering(new PPOrd::NestedImpl(nestedOrdering, GradingDim, NumIndets));
  }

//----   Constructors & destructor   ----//

  PPMonoidNestedImpl::PPMonoidNestedImpl(const PPMonoid& nestedPPM, const std::vector<symbol>& IndetNames, long GradingDim, ModOrdTypeForcing MOType):
    PPMonoidBase(NewNestedOrdering(ordering(nestedPPM), GradingDim, len(IndetNames)), IndetNames),
    MOType(MOType),
    nestedPPM(nestedPPM),
    numNestedIndets(len(nestedPPM->myIndets())),
    numExtraIndets(len(IndetNames) - numNestedIndets),
    myIndetVector()
{
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


PPMonoidNestedImpl::~PPMonoidNestedImpl()
{}

/////////////////////////////////////////////////////////////////////////////



const std::vector<PPMonoidElem>& PPMonoidNestedImpl::myIndets() const
{
  return myIndetVector;
}


const PPMonoidElem& PPMonoidNestedImpl::myOne() const
{
  return *myOnePtr;
}


PPMonoidElemRawPtr PPMonoidNestedImpl::myNew() const
{
  PPMonoidElemRawPtr rawpp(new PPMonoidNestedElem(nestedPPM->myNew(), numExtraIndets));
  myAssignOne(rawpp); // cannot throw
  return rawpp;
}

PPMonoidElemRawPtr PPMonoidNestedImpl::myNew(PPMonoidElemConstRawPtr rawcopypp) const
{
  PPMonoidElemRawPtr rawpp(new PPMonoidNestedElem(nestedPPM->myNew(), numExtraIndets));
  myAssign(rawpp, rawcopypp); // cannot throw
  return rawpp;
}


PPMonoidElemRawPtr PPMonoidNestedImpl::myNew(const std::vector<long>& expv) const
{
  CoCoA_ASSERT(myCheckExponents(expv));
  PPMonoidElemRawPtr rawpp(new PPMonoidNestedElem(nestedPPM->myNew(), numExtraIndets));
  myAssign(rawpp, expv); // cannot throw
  return rawpp;
}


void PPMonoidNestedImpl::myAssignOne(RawPtr rawpp) const
{
  PPMonoidNestedElem * const elem = myElem(rawpp);
  nestedPPM->myAssignOne(elem->nestedRawElem);
  for (long i = 0; i < numExtraIndets; ++i)
    elem->exponents[i] = 0;
}


void PPMonoidNestedImpl::myAssign(RawPtr rawpp, ConstRawPtr rawpp1) const
{
  if (rawpp == rawpp1) return;

  PPMonoidNestedElem * const elem = myElem(rawpp);
  const PPMonoidNestedElem * const elem1 = myElem(rawpp1);

  nestedPPM->myAssign(elem->nestedRawElem, elem1->nestedRawElem);
  elem->exponents = elem1->exponents;
}

void PPMonoidNestedImpl::myAssign(RawPtr rawpp, const vector<long>& expv1) const
{
  CoCoA_ASSERT(myCheckExponents(expv1));

  PPMonoidNestedElem * const elem = myElem(rawpp);
  vector<long> nested_expv(expv1);

  nested_expv.resize(numNestedIndets);
  nestedPPM->myAssign(elem->nestedRawElem, nested_expv);

  for (long i = 0; i < numExtraIndets; ++i)
    elem->exponents[i] = expv1[numNestedIndets + i];
}


void PPMonoidNestedImpl::myDelete(RawPtr rawpp) const
{
  PPMonoidNestedElem * const elem = myElem(rawpp);
  nestedPPM->myDelete(elem->nestedRawElem);
  delete elem;
}


void PPMonoidNestedImpl::mySwap(RawPtr rawpp1, RawPtr rawpp2) const
{
  if (rawpp1 == rawpp2) return;

  PPMonoidNestedElem * const elem1 = myElem(rawpp1);
  PPMonoidNestedElem * const elem2 = myElem(rawpp2);

  elem1->exponents.swap(elem2->exponents);
}


void PPMonoidNestedImpl::myMul(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
{
  // No worries about aliasing.
  PPMonoidNestedElem * const elem = myElem(rawpp);
  const PPMonoidNestedElem * const elem1 = myElem(rawpp1);
  const PPMonoidNestedElem * const elem2 = myElem(rawpp2);

  nestedPPM->myMul(elem->nestedRawElem, elem1->nestedRawElem, elem2->nestedRawElem);
  for (long i=0; i < numExtraIndets; ++i)
    {
      elem->exponents[i] = elem1->exponents[i] + elem2->exponents[i];
    }
}


void PPMonoidNestedImpl::myMulIndetPower(RawPtr rawpp, long indet, long exp) const  // assumes exp >= 0
{
  CoCoA_ASSERT(exp >= 0);
  CoCoA_ASSERT(0 <= indet && indet < myNumIndets);
  PPMonoidNestedElem * const elem = myElem(rawpp);
  if (indet < numNestedIndets) {
    nestedPPM->myMulIndetPower(elem->nestedRawElem, indet, exp);
  } else {
    elem->exponents[indet - numNestedIndets] += exp;
  }
}


void PPMonoidNestedImpl::myDiv(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
{
  // No worries about aliasing.
  PPMonoidNestedElem * const elem = myElem(rawpp);
  const PPMonoidNestedElem * const elem1 = myElem(rawpp1);
  const PPMonoidNestedElem * const elem2 = myElem(rawpp2);

  nestedPPM->myDiv(elem->nestedRawElem, elem1->nestedRawElem, elem2->nestedRawElem);

  for (long i=0; i < numExtraIndets; ++i)
    {
      CoCoA_ASSERT("Exponent Underflow" && elem1->exponents[i] >= elem2->exponents[i]);
      elem->exponents[i] = elem1->exponents[i] - elem2->exponents[i];
    }
}


void PPMonoidNestedImpl::myColon(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
{
  // No worries about aliasing.
  PPMonoidNestedElem * const elem = myElem(rawpp);
  const PPMonoidNestedElem * const elem1 = myElem(rawpp1);
  const PPMonoidNestedElem * const elem2 = myElem(rawpp2);

  nestedPPM->myColon(elem->nestedRawElem, elem1->nestedRawElem, elem2->nestedRawElem);

  for (long i=0; i < numExtraIndets; ++i)
    if (elem1->exponents[i] > elem2->exponents[i])
      elem->exponents[i] = elem1->exponents[i] - elem2->exponents[i];
    else
      elem->exponents[i] = 0;
}


void PPMonoidNestedImpl::myGcd(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
{
  // No worries about aliasing.
  PPMonoidNestedElem * const elem = myElem(rawpp);
  const PPMonoidNestedElem * const elem1 = myElem(rawpp1);
  const PPMonoidNestedElem * const elem2 = myElem(rawpp2);

  nestedPPM->myGcd(elem->nestedRawElem, elem1->nestedRawElem, elem2->nestedRawElem);

  for (long i = 0; i < numExtraIndets; ++i)
    elem->exponents[i] = min(elem1->exponents[i], elem2->exponents[i]);
}


void PPMonoidNestedImpl::myLcm(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
{
  // No worries about aliasing.
  PPMonoidNestedElem * const elem = myElem(rawpp);
  const PPMonoidNestedElem * const elem1 = myElem(rawpp1);
  const PPMonoidNestedElem * const elem2 = myElem(rawpp2);

  nestedPPM->myLcm(elem->nestedRawElem, elem1->nestedRawElem, elem2->nestedRawElem);

  for (long i = 0; i < numExtraIndets; ++i)
    elem->exponents[i] = max(elem1->exponents[i], elem2->exponents[i]);
}


void PPMonoidNestedImpl::myRadical(RawPtr rawpp, ConstRawPtr rawpp1) const
{
  PPMonoidNestedElem * const elem = myElem(rawpp);
  const PPMonoidNestedElem * const elem1 = myElem(rawpp1);

  nestedPPM->myRadical(elem->nestedRawElem, elem1->nestedRawElem);

  for (long i = 0; i < numExtraIndets; ++i)
    elem->exponents[i] = (elem1->exponents[i] > 0) ? 1 : 0;
}


void PPMonoidNestedImpl::myPowerSmallExp(RawPtr rawpp, ConstRawPtr rawpp1, long LongExp) const  // assumes exp >= 0
{
  CoCoA_ASSERT(LongExp >= 0);

  PPMonoidNestedElem * const elem = myElem(rawpp);
  const PPMonoidNestedElem * const elem1 = myElem(rawpp1);

  nestedPPM->myPowerSmallExp(elem->nestedRawElem, elem1->nestedRawElem, LongExp);

  for (long i = 0; i < numExtraIndets; ++i)
    elem->exponents[i] = LongExp * elem1->exponents[i];
}


bool PPMonoidNestedImpl::myIsOne(ConstRawPtr rawpp) const
{
  const PPMonoidNestedElem * const elem = myElem(rawpp);

  if (! nestedPPM->myIsOne(elem->nestedRawElem)) return false;

  for (long i = 0; i < numExtraIndets; ++i)
    if (elem->exponents[i] != 0) return false;

  return true;
}


bool PPMonoidNestedImpl::myIsIndet(long& index, ConstRawPtr rawpp) const
{
  const PPMonoidNestedElem * const elem = myElem(rawpp);
  long j = myNumIndets;

  if (! nestedPPM->myIsOne(rawpp)) {
    return nestedPPM->myIsIndet(j, elem->nestedRawElem);
  } else {
    for (long i = 0; i < numExtraIndets; ++i) {
      if (elem->exponents[i] == 0) continue;
      if (j != myNumIndets || (elem->exponents[i] != 1)) return false;
      j = i;
    }
    if (j == myNumIndets) return false;
    index = j;
    return true;
  }
}


bool PPMonoidNestedImpl::myIsCoprime(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
{
  const PPMonoidNestedElem * const elem1 = myElem(rawpp1);
  const PPMonoidNestedElem * const elem2 = myElem(rawpp2);

  if (! nestedPPM->myIsCoprime(elem1->nestedRawElem, elem2->nestedRawElem)) return false;

  for (long i = 0; i < numExtraIndets; ++i)
    if ((elem1->exponents[i] != 0) && (elem2->exponents[i] != 0)) return false;

  return true;
}


bool PPMonoidNestedImpl::myIsEqual(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
{
  const PPMonoidNestedElem * const elem1 = myElem(rawpp1);
  const PPMonoidNestedElem * const elem2 = myElem(rawpp2);

  if (! nestedPPM->myIsEqual(elem1->nestedRawElem, elem2->nestedRawElem)) return false;

  for (long i = 0; i < numExtraIndets; ++i)
    if (elem1->exponents[i] != elem2->exponents[i]) return false;
  return true;
}


bool PPMonoidNestedImpl::myIsDivisible(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
{
  const PPMonoidNestedElem * const elem1 = myElem(rawpp1);
  const PPMonoidNestedElem * const elem2 = myElem(rawpp2);

  if (! nestedPPM->myIsDivisible(elem1->nestedRawElem, elem2->nestedRawElem)) return false;

  for (long i = 0; i < numExtraIndets; ++i)
    if (elem1->exponents[i] < elem2->exponents[i]) return false;

  return true;
}


bool PPMonoidNestedImpl::myIsRadical(ConstRawPtr rawpp) const
{
  const PPMonoidNestedElem * const elem = myElem(rawpp);

  if (! nestedPPM->myIsRadical(elem->nestedRawElem)) return false;

  for (long i = 0; i < numExtraIndets; ++i)
    if (elem->exponents[i] > 1) return false;

  return true;
}


int PPMonoidNestedImpl::myCmp(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
{
  const PPMonoidNestedElem * const elem1 = myElem(rawpp1);
  const PPMonoidNestedElem * const elem2 = myElem(rawpp2);

  degree deg1(numExtraIndets - 1);
  degree deg2(numExtraIndets - 1);

  myWDeg(deg1, rawpp1);
  myWDeg(deg2, rawpp2);

  switch (MOType) {
  case PosWDegTO:
    if (elem1->exponents[numExtraIndets - 1] != elem2->exponents[numExtraIndets - 1]) {
      if (elem1->exponents[numExtraIndets - 1] > elem2->exponents[numExtraIndets - 1]) return 1; else return -1;
    }
    if (deg1 != deg2) {
      return cmp(deg1, deg2);
    }
    return nestedPPM->myCmp(elem1->nestedRawElem, elem2->nestedRawElem);

  case WDegTOPos:
    if (deg1 != deg2) {
      return cmp(deg1, deg2);
    }
    if (nestedPPM->myCmp(elem1->nestedRawElem, elem2->nestedRawElem) != 0) {
      return nestedPPM->myCmp(elem1->nestedRawElem, elem2->nestedRawElem);
    }
    if (elem1->exponents[numExtraIndets - 1] != elem2->exponents[numExtraIndets - 1]) {
      if (elem1->exponents[numExtraIndets - 1] > elem2->exponents[numExtraIndets - 1]) return 1; else return -1;
    }
    return 0;

  case WDegPosTO:; // This is the default
  default:;
    if (deg1 != deg2) {
      return cmp(deg1, deg2);
    }
    if (elem1->exponents[numExtraIndets - 1] != elem2->exponents[numExtraIndets - 1]) {
      if (elem1->exponents[numExtraIndets - 1] > elem2->exponents[numExtraIndets - 1]) return 1; else return -1;
    }
    if (nestedPPM->myCmp(elem1->nestedRawElem, elem2->nestedRawElem) != 0) {
      return nestedPPM->myCmp(elem1->nestedRawElem, elem2->nestedRawElem);
    }
    return 0;
  }

  CoCoA_ERROR(ERR::SERIOUS, "Unknown module ordering type in PPMonoidNested");
}


long PPMonoidNestedImpl::myStdDeg(ConstRawPtr rawpp) const
{
  const PPMonoidNestedElem * const elem = myElem(rawpp);
  long degree = nestedPPM->myStdDeg(elem->nestedRawElem);

  for (long i=0; i<numExtraIndets; ++i) {
    degree += elem->exponents[i];
  }

  return degree;
}


void PPMonoidNestedImpl::myWDeg(degree& d, ConstRawPtr rawpp) const
{
  const PPMonoidNestedElem * const elem = myElem(rawpp);

  //degree deg(GradingDim(nestedPPM));

  // Only two possibilities: we have shifts and numExtraIndets = GradingDim + 1, or we don't, and numExtraIndets = 1

  //if (numExtraIndets == GradingDim(nestedPPM) + 1) {
  if (numExtraIndets > 1) {
    degree shifts(GradingDim(nestedPPM));
    for (long i=0; i < numExtraIndets-1; ++i) {
      shifts.mySetComponent(i, elem->exponents[i]);
    }
    nestedPPM->myWDeg(d, elem->nestedRawElem);
    d += shifts;
  }

  // CoCoA_ERROR(ERR::NYI, "PPMonoid comparison in PPMonoidNested");
}


int PPMonoidNestedImpl::myCmpWDeg(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
{
  CoCoA_ERROR(ERR::NYI, "myCmpWDeg in PPMonoidNested");
}


int PPMonoidNestedImpl::myCmpWDegPartial(ConstRawPtr rawpp1, ConstRawPtr rawpp2, long i) const
{
  CoCoA_ERROR(ERR::NYI, "myCmpWDegPartial in PPMonoidNested");
}


long PPMonoidNestedImpl::myExponent(ConstRawPtr rawpp, long indet) const
{
  const PPMonoidNestedElem * const elem = myElem(rawpp);

  CoCoA_ASSERT(indet < myNumIndets);

  if (indet < numNestedIndets) {
    return nestedPPM->myExponent(elem->nestedRawElem, indet);
  }

  return elem->exponents[indet - numNestedIndets];
}

void PPMonoidNestedImpl::myBigExponent(BigInt& EXP, ConstRawPtr rawpp, long indet) const
{
  EXP = myExponent(rawpp, indet);
}


void PPMonoidNestedImpl::myExponents(std::vector<long>& v, ConstRawPtr rawpp) const
{
  const PPMonoidNestedElem * const elem = myElem(rawpp);
  std::vector<long> nested_exponents(numNestedIndets);

  CoCoA_ASSERT(len(v) == myNumIndets);

  nestedPPM->myExponents(nested_exponents, elem->nestedRawElem);

  for (long i=0; i < myNumIndets; ++i) {
    if (i < numNestedIndets) {
      v[i] = nested_exponents[i];
    } else {
      v[i] = elem->exponents[i - numNestedIndets];
    }
  }
}


void PPMonoidNestedImpl::myBigExponents(std::vector<BigInt>& expvector, ConstRawPtr rawpp) const
{
  const PPMonoidNestedElem * const elem = myElem(rawpp);
  std::vector<BigInt> nested_exponents(numNestedIndets);

  CoCoA_ASSERT(len(elem->exponents) == myNumIndets);

  nestedPPM->myBigExponents(nested_exponents, elem->nestedRawElem);

  for (long i=0; i < myNumIndets; ++i) {
    if (i < numNestedIndets) {
      expvector[i] = nested_exponents[i];
    } else {
      expvector[i] = elem->exponents[i - numNestedIndets];
    }
  }
}


void PPMonoidNestedImpl::myComputeDivMask(DivMask& dm, const DivMaskRule& DivMaskImpl, ConstRawPtr rawpp) const
{
  CoCoA_ERROR(ERR::NYI, "computeDivMask in PPMonoidNested");
  /* DivMaskImpl->myAssignFromExpv(dm, myElem(rawpp), myNumIndets); */
}


void PPMonoidNestedImpl::myOutputSelf(std::ostream& out) const
{
  out << "PPMonoidNested(" << numNestedIndets << "+" << numExtraIndets << ", " << myOrd <<")";
}


void PPMonoidNestedImpl::myDebugPrint(std::ostream& out, ConstRawPtr rawpp) const
{
  out << "DEBUG PP: myNumIndets=" << myNumIndets << ", exps=[";
  for (long i=0; i < myNumIndets; ++i)
    out << myExponent(rawpp, i) << " ";
  out << "]" << std::endl;
}

  //////////////////////////////////////////////////////////////////
  // Pseudo-ctors

  PPMonoid NewPPMonoidNested(const PPMonoid& nested, const std::vector<symbol>& IndetNames, long GradingDim, ModOrdTypeForcing MOType)
  {
    // Sanity check on the indet names given.

    if (!AreDistinct(IndetNames))
      CoCoA_ERROR(ERR::BadIndetNames, "NewPPMonoidNested(nestedPPM,IndetNames,GradingDim,MOType)");
    if (!AreArityConsistent(IndetNames))
      CoCoA_ERROR(ERR::BadIndetNames, "NewPPMonoidNested(nestedPPM,IndetNames,GradingDim,MOType)");

    return PPMonoid(new PPMonoidNestedImpl(nested, IndetNames, GradingDim, MOType));
  }

} // end namespace CoCoA
