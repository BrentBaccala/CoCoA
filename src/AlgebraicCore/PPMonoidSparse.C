//   Copyright (c)  2005,2007,2010  John Abbott

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


// Implementation of class PPMonoidSparseImpl

#include "CoCoA/PPMonoidSparse.H"
#include "CoCoA/assert.H"
#include "CoCoA/error.H"
#include "CoCoA/VectorOperations.H"
#include "CoCoA/symbol.H"
#include "CoCoA/utils.H"
#include "CoCoA/DivMask.H"
#include "CoCoA/MemPool.H"
#include "CoCoA/PPMonoid.H"
#include "CoCoA/PPOrdering.H"
#include "CoCoA/BigInt.H"

#include <algorithm>
using std::swap;
#include <iostream>
using std::ostream;
#include <list>
using std::list;
#include <memory>
#include <new>
//for placement new
#include <string>
using std::string;
#include <utility>
using std::pair;
#include <vector>
using std::vector;


namespace CoCoA
{

  namespace // anonymous namespace for file local definitions
  {
    struct IndexExp
    {
    public:
      IndexExp(long IndetIndex, int exp);
    public:
      long myIndetIndex;
      int myExp;
    };

    IndexExp::IndexExp(long IndetIndex, int exp):
        myIndetIndex(IndetIndex),
        myExp(exp)
    {}

    std::ostream& operator<<(std::ostream& out, const IndexExp& ie)
    {
      out << "IndexExp(" << ie.myIndetIndex << ", " << ie.myExp << ")";
      return out;
    }

  }


/////////////////////////////////////////////////////////////////////////////

  class PPMonoidSparseImpl: public PPMonoidBase
  {
  protected:
    friend PPMonoid NewPPMonoidSparse(const std::vector<symbol>& IndetNames, const PPOrdering& ord); // the only fn which calls the ctor
    PPMonoidSparseImpl(const std::vector<symbol>& IndetNames, const PPOrdering& ord);
    virtual ~PPMonoidSparseImpl() {};
  private: // disable copy ctor and assignment
    PPMonoidSparseImpl(const PPMonoidSparseImpl& copy);       ///< NEVER DEFINED -- disable default copy ctor
    PPMonoidSparseImpl& operator=(const PPMonoidSparseImpl&); ///< NEVER DEFINED -- disable assignment

  public:
    typedef PPMonoidElemRawPtr RawPtr;           ///< just to save typing
    typedef PPMonoidElemConstRawPtr ConstRawPtr; ///< just to save typing
  private:
    typedef std::list<IndexExp> value_t;
    static value_t& import(RawPtr rawpp);
    static const value_t& import(ConstRawPtr rawpp);

    const PPOrdering& myOrdering() const;
    virtual const std::vector<PPMonoidElem>& myIndets() const;            ///< std::vector whose n-th entry is n-th indet as PPMonoidElem

    virtual const PPMonoidElem& myOne() const;
    using PPMonoidBase::myNew;    // disable warnings of overloading
    virtual PPMonoidElemRawPtr myNew() const;                              ///< ctor from nothing
    virtual PPMonoidElemRawPtr myNew(PPMonoidElemConstRawPtr rawpp) const;    ///< ctor from another pp
    virtual PPMonoidElemRawPtr myNew(const std::vector<long>& expv) const; ///< ctor from exp vector
//NYI    virtual PPMonoidElemRawPtr myNew(const std::vector<BigInt>& v) const;     ///< ctor from exp vector

    virtual void myDelete(RawPtr rawpp) const;                                ///< dtor, frees pp
    virtual void mySwap(RawPtr rawpp1, RawPtr rawpp2) const;                     ///< swap(pp1, pp2)
    virtual void myAssignOne(RawPtr rawpp) const;                             ///< pp = 1
    virtual void myAssign(RawPtr rawpp, ConstRawPtr rawpp1) const;               ///< pp = pp1
    virtual void myAssign(RawPtr rawpp, const std::vector<long>& expv) const; ///< pp = v (assign from exp vector)
    virtual void myMul(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const; ///< pp = pp1*pp2
    using PPMonoidBase::myMulIndetPower;    // disable warnings of overloading
    virtual void myMulIndetPower(RawPtr rawpp, long indet, long exp) const;///< pp *= indet^exp
    virtual void myDiv(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const; ///< pp = pp1/pp2
    virtual void myColon(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const;///< pp = pp1/gcd(pp1,pp2)
    virtual void myGcd(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const; ///< pp = gcd(pp1,pp2)
    virtual void myLcm(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const; ///< pp = lcm(pp1,pp2)
    virtual void myRadical(RawPtr rawpp, ConstRawPtr rawpp1) const;           ///< pp = radical(pp1)
    virtual void myPowerSmallExp(RawPtr rawpp, ConstRawPtr rawpp1, long exp) const;// pp = pp1^exp (non-trivial), assumes exp >= 0
///NYI    virtual void myPowerBigExp(RawPtr rawpp, ConstRawPtr rawpp1, const BigInt& EXP) const;// pp = pp1^EXP (non-trivial), assumes EXP >= 0

    virtual bool myIsOne(ConstRawPtr rawpp) const;                            ///< is pp = 1?
    virtual bool myIsIndet(long& index, ConstRawPtr rawpp) const;             ///< true iff pp is an indet
    virtual bool myIsCoprime(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const;      ///< are pp1, pp2 coprime?
    virtual bool myIsEqual(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const;            ///< is pp1 equal to pp2?
    virtual bool myIsDivisible(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const;    ///< is pp1 is divisible by pp2?
    virtual bool myIsRadical(ConstRawPtr rawpp) const;                        ///< is pp equal to its radical?

    virtual int myCmp(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const;             ///< <0, =0, >0 as pp1 < = > pp2
    //    virtual int myHomogCmp(ConstRawPtr rawt1, ConstRawPtr rawt2) const;   ///< <0, =0, >0 as t1 < = > t2 assuming t1 and t2 have same multi-degree
    //    virtual int myHomogDegRevLex(ConstRawPtr rawt1, ConstRawPtr rawt2) const; ///< <0, =0, >0 as t1 < = > t2 ??? degrevlex assuming t1 and t2 have same multi-degree TO BE REMOVED

    virtual long myStdDeg(ConstRawPtr rawpp) const;                    ///< standard degree of pp
    virtual void myWDeg(degree& d, ConstRawPtr rawpp) const;                  ///< d = grading(pp)
    virtual int myCmpWDeg(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const;         ///< <0, =0, >0 as wdeg(pp1) < = > wdeg(pp2)
    virtual int myCmpWDegPartial(ConstRawPtr rawpp1, ConstRawPtr rawpp2, long n) const;      ///< <0, =0, >0 as wdeg(pp1) < = > wdeg(pp2) wrt the first n rows of weights
    virtual long myExponent(ConstRawPtr rawpp, long indet) const;           ///< degree of pp in indet
    virtual void myBigExponent(BigInt& E, ConstRawPtr rawpp, long indet) const;           ///< degree of pp in indet

    virtual void myExponents(std::vector<long>& expv, ConstRawPtr rawpp) const; ///< get exponents, SHOULD BE vector<BigInt> ????
    virtual void myBigExponents(std::vector<BigInt>& expv, ConstRawPtr rawpp) const; ///< get exponents, SHOULD BE vector<BigInt> ????
    virtual void myComputeDivMask(DivMask& dm, const DivMaskRule& DivMaskImpl, ConstRawPtr rawpp) const; ///< computes the DivMask for pp according to DivMaskImpl
    virtual void myOutputSelf(std::ostream& out) const;                    ///< print value of PPMonoid
///USE DEFAULT IMPL IN PPMONOIDBASE    virtual void myOutput(std::ostream& out, ConstRawPtr rawpp) const;   ///< NOT PURE!!
///???    virtual void myOutput(OpenMath::OutputChannel& OMOut, ConstRawPtr rawpp) const;///< NOT PURE!!

  protected: // Data members
    std::vector<PPMonoidElem> myIndetVector; ///< the indets as PPMonoidElems
    std::unique_ptr<PPMonoidElem> myOnePtr;
  };


  inline PPMonoidSparseImpl::value_t& PPMonoidSparseImpl::import(RawPtr rawpp)
  {
    return *static_cast<value_t*>(rawpp.myRawPtr());
  }

  inline const PPMonoidSparseImpl::value_t& PPMonoidSparseImpl::import(ConstRawPtr rawpp)
  {
    return *static_cast<const value_t*>(rawpp.myRawPtr());
  }


  PPMonoidSparseImpl::PPMonoidSparseImpl(const std::vector<symbol>& IndetNames, const PPOrdering& ord):
      PPMonoidBase(ord, IndetNames)
  {
    myRefCountInc();  // this is needed for exception cleanliness, in case one of the lines below throws
    myOnePtr.reset(new PPMonoidElem(PPMonoid(this)));
    myIndetVector.reserve(myNumIndets);
    vector<long> expv(myNumIndets);
    for (long i=0; i < myNumIndets; ++i)
    {
      expv[i] = 1;
      myIndetVector.push_back(PPMonoidElem(PPMonoid(this), myNew(expv)));
      expv[i] = 0;
    }
    myRefCountZero();  // ignore the "internal" references created above
  }


  const std::vector<PPMonoidElem>& PPMonoidSparseImpl::myIndets() const
  {
    return myIndetVector;
  }


  const PPMonoidElem& PPMonoidSparseImpl::myOne() const
  {
    return *myOnePtr;
  }


  PPMonoidElemRawPtr PPMonoidSparseImpl::myNew() const
  {
    std::unique_ptr<value_t> tmp;
    tmp.reset(new value_t());
    return PPMonoidElemRawPtr(tmp.release());
  }


  //  THIS IS NOT EXCEPTION CLEAN AS WRITTEN HERE: will leak if MemPool::alloc throws.
  PPMonoidElemRawPtr PPMonoidSparseImpl::myNew(PPMonoidElemConstRawPtr rawpp) const
  {
    std::unique_ptr<value_t> tmp;
    tmp.reset(new value_t(import(rawpp)));
    return PPMonoidElemRawPtr(tmp.release());
  }


  PPMonoidElemRawPtr PPMonoidSparseImpl::myNew(const std::vector<long>& expv) const
  {
    CoCoA_ASSERT(len(expv) == myNumIndets);
    std::unique_ptr<value_t> tmp;
    tmp.reset(new value_t());

    for (long i=0; i < myNumIndets; ++i)
      if (expv[i] != 0)
        tmp->push_back(IndexExp(i, expv[i]));

    return PPMonoidElemRawPtr(tmp.release());
  }

//NYI    virtual PPMonoidElemRawPtr myNew(const std::vector<BigInt>& v) const;     ///< ctor from exp vector

  void PPMonoidSparseImpl::myDelete(RawPtr rawpp) const
  {
    std::unique_ptr<value_t> tmp;
    tmp.reset(&import(rawpp));
  }


  void PPMonoidSparseImpl::mySwap(RawPtr rawpp1, RawPtr rawpp2) const
  {
    std::swap(import(rawpp1), import(rawpp2));
  }

  void PPMonoidSparseImpl::myAssignOne(RawPtr rawpp) const
  {
    import(rawpp).clear();
  }

  void PPMonoidSparseImpl::myAssign(RawPtr rawpp, ConstRawPtr rawpp1) const
  {
    import(rawpp) = import(rawpp1);
  }

  void PPMonoidSparseImpl::myAssign(RawPtr rawpp, const std::vector<long>& expv) const
  {
  }

  void PPMonoidSparseImpl::myMul(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    std::unique_ptr<value_t> ans;
    ans.reset(new value_t);
    value_t::const_iterator i = import(rawpp1).begin();
    const value_t::const_iterator endi = import(rawpp1).end();
    value_t::const_iterator j = import(rawpp2).begin();
    const value_t::const_iterator endj = import(rawpp2).end();
    while (i != endi && j != endj)
    {
      if (i->myIndetIndex < j-> myIndetIndex)
      {
        ans->push_back(*i++);
        continue;
      }
      if (i->myIndetIndex > j-> myIndetIndex)
      {
        ans->push_back(*j++);
        continue;
      }
      ans->push_back(IndexExp(i->myIndetIndex, i->myExp + j->myExp));
      ++i;
      ++j;
    }
    //    std::back_insert_iterator<value_t> k = back_inserter(*ans);
    //    copy(i, endi, k); // at most only one of these will actually copy anything
    //    copy(j, endj, k); //
    // at most only one of these two will actually do anything
    ans->insert(ans->end(), i, endi);
    ans->insert(ans->end(), j, endj);

    std::swap(import(rawpp), *ans);
  }

  void PPMonoidSparseImpl::myMulIndetPower(RawPtr rawpp, long indet, long exp) const
  {
    CoCoA_ERROR(ERR::NYI, "PPMonoidSparseImpl::myMulIndetPower");
  }

  void PPMonoidSparseImpl::myDiv(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    CoCoA_ERROR(ERR::NYI, "PPMonoidSparseImpl::myDiv");
  }

  void PPMonoidSparseImpl::myColon(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    CoCoA_ERROR(ERR::NYI, "PPMonoidSparseImpl::myColon");
  }

  void PPMonoidSparseImpl::myGcd(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    CoCoA_ERROR(ERR::NYI, "PPMonoidSparseImpl::myGcd");
  }

  void PPMonoidSparseImpl::myLcm(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    CoCoA_ERROR(ERR::NYI, "PPMonoidSparseImpl::myLcm");
  }

  void PPMonoidSparseImpl::myRadical(RawPtr rawpp, ConstRawPtr rawpp1) const
  {
    CoCoA_ERROR(ERR::NYI, "PPMonoidSparseImpl::myRadical");
  }

  void PPMonoidSparseImpl::myPowerSmallExp(RawPtr rawpp, ConstRawPtr rawpp1, long exp) const
  {
    CoCoA_ASSERT(0 <= exp);
    CoCoA_ERROR(ERR::NYI, "PPMonoidSparseImpl::myPower");
  }


  bool PPMonoidSparseImpl::myIsOne(ConstRawPtr rawpp) const
  {
    return import(rawpp).empty();
  }

  bool PPMonoidSparseImpl::myIsIndet(long& index, ConstRawPtr rawpp) const
  {
    const value_t T = import(rawpp);
    if (len(T) != 1) return false;
    index = T.front().myIndetIndex;
    return true;
  }


  bool PPMonoidSparseImpl::myIsCoprime(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    const value_t& T1 = import(rawpp1);
    const value_t& T2 = import(rawpp2);
    value_t::const_iterator i = T1.begin();
    const value_t::const_iterator endi = T1.end();
    value_t::const_iterator j = T2.begin();
    const value_t::const_iterator endj = T2.end();
    while (i != endi && j != endj)
    {
      if (i->myIndetIndex == j->myIndetIndex) return false;
      if (i->myIndetIndex < j->myIndetIndex) ++i;
      else ++j;
    }
    return true;
  }

  bool PPMonoidSparseImpl::myIsEqual(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    const value_t& T1 = import(rawpp1);
    const value_t& T2 = import(rawpp2);
    value_t::const_iterator i = T1.begin();
    const value_t::const_iterator endi = T1.end();
    value_t::const_iterator j = T2.begin();
    const value_t::const_iterator endj = T2.end();
    while (i != endi && j != endj)
    {
      if (i->myIndetIndex != j->myIndetIndex || i->myExp != j->myExp) return false;
      ++i;
      ++j;
    }
    return true;
  }

  bool PPMonoidSparseImpl::myIsDivisible(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    const value_t& T1 = import(rawpp1);
    const value_t& T2 = import(rawpp2);
    value_t::const_iterator i = T1.begin();
    const value_t::const_iterator endi = T1.end();
    value_t::const_iterator j = T2.begin();
    const value_t::const_iterator endj = T2.end();
    while (i != endi && j != endj)
    {
      if (i->myIndetIndex > j->myIndetIndex) return false;
      if (i->myIndetIndex < j->myIndetIndex) { ++i; continue; }
      if (i->myExp < j->myExp) return false;
      ++i; ++j;
    }
    return j == endj;
  }


  bool PPMonoidSparseImpl::myIsRadical(ConstRawPtr rawpp) const
  {
    const value_t& T = import(rawpp);
    value_t::const_iterator i = T.begin();
    const value_t::const_iterator endi = T.end();
    for (;i != endi; ++i)
    {
      if (i->myExp > 1) return false;
    }
    return true;
  }


  int PPMonoidSparseImpl::myCmp(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    CoCoA_ERROR(ERR::NYI, "PPMonoidSparseImpl::myCmp");
    // BUG!!! this impl INCORRECTLY ASSUMES a LEX ordering!
    const value_t& T1 = import(rawpp1);
    const value_t& T2 = import(rawpp2);
    value_t::const_iterator i = T1.begin();
    const value_t::const_iterator endi = T1.end();
    value_t::const_iterator j = T2.begin();
    const value_t::const_iterator endj = T2.end();
    while (i != endi && j != endj)
    {
      if (i->myIndetIndex != j->myIndetIndex)
      {
        if (i->myIndetIndex < j->myIndetIndex) return -1;
        return +1;
      }
      if (i->myExp == j->myExp) continue;
      ++i; ++j;
    }
    if (j != endj) return +1;
    if (i != endi) return -1;
    return 0;
  }

    //    virtual int myHomogCmp(ConstRawPtr rawt1, ConstRawPtr rawt2) const;   ///< <0, =0, >0 as t1 < = > t2 assuming t1 and t2 have same multi-degree
    //    virtual int myHomogDegRevLex(ConstRawPtr rawt1, ConstRawPtr rawt2) const; ///< <0, =0, >0 as t1 < = > t2 ??? degrevlex assuming t1 and t2 have same multi-degree TO BE REMOVED

  long PPMonoidSparseImpl::myStdDeg(ConstRawPtr rawpp) const
  {
    const value_t& T = import(rawpp);
    value_t::const_iterator i = T.begin();
    const value_t::const_iterator endi = T.end();
    long deg = 0;
    for (; i != endi; ++i)
    {
      deg += i->myExp;  //??? check for overflow????
    }
    return deg;
  }


  void PPMonoidSparseImpl::myWDeg(degree& d, ConstRawPtr rawpp) const
  {
    CoCoA_ERROR(ERR::NYI, "PPMonoidSparseImpl::myWDeg");
  }

  int PPMonoidSparseImpl::myCmpWDeg(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    CoCoA_ERROR(ERR::NYI, "PPMonoidSparseImpl::myCmpWDeg");
    return 0; // just to keep compiler quiet
  }

  int PPMonoidSparseImpl::myCmpWDegPartial(ConstRawPtr rawpp1, ConstRawPtr rawpp2, long n) const
  {
    CoCoA_ERROR(ERR::NYI, "PPMonoidSparseImpl::myCmpWDegPartial");
    return 0; // just to keep compiler quiet
  }


  long PPMonoidSparseImpl::myExponent(ConstRawPtr rawpp, long var) const
  {
    CoCoA_ASSERT(0 <= var && var <= myNumIndets);
    const value_t& T = import(rawpp);
    value_t::const_iterator i = T.begin();
    const value_t::const_iterator endi = T.end();
    for (; i != endi; ++i)
    {
      if (i->myIndetIndex == var) return i->myExp;
      if (i->myIndetIndex > var) return 0;
    }
    return 0;
  }

  void PPMonoidSparseImpl::myBigExponent(BigInt& E, ConstRawPtr rawpp, long indet) const
  {
    E = myExponent(rawpp, indet);
  }

  void PPMonoidSparseImpl::myExponents(std::vector<long>& expv, ConstRawPtr rawpp) const
  {
    CoCoA_ASSERT(len(expv) == myNumIndets);
    const value_t& T = import(rawpp);
    value_t::const_iterator i = T.begin();
    const value_t::const_iterator endi = T.end();
    for (long j=0; j < myNumIndets; ++j) expv[j] = 0;
    for (; i != endi; ++i)
    {
      expv[i->myIndetIndex] = i->myExp;
    }
  }

  void PPMonoidSparseImpl::myBigExponents(std::vector<BigInt>& expv, ConstRawPtr rawpp) const
  {
    CoCoA_ASSERT(len(expv) == myNumIndets);
    const value_t& T = import(rawpp);
    value_t::const_iterator i = T.begin();
    const value_t::const_iterator endi = T.end();
    for (long j=0; j < myNumIndets; ++j) expv[j] = 0;
    for (; i != endi; ++i)
    {
      expv[i->myIndetIndex] = i->myExp;
    }
  }

// CoCoA_ASSERT(v.size() == myNumIndets);
  void PPMonoidSparseImpl::myComputeDivMask(DivMask& dm, const DivMaskRule& DivMaskImpl, ConstRawPtr rawpp) const
  {
    CoCoA_ERROR(ERR::NYI, "PPMonoidSparseImpl::myComputeDivMask");
  }

  void PPMonoidSparseImpl::myOutputSelf(std::ostream& out) const
  {
    out << "PPMonoidSparse(" << myNumIndets << ")";
  }

/// USE GENERIC IMPL in PPMonoid.C  std::ostream& PPMonoidSparse::myOutput(std::ostream& out, ConstRawPtr rawpp) const;   ///< NOT PURE!!
//???    virtual OpenMath::OutputChannel& myOutput(OpenMath::OutputChannel& OMOut, ConstRawPtr rawpp) const;///< NOT PURE!!



  //----------------------------------------------------------------------
  // Pseudo-ctor for sparse PPMonoids

  PPMonoid NewPPMonoidSparse(const std::vector<symbol>& IndetNames, const PPOrdering& ord)
  {
  // Sanity check on the indet names given.
    const long nvars = NumIndets(ord);
  
    if (len(IndetNames) != nvars)
      CoCoA_ERROR(ERR::BadNumIndets, "NewPPMonoidSparse(IndetNames,ord)");
    if (!AreDistinct(IndetNames))
      CoCoA_ERROR(ERR::BadIndetNames, "NewPPMonoidSparse(IndetNames,ord)");
    if (!AreArityConsistent(IndetNames))
      CoCoA_ERROR(ERR::BadIndetNames, "NewPPMonoidSparse(IndetNames,ord)");
  
    return PPMonoid(new PPMonoidSparseImpl(IndetNames, ord));
  }

  PPMonoid NewPPMonoidSparse(const std::vector<symbol>& IndetNames, const PPOrderingCtor& ord)
  {
    return NewPPMonoidSparse(IndetNames, ord.myCtor(len(IndetNames)));
  }



}// end of namespace CoCoA




// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/PPMonoidSparse.C,v 1.17 2014/07/31 14:45:17 abbott Exp $
// $Log: PPMonoidSparse.C,v $
// Revision 1.17  2014/07/31 14:45:17  abbott
// Summary: Merged io.H and UtilsTemplate.H into new header VectorOperations.H
// Author: JAA
//
// Revision 1.16  2014/07/03 15:36:35  abbott
// Summary: Cleaned up impl of PPMonoids: moved myIndetSymbols & myNumIndets to base class
// Author: JAA
//
// Revision 1.15  2014/05/14 15:57:15  bigatti
// -- added "using" for clang with superpedantic flag
//
// Revision 1.14  2014/04/30 16:10:31  abbott
// Summary: Removed pointless include
// Author: JAA
//
// Revision 1.13  2012/04/02 15:03:02  bigatti
// -- using insert instead of copy+back_inserter
//
// Revision 1.12  2012/02/08 16:13:18  bigatti
// -- using  insert  instead of  back_inserter
//
// Revision 1.11  2011/08/14 15:52:17  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.10  2011/06/23 16:07:13  abbott
// Added incomplete but compiling PPMonoidSparse: first prototype,
// simple rather than fast!
//
// Revision 1.9  2010/11/30 11:18:11  bigatti
// -- renamed IndetName --> IndetSymbol
//
// Revision 1.8  2010/10/06 14:10:24  abbott
// Added increments to the ref count in ring and PPMonoid ctors to make
// them exception safe.
//
// Revision 1.7  2008/06/05 14:57:58  bigatti
// -- added: include for MemPool and DivMask
//
// Revision 1.6  2007/12/05 11:06:24  bigatti
// -- changed "size_t StdDeg/myStdDeg(f)" into "long"  (and related functions)
// -- changed "log/myLog(f, i)" into "MaxExponent/myMaxExponent(f, i)"
// -- fixed bug in "IsOne(ideal)" in SparsePolyRing.C
//
// Revision 1.5  2007/12/04 14:27:06  bigatti
// -- changed "log(pp, i)" into "exponent(pp, i)"
//
// Revision 1.4  2007/10/30 17:14:07  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.3  2007/05/31 14:54:31  bigatti
// -- now using AreDistinct and AreArityConsistent for sanity check on
//    indet names
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.10  2007/03/08 18:22:29  cocoa
// Just whitespace cleaning.
//
// Revision 1.9  2007/03/08 17:43:11  cocoa
// Swapped order of args to the NewPPMonoid pseudo ctors.
//
// Revision 1.8  2007/02/10 18:44:03  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.7  2007/01/18 14:22:44  cocoa
// -- added "raw" to RawPtr arguments
//
// Revision 1.6  2007/01/15 13:39:54  cocoa
// -- added prefix "raw" to RawPtr arguments names
//
// Revision 1.5  2006/11/24 17:04:32  cocoa
// -- reorganized includes of header files
//
// Revision 1.4  2006/10/19 13:56:24  cocoa
// Added #include<new> whenever needed (i.e. for placement new).
//
// Revision 1.3  2006/10/06 14:04:14  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.2  2006/08/07 21:23:25  cocoa
// Removed almost all publicly visible references to SmallExponent_t;
// changed to long in all PPMonoid functions and SparsePolyRing functions.
// DivMask remains to sorted out.
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.5  2006/04/28 16:33:51  cocoa
// Used SmartPtrIRC for PPOrderings.
//
// Revision 1.4  2006/02/20 22:41:20  cocoa
// All forms of the log function for power products now return SmallExponent_t
// (instead of int).  exponents now resizes the vector rather than requiring
// the user to pass in the correct size.
//
// Revision 1.3  2006/01/17 10:23:08  cocoa
// Updated DivMask; many consequential changes.
// A few other minor fixes.
//
// Revision 1.2  2005/10/24 15:22:12  cocoa
// Fixed a couple of buglets
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.1  2005/09/22 18:04:17  cocoa
// It compiles; the tests run OK.  The examples compile.
// No documentation -- the mindless eurocrats have rendered
// me mindless too.
//
