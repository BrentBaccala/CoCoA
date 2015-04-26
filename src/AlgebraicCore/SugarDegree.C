//   Copyright (c)  2009  Anna M. Bigatti

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


#include "CoCoA/SugarDegree.H"
#include "CoCoA/IntOperations.H"
#include "CoCoA/PPMonoid.H"
#include "CoCoA/ReductionCog.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/TmpGPoly.H"
#include "CoCoA/assert.H"

#include <algorithm>
using std::max; // for std sugar
#include <iostream>
//using std::ostream;
//#include <memory> // already included by SugarDegree.H
using std::auto_ptr;


namespace CoCoA
{
  std::ostream& operator<<(std::ostream& out, const SugarDegree& s)
  {
    if (IsInitialized(s)) return s->myOutput(out);
    return out << "uninitialized";
  }


  SugarDegree::SugarDegree(const SugarDegree& rhs)
  {
    if (IsInitialized(rhs)) myPtr.reset((rhs.myPtr)->myClone());
  }


  SugarDegree& SugarDegree::operator=(const SugarDegree& rhs)
  {
    if (this == &rhs) return *this;
    if (IsInitialized(rhs)) myPtr.reset((rhs.myPtr)->myClone());
    else                    myPtr.reset(0);
    return *this;
  }
  

  void SugarDegreeBase::myUpdate(ReductionCog F, const GPoly& g)
  { myUpdate(ActiveLPP(F)/LPPForOrd(g), g); }
  

  namespace SugarTypes
  {
    //---  StdDegBase and WDegBase abstract classes -------------

    class StdDegBase: public SugarDegreeBase
    {
    public:
      StdDegBase(long s): SugarDegreeBase() { myValue = s; }
      virtual const degree& myWSugar() const;  ///< returns error
      virtual long myStdSugar() const { return myValue; }
      virtual int myCmp(const SugarDegreeBase& s) const;  // this <=> s ? <0,=0,>0

      virtual std::ostream& myOutput(std::ostream& out) const;
      virtual void myMul(ConstRefPPMonoidElem pp);
      using SugarDegreeBase::myUpdate; // disables warnings of overloading
      virtual void myUpdate(ConstRefPPMonoidElem CofactorPP, const GPoly& g);
      virtual long myDeg(ConstRefPPMonoidElem pp) const = 0;
      virtual long myDeg0(ConstRefRingElem f) const; // based on myDeg; also defined on 0
    protected:
      long myValue;
    };


    class WDegBase: public SugarDegreeBase
    {
    public:
      WDegBase(const degree& s): SugarDegreeBase(), myValue(s) {}
      virtual const degree& myWSugar() const { return myValue; }
      virtual long myStdSugar() const;   ///< returns error
      virtual int myCmp(const SugarDegreeBase& s) const;  // this <=> s ? <0,=0,>0

      virtual std::ostream& myOutput(std::ostream& out) const;
      virtual void myMul(ConstRefPPMonoidElem pp);
      using SugarDegreeBase::myUpdate; // disables warnings of overloading
      virtual void myUpdate(ConstRefPPMonoidElem CofactorPP, const GPoly& g);
      virtual degree myDeg(ConstRefPPMonoidElem pp) const = 0;
      degree myDeg0(ConstRefRingElem f) const; // based on myDeg; also defined on 0
    protected:
      degree myValue;
    };


    //---  StdDegBase concrete classes -------------

    class StdDegPlain: public StdDegBase
    {
    public:
      StdDegPlain(long s): StdDegBase(s) {}
      StdDegPlain(ConstRefRingElem f): StdDegBase(0) { myValue=myDeg0(f); }
      StdDegPlain* myClone() const { return new StdDegPlain(myValue); }
      virtual long myDeg(ConstRefPPMonoidElem pp) const;
    };


    class StdDegNoIdx: public StdDegBase
    {
    public:
      StdDegNoIdx(long s, long NoIdx): StdDegBase(s) { myNoIdxValue=NoIdx;}
      StdDegNoIdx(ConstRefRingElem f, long NoIdx): StdDegBase(0) { myNoIdxValue=NoIdx; myValue=myDeg0(f); }
      StdDegNoIdx* myClone() const { return new StdDegNoIdx(myValue, myNoIdxValue); }
      virtual long myDeg(ConstRefPPMonoidElem pp) const;
    private:
      long myNoIdxValue;
    };


    class StdDegNoIdxSat: public StdDegBase
    {
    public:
      StdDegNoIdxSat(long s, long NoIdx): StdDegBase(s) { myNoIdxValue=NoIdx;}
      StdDegNoIdxSat(ConstRefRingElem f, long NoIdx): StdDegBase(0) { myNoIdxValue=NoIdx; myValue=myDeg0(f); }
      StdDegNoIdxSat* myClone() const { return new StdDegNoIdxSat(myValue, myNoIdxValue); }
      virtual long myDeg(ConstRefPPMonoidElem pp) const;
    private:
      long myNoIdxValue;
    };


    class StdDegSat: public StdDegBase
    {
    public:
      StdDegSat(long s): StdDegBase(s) {}
      StdDegSat(ConstRefRingElem f): StdDegBase(0) { myValue=myDeg0(f); }
      StdDegSat* myClone() const { return new StdDegSat(myValue); }
      virtual long myDeg(ConstRefPPMonoidElem pp) const;
    };


    //---  WDegBase concrete classes -------------

    class WDegPlain: public WDegBase
    {
    public:
      WDegPlain(const degree& s): WDegBase(s) {}
      WDegPlain(ConstRefRingElem f);
      WDegPlain* myClone() const { return new WDegPlain(myValue); }
      virtual degree myDeg(ConstRefPPMonoidElem pp) const;
    };


    class WDeg1CompTmp: public WDegBase
    {
    public:
      WDeg1CompTmp(const degree& s): WDegBase(s) {}
      WDeg1CompTmp(ConstRefRingElem f);
      WDeg1CompTmp* myClone() const { return new WDeg1CompTmp(myValue); }
      virtual degree myDeg(ConstRefPPMonoidElem pp) const;
      using SugarDegreeBase::myUpdate; // disables warnings of overloading
      virtual void myUpdate(ConstRefPPMonoidElem CofactorPP, const GPoly& g); ///< no default
    };


    class WDegConst: public WDegBase
    {
    public:
      WDegConst(const degree& s): WDegBase(s) {}
      WDegConst(ConstRefRingElem f);
      WDegConst* myClone() const { return new WDegConst(myValue); }
      virtual degree myDeg(ConstRefPPMonoidElem pp) const;
      using SugarDegreeBase::myUpdate; // disables warnings of overloading
      virtual void myUpdate(ConstRefPPMonoidElem /*CofactorPP*/, const GPoly& /*g*/) {}; ///< no default
    };


    class WDegSat: public WDegBase ///< wip
    {
    public:
      WDegSat(const degree& s): WDegBase(s) {}
      WDegSat(ConstRefRingElem f);
      WDegSat* myClone() const { return new WDegSat(myValue); }
      virtual degree myDeg(ConstRefPPMonoidElem pp) const;
      using SugarDegreeBase::myUpdate; // disables warnings of overloading
      virtual void myUpdate(ConstRefPPMonoidElem CofactorPP, const GPoly& g);
    };



    //--  StdDegBase  ----------------------------------

    const degree& StdDegBase::myWSugar() const
    {
      CoCoA_ERROR("Wrong type for this sugar", "SugarTypes::StdDegBase");
      static degree NeverUsed(0); // just to keep the compiler quiet
      return NeverUsed;
    }

    std::ostream& StdDegBase::myOutput(std::ostream& out) const
    { return out << myValue; }

    void StdDegBase::myMul(ConstRefPPMonoidElem pp)
    { myValue += myDeg(pp); }

    void StdDegBase::myUpdate(ConstRefPPMonoidElem CofactorPP, const GPoly& g)
    { myValue = max(myValue, myDeg(CofactorPP) + sugar(g)->myStdSugar()); }


    int StdDegBase::myCmp(const SugarDegreeBase& s) const
    {
      CoCoA_ASSERT((dynamic_cast<const StdDegBase*>(&s)) != 0);
      const StdDegBase& ss = static_cast<const StdDegBase&>(s);
      if (myValue != ss.myValue)
      {
        if (myValue > ss.myValue) return 1; else return -1;
      }
      return 0;
    }

    long StdDegBase::myDeg0(ConstRefRingElem f) const
    {
      if (IsZero(f)) return 0;
      return myDeg(LPP(f));  // ANNA: wrong! the right one follows
      // after finishing the testing phase use the correct definition
      long degf = 0;
      for (SparsePolyIter it = BeginIter(f); !IsEnded(it); ++it)
        degf = max(degf, myDeg(PP(it)));
      return degf;
    }
    
    //--  WDegBase  ----------------------------------

    long WDegBase::myStdSugar() const
    {
      CoCoA_ERROR("Wrong type for this sugar", "SugarTypes::WDegBase");
      return 0;  // just to keep the compiler quiet
    }

    std::ostream& WDegBase::myOutput(std::ostream& out) const
    { return out << myValue; }

    void WDegBase::myMul(ConstRefPPMonoidElem pp)
    { myValue += myDeg(pp); }

    void WDegBase::myUpdate(ConstRefPPMonoidElem CofactorPP, const GPoly& g)
    { myValue = top(myValue, myDeg(CofactorPP) + sugar(g)->myWSugar()); }

    int WDegBase::myCmp(const SugarDegreeBase& s) const
    {
      CoCoA_ASSERT(dynamic_cast<const WDegBase*>(&s) != 0);
      const WDegBase& ws = static_cast<const WDegBase&>(s);
      return FastCmp(myValue, ws.myValue); // cmp = compatibility test + FastCmp
    }

    degree WDegBase::myDeg0(ConstRefRingElem f) const
    {
      if (IsZero(f)) return myDeg(one(PPM(owner(f))));
      return myDeg(LPP(f));
      // ANNA: BUG BUG BUG return above is wrong!!!! the right one follows
      SparsePolyIter it = BeginIter(f);
      degree degf = myDeg(PP(it));
      for (++it; !IsEnded(it); ++it)
        degf = top(degf, myDeg(PP(it)));
      return degf;
    }

    //--  end abstract classes  ----------------------------------

    //--  StdDegPlain  ----------------------------------

    long StdDegPlain::myDeg(ConstRefPPMonoidElem pp) const
    { return StdDeg(pp); }

    //--  StdDegSat  ----------------------------------

    long StdDegSat::myDeg(ConstRefPPMonoidElem pp) const
    { return StdDeg(pp) - exponent(pp, NumIndets(owner(pp))-1); }

    //--  StdDegNoIdx  ----------------------------------

    long StdDegNoIdx::myDeg(ConstRefPPMonoidElem pp) const
    { return StdDeg(pp) - exponent(pp, myNoIdxValue); }

    //--  StdDegNoIdxSat  ----------------------------------

    long StdDegNoIdxSat::myDeg(ConstRefPPMonoidElem pp) const
    {
      return StdDeg(pp)
        - exponent(pp, myNoIdxValue)
        - exponent(pp, NumIndets(owner(pp))-1);
    }

    //--  WDegPlain  ----------------------------------

    WDegPlain::WDegPlain(ConstRefRingElem f): WDegBase(degree(GradingDim(owner(f)))) // initialize
    { myValue = myDeg0(f); }

    degree WDegPlain::myDeg(ConstRefPPMonoidElem pp) const
    { return wdeg(pp); }

    //--  WDegConst  ----------------------------------

    WDegConst::WDegConst(ConstRefRingElem f): WDegBase(degree(GradingDim(owner(f))))
    { if (!IsZero(f)) myValue=wdeg(LPP(f)); }

    degree WDegConst::myDeg(ConstRefPPMonoidElem pp) const
    { return wdeg(pp); }

    // inline myUpdate(pp, g) does *nothing* (unlike default implementation)

    //--  WDeg1CompTmp  ----------------------------------

    WDeg1CompTmp::WDeg1CompTmp(ConstRefRingElem f): WDegBase(degree(GradingDim(owner(f))))
    { myValue = myDeg0(f); }

    degree WDeg1CompTmp::myDeg(ConstRefPPMonoidElem pp) const
    { return wdeg(pp); }

    void WDeg1CompTmp::myUpdate(ConstRefPPMonoidElem CofactorPP, const GPoly& g)
    { // very temporary!!  only for checking backward compatibility
      // moreover it uses max instead of top (a>b ? a : b)
      // could it be useful with only one component?
      if (g.myGRingInfo().IsMyGradingPosPlus())
        CoCoA_ERROR("TMP: use StdDeg for PosTO", "WDeg1CompTmp::myUpdate");
      myValue = max(myValue, wdeg(CofactorPP) + sugar(g)->myWSugar());
    }

    //--  WDegSat  ----------------------------------

    WDegSat::WDegSat(ConstRefRingElem f): WDegBase(degree(GradingDim(owner(f))))
    { myValue = myDeg0(f); }


    degree WDegSat::myDeg(ConstRefPPMonoidElem pp) const
    {
      degree tmp = wdeg(pp);
      // bi-grading trick
      // tmp.mySetComponent(0, tmp[1]);//WARNING: this works only for GrDim=1
      // TEMPORARY
      const BigInt d = tmp[0] - exponent(pp, NumIndets(owner(pp))-1); // only for GrDim=1
      SetComponent(tmp, 0, d);
      return tmp;
    }


    void WDegSat::myUpdate(ConstRefPPMonoidElem CofactorPP, const GPoly& g)
    {
      if (g.myGRingInfo().IsMyGradingPosPlus())
        CoCoA_ERROR("TMP: use StdDeg for PosTO", "WDeg1CompTmp::myUpdate");
      myValue = max(myValue, myDeg(CofactorPP) + sugar(g)->myWSugar());
    }

    //-------

  }  // namespace SugarTypes

  //----------------------------------------------------------------------
  // Here are the pseudo-constructors:

  SugarDegree NewStdSugar(ConstRefRingElem f)
  { return SugarDegree(new SugarTypes::StdDegPlain(f)); }

  SugarDegree NewStdSugarNoIdx(ConstRefRingElem f, long NoIdx)
  { return SugarDegree(new SugarTypes::StdDegNoIdx(f, NoIdx)); }

  SugarDegree NewStdSugarSat(ConstRefRingElem f)
  { return SugarDegree(new SugarTypes::StdDegSat(f)); }

  SugarDegree NewStdSugarNoIdxSat(ConstRefRingElem f, long NoIdx)
  { return SugarDegree(new SugarTypes::StdDegNoIdxSat(f, NoIdx)); }


  SugarDegree NewWSugar(ConstRefRingElem f)
  { return SugarDegree(new SugarTypes::WDegPlain(f)); }

  SugarDegree NewWDeg1CompTmp(ConstRefRingElem f)
  { return SugarDegree(new SugarTypes::WDeg1CompTmp(f)); }

  SugarDegree NewWSugarConst(ConstRefRingElem f)
  { return SugarDegree(new SugarTypes::WDegConst(f)); }

  SugarDegree NewWSugarSat(ConstRefRingElem f)
  { return SugarDegree(new SugarTypes::WDegSat(f)); }



} // namespace CoCoA



//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/SugarDegree.C,v 1.18 2014/07/07 12:45:30 abbott Exp $
// $Log: SugarDegree.C,v $
// Revision 1.18  2014/07/07 12:45:30  abbott
// Summary: Removed AsPolyRing
// Author: JAA
//
// Revision 1.17  2014/05/09 09:29:03  bigatti
// added  using SugarDegreeBase::myUpdate; // disables warnings of overloading
//
// Revision 1.16  2012/05/28 09:18:20  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.15  2011/08/14 15:52:16  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.14  2010/05/21 16:11:37  bigatti
// -- removed CheckSugar
//
// Revision 1.13  2010/05/14 09:53:09  bigatti
// -- removed empty ctor for SugarDegree
// -- added marker for SugarDegree(uninitialized)
// -- SugarDegree for GBasis input is initialized by myPrepareGBasis
//
// Revision 1.12  2009/10/27 17:15:55  bigatti
// -- fixed: using sugar(g)->myWSugar() insted of wsugar(g)
//
// Revision 1.11  2009/09/25 12:36:47  bigatti
// -- added CheckSugar (temporary)
//
// Revision 1.10  2009/07/24 15:13:50  bigatti
// -- some fixes
//
// Revision 1.9  2009/07/20 14:28:44  bigatti
// -- modified interface for pseudo constructors (now with RingElem argument)
//
// Revision 1.8  2009/04/24 16:37:41  bigatti
// -- added myDeg and pseudo-ctor
// -- improved inheritance
//
// Revision 1.7  2009/03/20 14:01:34  bigatti
// -- minor changes (names)
//
// Revision 1.6  2009/03/18 16:37:13  bigatti
// -- "almost-final" cleanup, tested against some current sugar definitions
//
// Revision 1.5  2009/03/18 15:12:02  abbott
// Minor cosmetic changes.
//
// Revision 1.4  2009/03/16 16:39:48  bigatti
// -- now use auto_ptr instead of SmartPtrIRC
// -- added copy ctor, assignment, cmp
//
// Revision 1.3  2009/03/16 07:27:20  bigatti
// -- added necessary "const"
// -- added WDeg1CompTmp (temporary, for testing)
//
// Revision 1.2  2009/02/20 15:45:25  bigatti
// -- sugar.H --> SugarDegree.H
//
// Revision 1.1  2009/02/20 13:27:15  bigatti
// -- renamed from "sugar.[CH]"
//
// Revision 1.3  2009/02/20 11:01:15  bigatti
// -- added NewHomogWSugar for graded=homoeneous case (just constant)
//
// Revision 1.2  2009/02/20 09:53:27  bigatti
// -- changed name: sweetener --> SugarDegree
// -- introduced abstract classes SugarTypes::StdDegBase, SugarTypes::WDegBase
//
// Revision 1.1  2009/02/09 13:57:05  bigatti
// -- first import
//
