//   Copyright (c)  2005  Massimo Caboara

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


#include "CoCoA/TmpGPair.H"
#include "CoCoA/IntOperations.H"

#include <iostream>
using std::ostream;
using std::endl;
#include <algorithm>
using std::find;
using std::find_if;
using std::max;
#include <functional>
using std::less;
using std::bind1st;
using std::bind2nd;
using std::binary_function;

namespace CoCoA
{

  // -- ANNA: temporary
  SugarDegree NewSugar(const GPair& gp)
  {
    PPMonoidElem cofactor1 = LPP(OrdPoly(gp))/LPPForOrd(gp.myFirstGPoly());
    //PP(myLCMwMask)/LPPForDiv(myFirstGPoly());
    PPMonoidElem cofactor2 = LPP(OrdPoly(gp))/LPPForOrd(gp.mySecondGPoly());
    SugarDegree s(sugar(gp.myFirstGPoly()));
    s->myMul(cofactor1);
    s->myUpdate(cofactor2, gp.mySecondGPoly());
    return s;
  }



// GPair ////////////////////////////////////////////////////////////////////
 // Special pair: it represents an input polynomial 
  GPair::GPair(GPoly& the_p, unsigned int i):
    myLCMwMask(the_p.myGRingInfo().myPPM(), the_p.myGRingInfo().myDivMaskRule()),
    myOrdPoly(monomial(owner(the_p), 1, LPPForOrd(the_p))),
    myWDeg(wdeg(myOrdPoly)),
    mySugar(sugar(the_p))
  {
    myFirstGPolyPtr=&the_p;
    mySecondGPolyPtr=0;	
    myLCMwMask = LPPForDivwMask(the_p);
    IamCoprimeFlag = false;
    //ANNA    std::cout << "GPair: i = " << i << " " << Age(the_p) << " = Age" << std::endl;
    myFirstIndexValue=i;
    mySecondIndexValue=-1;
    myComponent=Component(the_p);
    //    myStdDeg=0;
 }//GPair::GPair


  GPair::GPair(GPoly& p, GPoly& q,  unsigned int i, unsigned int j):
    myLCMwMask(p.myGRingInfo().myPPM(), p.myGRingInfo().myDivMaskRule()),
    myOrdPoly(one(owner(p))),// Fake, filled by myComplete
    myWDeg(GradingDim(owner(p))),// Fake, filled by myComplete
    mySugar(uninitialized)// Fake, filled by myComplete
  {
    (void)(i); (void)(j); // to avoid compiler warning about unused parameter
    myFirstGPolyPtr=&p;
    mySecondGPolyPtr=&q;
    myLCMwMask = lcm(LPPForDiv(p), LPPForDiv(q));
    //    IamCoprimeFlag = (PP(myLCMwMask) == LPPForDiv(p)*LPPForDiv(q) );// MOD Aggiustare
    IamCoprimeFlag = IsCoprime(LPPForDiv(p), LPPForDiv(q));
    //ANNA    std::cout << "GPair: i = " << i << " " << Age(p) << " = Age  &&  " << 
    //ANNA    "j = " << j << " " << Age(q) << " = Age" << std::endl;
//     myFirstIndexValue=i;
//     mySecondIndexValue=j;	
    myFirstIndexValue = 0; // testing 0
    mySecondIndexValue = 0;	 // testing 0
    myComponent=Component(p);//had to be : p.Component()==q.Component()
    //    myStdDeg=0;// Fake, filled by myComplete
  }

  // ANNA: not called
//   GPair::GPair(const GPoly* p,const GPoly* q,  unsigned int i, unsigned int j):
//     myLCMwMask(p->myGRingInfo().myPPM(), p->myGRingInfo().myDivMaskRule()),
//     myOrdPoly(one(owner(*p))),// Fake, filled by myComplete
//     myWDeg(GradingDim(owner(*p))),// Fake, filled by myComplete
//     mySugar(uninitialized)// Fake, filled by myComplete
//   {
//     myFirstGPolyPtr=p;
//     mySecondGPolyPtr=q;
//     myLCMwMask = lcm(LPPForDiv(*p), LPPForDiv(*q));
//     //    IamCoprimeFlag = (PP(myLCMwMask) == LPP(*p)*LPP(*q) );// MOD Aggiustare
//     IamCoprimeFlag = IsCoprime(LPPForDiv(*p), LPPForDiv(*q));
//     std::cout << "GPair*: i = " << i << " " << Age(*p) << " = Age" << std::endl;
//     std::cout << "GPair*: j = " << j << " " << Age(*q) << " = Age" << std::endl;
//     myFirstIndexValue=i;
//     mySecondIndexValue=j;	
//     myComponent=Component(*p);//had to be : p->Component()==q->Component()
//     //    myStdDeg=0;// Fake, filled by myComplete
//  }


  void GPair::myComplete()
  {
    const SparsePolyRing P = owner(myOrdPoly);
    std::vector<long> expv;
    exponents(expv, PP(myLCMwMask));
    myOrdPoly = monomial(P, 1, expv);
    myWDeg = wdeg(PP(myLCMwMask));
    mySugar = NewSugar(*this);
    //  myStdDeg=GRI.TmpStdDeg(PP(myLCMwMask));  // ANNA is this used??
  }//myComplete


bool GPair::operator==(const GPair& P)const
{
//   return (this->myFirstIndexValue==P.myFirstIndexValue
//           &&
//           this->mySecondIndexValue==P.mySecondIndexValue);
// check on pointers
   return (myFirstGPolyPtr==P.myFirstGPolyPtr
           &&
           mySecondGPolyPtr==P.mySecondGPolyPtr);
}//operator==


ostream& operator<<(ostream& out, const GPair& P)
{
  out<<"<"<<P.myFirstIndexValue<<", ";
  if (P.mySecondIndexValue==0) out<<0;else out<<P.mySecondIndexValue;
  out<<",Is InputPoly="<<P.IsInputPoly();
  if (P.myComponent!=0) out<<", Comp="<<P.myComponent;out<<",";
  out<<PP(P.myLCMwMask)<<",COP="<<P.IamCoprimeFlag;
  out<<",Deg="<< P.myWDeg;
  out<<",Sugar="<< P.mySugar;
//out<<" FirstPoly ="<< *(P.myFirstGPolyPtr);
//out<<" FirstPoly Deg="<< wdeg(*(P.myFirstGPolyPtr));
//out<<"OP Type "<<GetGRingInfo(P).myInputAndGrading();
  return out;
}//operator<<


void Ordered_Insert(GPairList& L,GPair P)
{
  GPairList::iterator it=find_if(L.begin(),L.end(),bind1st(less<GPair>(),P));
  L.insert(it,P);
}//Ordered_Insert

void RemoveFromGPList(GPairList& L,GPair& P)
{
  GPairList::iterator it=find(L.begin(),L.end(),P);
  if (it!=L.end())
    L.erase(it);
}//RemoveFromGPList


const ring& CoeffRing(const GPair& P)
{
  return CoeffRing(*(P.myFirstGPolyPtr));
}//CoeffRing

const SparsePolyRing& owner(const GPair& P)
{
  return owner(*(P.myFirstGPolyPtr));
}//owner


// ***************************************************  ModuleGPairsList

// *********** These 3 function objects are stricly for use in the procdures below
template <class T1, class T2>
struct GPairDividesLCM:public binary_function<T1,T2,bool>
{
  bool operator()(const GPair& P1, const GPair& P2) const
    {
      return IsDivisibleFast(LCMwMask(P2), LCMwMask(P1)) && GPairComponent(P1)==GPairComponent(P2);
    }
};

template <class T1, class T2>
struct GPairDividesLCMProperly:public binary_function<T1,T2,bool>
{
  bool operator()(const GPair& P1, const GPair& P2) const
    {
      return IsDivisibleFast(LCMwMask(P2), LCMwMask(P1))
        && (LCMwMask(P1)!=LCMwMask(P2)) &&
        GPairComponent(P1)==GPairComponent(P2);
    }
};

template <class T1, class T2>
struct GPairEqualLCM:public binary_function<T1,T2,bool>
{
  bool operator()(const GPair& P1, const GPair& P2) const
    {
      return LCMwMask(P1)==LCMwMask(P2) && GPairComponent(P1)==GPairComponent(P2);
    }
};
// ***********

ModuleGPairList::ModuleGPairList()
{
  myMGPList.resize(10000);
}

void ModuleGPairList::Insert(GPair& P)
{
  Ordered_Insert(myMGPList[P.mySecondIndex()],P);
  //make a copy of en element that will be destroyed - think about using splice
}

GPairList::iterator ModuleGPairList::FindDivides(GPair& P,bool& found)
{
  GPairList::iterator it;
  it=find_if(myMGPList[P.mySecondIndex()].begin(),
             myMGPList[P.mySecondIndex()].end(),
             bind2nd(GPairDividesLCM<GPair,GPair>(),P));
  found=(it!=myMGPList[P.mySecondIndex()].end());
  return it;
}//FindDivides



bool ModuleGPairList::IsIn(const GPair& P)
{
  return find(myMGPList[P.mySecondIndex()].begin(),
              myMGPList[P.mySecondIndex()].end(),
              P)
    !=
    myMGPList[P.mySecondIndex()].end();
}//IsIn


GPairList::iterator ModuleGPairList::FindSameLCMAndSecondInd(GPair& P,bool& found)
{
  GPairList::iterator it;
  it=find_if(myMGPList[P.mySecondIndex()].begin(),
             myMGPList[P.mySecondIndex()].end(),
             bind2nd(GPairEqualLCM<GPair,GPair>(),P));
  found=(it!=myMGPList[P.mySecondIndex()].end());
  return it;
}//FindSameLCMAndSecondInd


long ModuleGPairList::size() const
{
  long S=0;
  for (int i=0; i!=len(myMGPList);++i)
    S += len(myMGPList[i]);
  return S;
}



ostream& operator<<(ostream& out,const ModuleGPairList& MGPL)
{
  MGPairList::const_iterator it;
  //out<<endl<<"SIZE"<<len(MGPL.myMGPList)<<endl;
  for (it=MGPL.myMGPList.begin();it!=MGPL.myMGPList.end();++it)
  {
    //out<<".";
    if (!it->empty())
    {
      GPairList::const_iterator it1=it->begin();
      out<<endl<<"Index="<<it1->mySecondIndex()<<endl;
      for (;it1!=it->end();++it1) {out<<*it1;};
    }
  }
  return out;
}//operator<<


/********** End Module GPLists ********************************************/

// ex inline //

bool GPair::BCriterion_OK(const PPWithMask& NewPP)const
{
  if (IsInputPoly()) return true;	
  return  !(IsDivisibleFast(myLCMwMask, NewPP)
            &&(lcm(LPPForDiv(*mySecondGPolyPtr),PP(NewPP)) != PP(myLCMwMask))
            &&(lcm(LPPForDiv(*myFirstGPolyPtr),PP(NewPP)) != PP(myLCMwMask)));
}


///*
//  You have the choice of DegRevLex (Default) and the ring ordering
//  GPair ordering is SparsePolyRing (owner(myOrdPoly)) ordering
bool GPair::operator<(const GPair& the_gp)const
{
  int CMP;
  if ( (CMP=cmp(mySugar, the_gp.mySugar))!=0 )   return CMP<0;
  //  if (! IsConstWSugar(mySugar) )  // Anna: avoid useless computation in homog case
  if ( (CMP=FastCmp(myWDeg, the_gp.myWDeg))!=0 ) return CMP<0;
  if (IsInputPoly()!=the_gp.IsInputPoly())
    return the_gp.IsInputPoly();
  return SparsePolyRingPtr(owner(myOrdPoly))->myCmpLPP(raw(myOrdPoly), raw(the_gp.myOrdPoly))<0;
}

}// end namespace cocoa

// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/TmpGPair.C,v 1.15 2014/07/07 13:03:19 abbott Exp $
// $Log: TmpGPair.C,v $
// Revision 1.15  2014/07/07 13:03:19  abbott
// Summary: Removed AsSparsePolyRing
// Author: JAA
//
// Revision 1.14  2014/06/17 10:14:24  abbott
// Summary: Added (void)(var_name) to avoid compiler warning about unused param
// Author: JAA
//
// Revision 1.13  2014/04/30 16:15:13  abbott
// Summary: Replaced X.size() by len(X)
// Author: JAA
//
// Revision 1.12  2014/01/28 11:01:40  bigatti
// -- changed names myFirstGPolyValue --> myFirstGPolyPtr
// -- equality test based on ptr instead of indices
//
// Revision 1.11  2013/10/28 13:12:40  bigatti
// -- IsSpecial --> IsInputPoly
//
// Revision 1.10  2012/05/28 09:18:20  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.9  2010/05/14 09:53:09  bigatti
// -- removed empty ctor for SugarDegree
// -- added marker for SugarDegree(uninitialized)
// -- SugarDegree for GBasis input is initialized by myPrepareGBasis
//
// Revision 1.8  2009/11/20 16:00:52  bigatti
// -- removed unused myStdDeg
//
// Revision 1.7  2009/10/27 17:15:14  bigatti
// -- fixed: using sugar(g)->myWSugar() insted of wsugar(g)
//
// Revision 1.6  2009/04/27 13:23:38  bigatti
// -- changed (faster check in operator<)
//
// Revision 1.5  2009/04/27 12:28:22  bigatti
// -- added   case SaturatingAlgNoDRL
//
// Revision 1.4  2008/09/19 13:33:42  bigatti
// -- added: Sat algorithm (M.Caboara)
//
// Revision 1.3  2008/09/16 15:03:42  bigatti
// -- added LPPForDiv
// -- changed LPP into LPPForOrd
//
// Revision 1.2  2007/10/30 17:14:07  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1  2007/03/09 18:56:56  bigatti
// -- added Tmp prefix to Groebner related files
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.18  2007/03/08 18:22:29  cocoa
// Just whitespace cleaning.
//
// Revision 1.17  2007/03/08 16:38:43  bigatti
// -- fix: workaround for bug with max(ZZ:rtn, ZZ:rtn)
//
// Revision 1.16  2007/03/07 17:04:31  cocoa
// -- several changes by M.Caboara: more operations on ideals,
//    exception cleaner, coding conventions, WSugar, dynamic
//
// Revision 1.15  2007/02/10 18:44:03  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.14  2006/12/21 13:48:33  cocoa
// Made all increment/decrement calls prefix (except where the must be postfix).
//
// Revision 1.13  2006/12/04 13:55:54  cocoa
// -- added: sugar for GradingDim > 0  (called wsugar)
//
// Revision 1.12  2006/11/24 17:01:42  cocoa
// -- reorganized includes of header files
//
// Revision 1.11  2006/11/20 14:57:17  cocoa
// -- added: (standard) sugar for modules
// -- fixed: non-homogeneous sysygies
// -- minor fixes     [M.Caboara]
//
// Revision 1.10  2006/10/11 13:33:06  cocoa
// -- rearranged code for sugar in reduce.C
// -- activated sugar in GPair.C
// -- removed operator<< for GPairList (template in io)
//
// Revision 1.9  2006/10/06 16:36:21  cocoa
// -- (minor)
//
// Revision 1.8  2006/10/06 15:14:28  cocoa
// -- minor change: removed useless computation of StdDeg
//
// Revision 1.7  2006/10/06 14:04:15  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.6  2006/08/17 09:24:17  cocoa
// -- added: sugar
// -- changed: coding conventions
//
// Revision 1.5  2006/06/20 17:19:45  cocoa
// -- moved  GPair::operator<  into .C file
//
// Revision 1.4  2006/06/13 16:42:32  cocoa
// -- alternative way to compute coprimality in GPair constructors
//
// Revision 1.3  2006/06/09 16:21:47  cocoa
// -- small adjustments on StdDeg computations
//
// Revision 1.2  2006/06/09 16:06:04  cocoa
// -- myStdDeg computed only if GradingDim==0
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.9  2006/05/02 14:37:12  cocoa
// -- only changes to logging info (by M.Abshoff)
//
// Revision 1.8  2006/04/21 16:12:51  cocoa
// -- just commented code
//
// Revision 1.7  2006/04/11 14:54:21  cocoa
// -- added: myStdDeg
//
// Revision 1.6  2006/04/10 17:02:47  cocoa
// -- changed: BCriterion_OK now uses PPWithMask instead of ConstRefPPMonoidElem
//
// Revision 1.5  2006/03/02 13:44:39  cocoa
// -- just comments
//
// Revision 1.4  2006/01/20 15:43:30  cocoa
// -- fixed: use of RefPPMonoidElem and ConstRefPPMonoidElem
//
// Revision 1.3  2006/01/17 15:44:56  cocoa
// -- chamges by Max for operations with modules
//
// Revision 1.2  2005/12/31 12:22:18  cocoa
// Several minor tweaks to silence the Microsoft compiler:
//  - added some missing #includes and using directives
//  - moved some function defns into the right namespace
//  - etc.
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.2  2005/07/01 16:08:15  cocoa
// Friday check-in.  Major change to structure under PolyRing:
// now SparsePolyRing and DUPolyRing are separated (in preparation
// for implementing iterators).
//
// A number of other relatively minor changes had to be chased through
// (e.g. IndetPower).
//
// Revision 1.1.1.1  2005/05/03 15:47:31  cocoa
// Imported files
//
// Revision 1.2  2005/04/19 14:06:04  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.10  2004/11/19 15:44:27  cocoa
// Changed names of "casting" functions which convert a ring into
// one with a more special structure (e.g. FractionField).  These
// functions now have names starting with "As".  There were several
// consequential changes.
//
// Revision 1.9  2004/11/18 18:33:41  cocoa
// Now every ring know its own "one" element (as well as "zero").
// Several consequential changes.
//
// Revision 1.8  2004/11/09 16:30:51  cocoa
// Removed references to cout.
//
// Revision 1.7  2004/10/29 15:51:37  cocoa
// -- changed myLCM into myLCMwMask (PPMonoidElem --> PPWithMask)
// -- function IsDivisible had wrong semantics --> swapped arguments everywhere
//
// Revision 1.6  2004/06/16 16:13:41  cocoa
// Improved I/O facilities with knock-on changes
//
// Revision 1.5  2004/05/27 16:27:55  cocoa
// -- removed ";" at the end of function bodies (g++ gives error on them)
//
// Revision 1.4  2003/11/14 14:09:53  cocoa
// Max - clean up
//
// Revision 1.3  2003/10/09 12:48:17  cocoa
// New coding convention for rings.
//
// Revision 1.2  2003/10/01 10:35:32  cocoa
// - applied "my" coding convention to PPMonoid and PPOrdering
//
// Revision 1.1.1.1  2003/09/24 12:55:43  cocoa
// Imported files
//
// Revision 1.11  2003/09/22 17:23:59  bigatti
// - new field myOrdPoly to order the GPairs
//
// Revision 1.10  2003/06/23 17:10:41  abbott
// Minor cleaning prior to public release.
// Improved the include directives,
//
// Revision 1.9  2003/05/29 16:14:12  bigatti
// - moved Component(..) to GPair.H
//
// Revision 1.8  2003/05/28 14:21:29  bigatti
// - new code for modules
//
// Revision 1.7  2003/05/14 16:40:13  bigatti
// - myDeg is now of type degree
// - new ring/PPMonoid syntax
// - divides --> IsDivisible
//
// Revision 1.6  2002/11/15 17:15:00  bigatti
// - changed GPair::GPair to avoid copy of lcm
//
// Revision 1.5  2002/09/19 17:20:30  bigatti
// - Cleaner code
//
// Revision 1.4  2002/04/15 17:13:31  bigatti
// - Max's new code
//
// Revision 1.3  2002/04/09 14:15:21  bigatti
// - SPoly now takes a GPair as argument (in GPoly.C)
//
// Revision 1.2  2001/12/12 18:21:32  bigatti
// - new structure of reduction
//
// Revision 1.1  2001/12/05 13:01:22  bigatti
// Initial revision
//

