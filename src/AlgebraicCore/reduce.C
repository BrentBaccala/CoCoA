//   Copyright (c)  2005  Anna Bigatti

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


#include "CoCoA/TmpGPoly.H"

#include "CoCoA/IntOperations.H"
#include "CoCoA/ReductionCog.H"
#include "CoCoA/SugarDegree.H"
#include "CoCoA/assert.H"
#include "CoCoA/VectorOperations.H" // for debugging only

// #include <vector>
using std::vector;

namespace CoCoA
{

// Horrible hack, very temporary - profiling only
degree HereForProfilingOnlyWDeg(ConstRefPPMonoidElem cofactor1)
{return wdeg(cofactor1);}



  //-------- headers ----------------------------------------

  // ANNA: where should these go???
  void reduce(ReductionCog& F, const Reductors& v);
  void reduce(ReductionCog& F, const GPoly& g);
  ReductionCog ChooseReductionCogPoly(const GRingInfo& GRI);
  ReductionCog ChooseReductionCogGeobucket(const GRingInfo& GRI);

  //---------------------------------------------------------

  bool ReductorData::myBorelUpdate(ConstRefPPMonoidElem pp, const Reductors& theReductors)
  {
    if (myCount < 100) return false;

    GPoly* f = myGPolyPtr;
    const PPMonoid thePPM(PPM(*f));
    const long BorelIndetIndex = NumIndets(owner(*f))-1;
    PPMonoidElem quot(thePPM);

    thePPM->myDiv(raw(quot), raw(pp), raw(PP(myLPPForDivwMask)));
    if (quot != IndetPower(thePPM, BorelIndetIndex, exponent(quot, BorelIndetIndex)))
      return false;

    f->MultiplyByPP(quot);
    myLPPForDivwMask.myAssign(pp);
    IamBorelUpdated = true;
    f->myReduceTail(theReductors);
    return true;
  }


  inline GPoly* FindReducer(const PPWithMask& pm,
                            const long comp,
                            const Reductors& theReductors)
  {
    for (vector<ReductorData>::iterator it = theReductors.myBorelReductors.begin();
         it != theReductors.myBorelReductors.end();
         ++it)
      if ( comp==it->myComponent && IsDivisibleFast(pm, it->myLPPForDivwMask))
        if (it->IamBorelUpdated || it->myBorelUpdate(PP(pm), theReductors) )
        {
          (it->myCount)++;
          return it->myGPolyPtr;
        }
    for (vector<ReductorData>::const_iterator it = theReductors.myReductors.begin();
         it != theReductors.myReductors.end();
         ++it)
      if (comp==it->myComponent && IsDivisibleFast(pm, it->myLPPForDivwMask)&&!it->IamNotToBeUsed())
      {
        (it->myCount)++;
        return it->myGPolyPtr;
      }
    return 0;
  }


  //---- ReductionCog code --------------------------------------------------

  inline GPoly* FindReducer(const ReductionCog& F, const Reductors& theReductors)
  {
    if ( IsActiveZero(F) ) return 0;
    GRingInfo GRI(theReductors.myGRingInfo());
    //    const PPWithMask pm(ActiveLPP(F), GRI.myDivMaskRule());
    PPWithMask pm(GRI.myPPM(), GRI.myDivMaskRule());
    vector<long> expv;
    exponents(expv, ActiveLPP(F));
    pm = expv;
    return FindReducer(pm, GRI.myComponent(PP(pm)), theReductors);
  }


//   inline int FindReducerIndex(const ReductionCog& F, const vector<RingElem>& v)
//   {
//     if ( IsActiveZero(F) ) return -1;
//     return FindReducerIndex(ActiveLPP(F), v);
//   }


  void ReduceActiveLM(ReductionCog& F, const Reductors& v)
  {
    GPoly* g;
    while ( (g = FindReducer(F, v)) != 0 )
      F->myReduce(poly(*g), len(*g));
  }


  void ReduceActiveLM(ReductionCog& F, SugarDegree& s, const Reductors& v)
  {
    GPoly* g;
    while ( (g = FindReducer(F, v)) != 0 )
    {
      CoCoA_ASSERT( !IsZero(*g));
      s->myUpdate(F, *g);
      F->myReduce(poly(*g), len(*g));
    }//while
  }//ReduceActiveLM


  void reduce(ReductionCog& F, const Reductors& v)
  {
    ReduceActiveLM(F, v);
    while ( !IsActiveZero(F) )
    {
      F->myMoveToNextLM();
      ReduceActiveLM(F, v);
    }
  }


  void reduce(ReductionCog& F, SugarDegree& s, const Reductors& v)
  {
    ReduceActiveLM(F, s, v);
    while ( !IsActiveZero(F) )
    {
      F->myMoveToNextLM();
      ReduceActiveLM(F, s, v);
    }
  }

  void reduce(ReductionCog& F, SugarDegree& s, const GPoly& g)
  {
    if ( IsActiveZero(F) || ActiveLPP(F) < LPPForOrd(g) ) return;
    const GRingInfo& GRI(g.myGRingInfo());
    const PPWithMask& PMg(LPPForDivwMask(g));
    const long Componentg = GRI.myComponent(PP(PMg));
    //    PPWithMask PMF(ActiveLPP(F), GRI.myDivMaskRule());
    PPWithMask PMF(GRI.myPPM(), GRI.myDivMaskRule());
    vector<long> expv;
    exponents(expv, ActiveLPP(F));
    PMF = expv;
    long ComponentF = GRI.myComponent(PP(PMF));
    while ( !IsActiveZero(F) )
    {
      if ( (ComponentF==Componentg && IsDivisibleFast(PMF, PMg) ) )
      {
        s->myUpdate(F, g);
        F->myReduce(poly(g), len(g));// MAX: here the real work is done
      }
      else
        F->myMoveToNextLM();
      if ( IsActiveZero(F) || ActiveLPP(F) < LPPForOrd(g) ) return;
      //      PMF.myAssign(ActiveLPP(F));
      vector<long> expv;
      exponents(expv, ActiveLPP(F));
      PMF = expv;
      ComponentF = GRI.myComponent(PP(PMF));
    }//while
  }//reduce


  void reduce(ReductionCog& F, const GPoly& g)
  {
    // SugarDegree s(sugar(g)); // not used --> use WSugarConst
    SugarDegree s(NewWSugarConst(zero(owner(g))));
    reduce(F, s, g);
  }


  //---- ReductionCog code end ----------------------------------------------


  //-------- ChooseReductionCog... ----------------------------------------

  //ANNA where should these go???
  ReductionCog ChooseReductionCogPoly(const GRingInfo& GRI)
  {
    if ( GRI.myCoeffRingType() == CoeffEncoding::Field )
      return NewRedCogPolyField(GRI.myNewSPR());
    else if ( GRI.myCoeffRingType() == CoeffEncoding::FrFldOfGCDDomain )
      return NewRedCogPolyGCD(GRI.myNewSPR());
    else CoCoA_ERROR("Don't know what to do with these coefficients", "ChooseReductionCog");
    return ReductionCog(0);  // just to keep the compiler quiet
  }


  ReductionCog ChooseReductionCogGeobucket(const GRingInfo& GRI)
  {
    if ( GRI.myCoeffRingType() == CoeffEncoding::Field )
      return NewRedCogGeobucketField(GRI.myNewSPR());
    else if ( GRI.myCoeffRingType() == CoeffEncoding::FrFldOfGCDDomain )
      return NewRedCogGeobucketGCD(GRI.myNewSPR());
    else CoCoA_ERROR("Don't know what to do with these coefficients", "ChooseReductionCog");
    return ReductionCog(0);  // just to keep the compiler quiet
  }


  //-------- GPoly member functions ----------------------------------------

  void GPoly::myPolySetSPoly(const GPoly& f, const GPoly& g)
  {
    myPolyValue = poly(f);
    (owner(f))->myMulByPP(raw(myPolyValue), raw(colon(LPPForOrd(g), LPPForOrd(f))));
    ReductionCog F = ChooseReductionCogPoly(myGRingInfo());
    F->myAssignReset(myPolyValue, f.myLen);
    F->myReduce(poly(g), len(g));
    F->myRelease(myPolyValue);
    // sugar will be set as the sugar of the originating pair
  }


  void GPoly::myReduceTail(const GPoly& g)
  {
    CoCoA_ASSERT( !IsZero(*this) );
    if ( LPPForOrd(*this) < LPPForOrd(g) ) return;
    // geobucket because reduce(F, g) might make many reductions by g
    ReductionCog F = ChooseReductionCogGeobucket(myGRingInfo());
    F->myAssignReset(myPolyValue, myLen);
    F->myMoveToNextLM();
    reduce(F, g);
    F->myRelease(myPolyValue);
    myLen = NumTerms(myPoly()); // LPP,Deg,Comp are the same
  }


  void GPoly::myReduceTail(const Reductors& theReductors)
  {
    CoCoA_ASSERT( !IsZero(*this) );
    ReductionCog F = ChooseReductionCogGeobucket(myGRingInfo());
    F->myAssignReset(myPolyValue, myLen);
    F->myMoveToNextLM();
    reduce(F, mySugar, theReductors);
    F->myRelease(myPolyValue);
    myLen = NumTerms(myPoly()); // LPP,Deg,Comp are the same
  }


  void GPoly::myReduce(const Reductors& theReductors)
  {
    if ( IsZero(*this) ) return;
    ReductionCog F = ChooseReductionCogGeobucket(myGRingInfo());
    F->myAssignReset(myPolyValue, myLen); // myPolyValue gets 0
    reduce(F, mySugar, theReductors); // mySugar updated
    F->myRelease(myPolyValue); // myPolyValue gets the value of F
    myUpdateLenLPPLCDegComp(); // myLen, myWDeg, myComponent are updated

    if ( !IsZero(*this) && !IsOne(myLCValue) ) // makes myPolyValue monic
      if ( myGRingInfo().myCoeffRingType()==CoeffEncoding::Field )
        myGRingInfo().myNewSPR()->myDivByCoeff(raw(myPolyValue), raw(myLCValue));
    // if CoeffEncoding::Field myRelease does NOT make poly monic
    // if CoeffEncoding::FrFldOfGCDDomain myRelease makes poly content free
  }




} // end of namespace CoCoA


// RCS Header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/reduce.C,v 1.23 2014/07/31 14:45:19 abbott Exp $
// $Log: reduce.C,v $
// Revision 1.23  2014/07/31 14:45:19  abbott
// Summary: Merged io.H and UtilsTemplate.H into new header VectorOperations.H
// Author: JAA
//
// Revision 1.22  2014/04/30 16:28:39  abbott
// Summary: Replaced size_t by long
// Author: JAA
//
// Revision 1.21  2012/05/28 09:18:20  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.20  2011/03/11 17:42:21  bigatti
// -- changed  unsigned int --> long
//
// Revision 1.19  2009/11/26 17:24:23  bigatti
// -- just sorted includes
//
// Revision 1.18  2009/11/20 16:06:53  bigatti
// -- just a comment
//
// Revision 1.17  2009/10/27 17:32:40  bigatti
// -- removed temp logging
//
// Revision 1.16  2009/10/27 17:17:51  bigatti
// -- using SugarDegree instead of sugar/wsugar
//
// Revision 1.15  2009/09/25 14:20:01  bigatti
// -- removed temporary code CheckSugar (moved to SugarDegree.H)
//
// Revision 1.14  2009/04/27 13:59:28  bigatti
// -- added SaturatingAlgNoDRL
// -- changed Saturating sugar
//
// Revision 1.13  2009/03/20 17:36:06  bigatti
// -- typo in comment
//
// Revision 1.12  2009/03/20 14:17:24  bigatti
// -- minor fix on (non-tested) SugarSat
//
// Revision 1.11  2009/03/18 17:16:05  bigatti
// -- minor changes, getting ready for SugarDegree checkin
//
// Revision 1.10  2009/02/09 08:17:28  bigatti
// -- just comments and 1 error code
//
// Revision 1.9  2008/09/19 13:46:30  bigatti
// -- added SaturatingAlg (M.Caboara)
//
// Revision 1.8  2008/09/16 15:03:42  bigatti
// -- added LPPForDiv
// -- changed LPP into LPPForOrd
//
// Revision 1.7  2007/12/05 11:06:24  bigatti
// -- changed "size_t StdDeg/myStdDeg(f)" into "long"  (and related functions)
// -- changed "log/myLog(f, i)" into "MaxExponent/myMaxExponent(f, i)"
// -- fixed bug in "IsOne(ideal)" in SparsePolyRing.C
//
// Revision 1.6  2007/12/04 14:27:06  bigatti
// -- changed "log(pp, i)" into "exponent(pp, i)"
//
// Revision 1.5  2007/11/20 08:59:21  bigatti
// -- changed way to compute sugar (uses myDiv instead of operator/)
// -- added comments for naming algorithms
// -- added default in "switch"
//
// Revision 1.4  2007/10/30 17:14:06  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.3  2007/04/27 14:24:15  bigatti
// -- added commented include for <vector>
//
// Revision 1.2  2007/03/09 18:56:56  bigatti
// -- added Tmp prefix to Groebner related files
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.8  2007/03/08 18:22:28  cocoa
// Just whitespace cleaning.
//
// Revision 1.7  2007/03/07 17:04:31  cocoa
// -- several changes by M.Caboara: more operations on ideals,
//    exception cleaner, coding conventions, WSugar, dynamic
//
// Revision 1.6  2007/02/10 18:44:03  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.5  2006/12/04 13:55:54  cocoa
// -- added: sugar for GradingDim > 0  (called wsugar)
//
// Revision 1.4  2006/10/11 13:33:06  cocoa
// -- rearranged code for sugar in reduce.C
// -- activated sugar in GPair.C
// -- removed operator<< for GPairList (template in io)
//
// Revision 1.3  2006/10/06 16:32:06  cocoa
// -- changed: GPoly::SPoly --> GPoly::myAssignSPoly
// -- changed: Len(const GPoly&) --> len(const GPoly&)
// -- added:   poly(const GPoly&)
// -- added:   GPoly::myUpdateLenLPPDegComp()
// -- in reduce.C added functions for computing sugar during reduction
//
// Revision 1.2  2006/10/06 14:04:14  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
