//   Copyright (c)  2005,2009  John Abbott

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


// Source code for abstract class SparsePolyRing and friends

#include "CoCoA/TmpMonomialIdeal.H"
#include "CoCoA/IntOperations.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/PPWithMask.H"  // for monomial ideals
#include "CoCoA/QuotientRing.H" // for IsQuotientRing
#include "CoCoA/RingDistrMPolyClean.H" // for NewPolyRing_DMP
#include "CoCoA/RingFp.H" // for IsRingFp
#include "CoCoA/RingTwinFloat.H" // for IsRingTwinFloat
#include "CoCoA/TmpGOperations.H"  // for myGcd
#include "CoCoA/TmpPPVector.H"  // for interreduce(PPs)
#include "CoCoA/assert.H"
#include "CoCoA/convert.H"
#include "CoCoA/degree.H"
#include "CoCoA/error.H"
#include "CoCoA/ideal.H"     // for myGcd
#include "CoCoA/symbol.H"


#include <algorithm>
using std::swap;
#include <functional>
using std::not1;    // for IsMonomial
using std::ptr_fun; // for IsMonomial
#include <iostream>
// using std::ostream in SparsePolyRingBase::myOutput
#include <list>
// using std::list;
#include <string>
// using std::string;
//#include <vector>
using std::vector;

namespace CoCoA
{

  const std::vector<RingElem>& SparsePolyRingBase::IdealImpl::myGBasisMonId() const
  //  void MonomialInterreduce(std::vector<RingElem>& res, const std::vector<RingElem>& l)
  {
    CoCoA_ASSERT(IhaveMonomialGens());
    if (myGBasisIsValid) return myGBasisValue;
    CoCoA_ASSERT(myGBasisValue.empty());
    if (IamZero()) return myGBasisValue;
    // convert input into PPVector, operate, convert back
    const SparsePolyRing P = myRing();
    PPVector g(PPM(P), NewDivMaskEvenPowers());
//     const std::vector<RingElem>& l = myGens();
//     for (vector<RingElem>::const_iterator it=l.begin(); it!=l.end() ; ++it)
//       PushBack(g, LPP(*it));
    convert(g, myGens());
    interreduce(g);
    convert(myGBasisValue, P, g);
    myGBasisIsValid = true;
    // myMinGensValue = myGBasisValue; // not necessary
    return myGBasisValue;
  }


  void SparsePolyRingBase::IdealImpl::myIntersectMonId(const ideal& J)
//  void MonomialIntersection(std::vector<RingElem>& res, const std::vector<RingElem>& l1,  const std::vector<RingElem>& l2)
  {
    CoCoA_ASSERT(IhaveMonomialGens());
    CoCoA_ASSERT(AreGensMonomial(J));
    if (IamZero()) return;
    std::vector<RingElem> res;
    if (!gens(J).empty())
    {
      // convert input into PPVector, operate, convert back
      const SparsePolyRing P = myRing();
      PPVector g1(PPM(P), NewDivMaskEvenPowers());
      const std::vector<RingElem>& l1 = myGens();
      for (vector<RingElem>::const_iterator it=l1.begin(); it!=l1.end() ; ++it)
        PushBack(g1, LPP(*it));
      PPVector g2(PPM(P), DMR(g1));
      const std::vector<RingElem>& l2 = gens(J);
      for (vector<RingElem>::const_iterator it=l2.begin(); it!=l2.end() ; ++it)
        PushBack(g2, LPP(*it));
      PPVector g(PPM(P), DMR(g1));
      g.myLcms(g1, g2);
      interreduce(g);
      convert(res, P, g);
    }
    // assign into this ideal
    // IhaveMonomialGens3Value = true; // stays constant
    if ((!IsTrue3(IhaveSquareFreeMonomialGens3Value))
        || (!IsTrue3(ourGetPtr(J)->IhaveSquareFreeMonomialGens3Value)))
      IhaveSquareFreeMonomialGens3Value = uncertain3; // can do better than this..
    std::swap(myGensValue, res); // assignment
    myGBasisValue = myGensValue;
    myGBasisIsValid = true;
    // myMinGensValue = myGBasisValue; // not necessary
  }


void SparsePolyRingBase::IdealImpl::myColonMonId(const ideal& J)
{
    CoCoA_ASSERT(IhaveMonomialGens());
    CoCoA_ASSERT(AreGensMonomial(J));
    if (IamZero()) return;
    std::vector<RingElem> res;
    if (!gens(J).empty())
    {
      // convert input into PPVector, operate, convert back
      const SparsePolyRing P = myRing();
      PPVector g1(PPM(P), NewDivMaskEvenPowers());
      const std::vector<RingElem>& l1 = myGens();
      for (vector<RingElem>::const_iterator it=l1.begin(); it!=l1.end() ; ++it)
        PushBack(g1, LPP(*it));
      PPVector g2(PPM(P), DMR(g1));
      const std::vector<RingElem>& l2 = gens(J);
      for (vector<RingElem>::const_iterator it=l2.begin(); it!=l2.end() ; ++it)
        PushBack(g2, LPP(*it));
      PPVector g(PPM(P), DMR(g1));
      PPVector tmp(PPM(g), DMR(g));
      const long len1 = len(g1);
      const long len2 = len(g2);
      for(long i=0; i<len2; ++i)
      {
        tmp.myClear();
        for (long j=0; j<len1; ++j)
          PushBack(tmp, colon(PP(g1[j]),PP(g2[i])));
        interreduce(tmp);
        if (i==0) swap(g,tmp);
        else lcms(g, g, tmp);
        interreduce(g);
      }
      convert(res, P, g);
    }
    // assign into this ideal
    // IhaveMonomialGens3Value = true; // stays constant
    std::swap(myGensValue, res); // assignment
    if (!IsTrue3(IhaveSquareFreeMonomialGens3Value))
      IhaveSquareFreeMonomialGens3Value = uncertain3; // can do better than this..
    myGBasisValue = myGensValue;
    myGBasisIsValid = true;
    // myMinGensValue = myGBasisValue; // not necessary
}
  

  void SparsePolyRingBase::IdealImpl::myMulMonId(const ideal& J)
  {
    CoCoA_ASSERT(IhaveMonomialGens());
    CoCoA_ASSERT(AreGensMonomial(J));
    if (IamZero()) return;
    std::vector<RingElem> res;
    if (!gens(J).empty())
    {
      // convert input into PPVector, operate, convert back
      const SparsePolyRing P = myRing();
      PPVector g1(PPM(P), NewDivMaskEvenPowers());
      const std::vector<RingElem>& l1 = myGens();
      for (vector<RingElem>::const_iterator it=l1.begin(); it!=l1.end() ; ++it)
        PushBack(g1, LPP(*it));
      PPVector g2(PPM(P), DMR(g1));
      const std::vector<RingElem>& l2 = gens(J);
      for (vector<RingElem>::const_iterator it=l2.begin(); it!=l2.end() ; ++it)
        PushBack(g2, LPP(*it));
      PPVector g(PPM(P), DMR(g1));
      for(long i=len(g2)-1; i>=0; --i)
        for (long j=len(g1)-1; j>=0; --j)
          PushBack(g, PP(g1[j])*PP(g2[i]));
      interreduce(g);
      convert(res, P, g);
    }
    // assign into this ideal
    // IhaveMonomialGens3Value = true; // stays constant
    IhaveSquareFreeMonomialGens3Value = uncertain3; // can do better than this..
    std::swap(myGensValue, res); // assignment
    myGBasisValue = myGensValue;
    myGBasisIsValid = true;
    // myMinGensValue = myGBasisValue; // not necessary
  }

  
  void SparsePolyRingBase::IdealImpl::myElimMonId(const std::vector<RingElem>& ELimIndets)
  {
    CoCoA_ASSERT(IhaveMonomialGens());
    if (IamZero()) return;
    std::vector<RingElem> res;
    {
      // convert input into PPVector, operate, convert back
      const SparsePolyRing P = myRing();
      PPVector g1(PPM(P), NewDivMaskEvenPowers());
      const std::vector<RingElem>& l1 = myGens();
      for (vector<RingElem>::const_iterator it=l1.begin(); it!=l1.end() ; ++it)
        PushBack(g1, LPP(*it));
      PPVector g2(PPM(P), DMR(g1));
      const std::vector<RingElem>& l2 = ELimIndets;
      for (vector<RingElem>::const_iterator it=l2.begin(); it!=l2.end() ; ++it)
        PushBack(g2, LPP(*it));
      PPVector g(PPM(P), DMR(g1));
      const long len1 = len(g1);
      for(long i=0; i<len1; ++i)
        if (!IsDivisible(g1[i], g2)) PushBack(g, g1[i]);
      interreduce(g);
      convert(res, P, g);
    }
    // assign into this ideal
    // IhaveMonomialGens3Value = true; // stays constant
    if (!IsTrue3(IhaveSquareFreeMonomialGens3Value))
      IhaveSquareFreeMonomialGens3Value = uncertain3;
    std::swap(myGensValue, res); // assignment
    myGBasisValue = myGensValue;
    myGBasisIsValid = true;
  }  


  ideal IndetsIdeal(const PolyRing& P, ConstRefPPMonoidElem pp)
  {
    vector<RingElem> g;
    for (long i=0 ; i < NumIndets(owner(pp)) ; ++i )
      if ( exponent(pp, i) != 0 )
      {
        if ( exponent(pp, i) != 1 )
          CoCoA_ERROR("input must be square free", "IndetsIdeal");
        g.push_back(indet(P, i));
      }
    return ideal(P, g);
    // IhaveMonomialGens3Value = true;
    // IhaveSquareFreeMonomialGens = true;
  }


  ideal AlexanderDual(const ideal& I)
  {
    const SparsePolyRing P = RingOf(I);
    if (!AreGensMonomial(I))
      CoCoA_ERROR("not monomial", "AlexanderDual");
    if (!AreGensSquareFreeMonomial(I))
      CoCoA_ERROR(ERR::NYI, "AlexanderDual");
    DivMaskRule DMR = NewDivMaskEvenPowers();
    PPVector g(PPM(P), DMR);
    PPVector AD(PPM(P), DMR);
    const std::vector<RingElem>& l1 = gens(I);
    for (vector<RingElem>::const_iterator it=l1.begin(); it!=l1.end() ; ++it)
      PushBack(g, LPP(*it));
    AD.myAlexanderDual(g);
    vector<RingElem> res;
    convert(res, P, AD);
    return ideal(P, res);
  }


  vector<ideal> PrimaryDecompositionMonId(const ideal& I)
  {
    const SparsePolyRing P = RingOf(I);
    if (!AreGensSquareFreeMonomial(I))
      CoCoA_ERROR(ERR::NYI, "PrimaryDecompositionMonId");
    vector<RingElem> g = gens(AlexanderDual(I));
    vector<ideal> PD;
    for (vector<RingElem>::const_iterator it=g.begin(); it!=g.end() ; ++it)
      PD.push_back(IndetsIdeal(P, LPP(*it)));
    return PD;
  }

} // end of namespace CoCoA

// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/TmpMonomialIdeal.C,v 1.13 2014/07/30 14:11:41 abbott Exp $
// $Log: TmpMonomialIdeal.C,v $
// Revision 1.13  2014/07/30 14:11:41  abbott
// Summary: Changed name AmbientRing --> RingOf; myAmbientRing --> myRing
// Author: JAA
//
// Revision 1.12  2014/07/14 15:08:59  abbott
// Summary: Removed include of tmp.H (no longer needed?)
// Author: JAA
//
// Revision 1.11  2014/07/08 15:23:53  abbott
// Summary: Updated comment
// Author: JAA
//
// Revision 1.10  2014/07/07 13:14:25  abbott
// Summary: Removed AsSparsePolyRing
// Author: JAA
//
// Revision 1.9  2014/03/27 14:58:05  bigatti
// -- just comments to remind that GBasis is minimal generators
//
// Revision 1.8  2012/05/30 13:44:45  bigatti
// -- renamed IhaveMonomialGensB3Value --> IhaveMonomialGens3Value
//
// Revision 1.7  2012/05/29 07:45:23  abbott
// Implemented simplification change to bool3:
//  changed names of the constants,
//  changed names of the testing fns.
//
// Revision 1.6  2012/05/28 09:18:20  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.5  2012/02/10 10:29:07  bigatti
// -- changed RingZ.H, RingQ.H --> RingZZ.H, RingQQ.H
//
// Revision 1.4  2011/11/07 11:06:14  bigatti
// -- myGBasisMonId: using now convert(g, myGens()); instead of extended code
//
// Revision 1.3  2011/07/27 15:50:56  bigatti
// -- improved use of len in for loops
//
// Revision 1.2  2011/07/05 15:02:17  bigatti
// -- added AlexanderDual
// -- added ad-hoc functions for colon, elim on monomial ideals
//
// Revision 1.1  2011/06/27 13:26:51  bigatti
// -- first import (soem functions were in SparsePolyRing)
//
