//   Copyright (c)  2005-2007  Massimo Caboara

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


#include "CoCoA/TmpGReductor.H"
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/FractionField.H"
#include "CoCoA/FreeModule.H"
#include "CoCoA/IntOperations.H"
#include "CoCoA/MatrixOperations.H"
#include "CoCoA/MatrixForOrdering.H"
#include "CoCoA/MatrixView.H"
#include "CoCoA/ModuleOrdering.H"
#include "CoCoA/RingDistrMPolyInlFpPP.H"
#include "CoCoA/RingDistrMPolyInlPP.H"
#include "CoCoA/RingFp.H" // for dynamic_cast<RingFpImpl*>(CoeffRing.myRingPtr())
#include "CoCoA/RingQQ.H"
#include "CoCoA/apply.H"
#include "CoCoA/assert.H"
#include "CoCoA/convert.H"
#include "CoCoA/VectorOperations.H"  // for templates which print lists/vectors (in myStampaPairs)
#include "CoCoA/matrix.H"
#include "CoCoA/symbol.H"
#include "CoCoA/time.H"
#include "CoCoA/VectorOperations.H" // for HasUniqueOwner
//#include "CoCoA/utils.H" // for LongRange

#include <algorithm>
using std::for_each;
using std::find_if;  // for IsHomog
#include <list>
using std::list;
#include <functional>
using std::binary_function;
using std::less;
using std::mem_fun_ref; // for calling GPair::complete on GPairList
using std::not1;    // for IsHomog
using std::ptr_fun; // for IsHomog
#include <iostream>
using std::ostream;
using std::clog;
using std::endl;
using std::flush;
#include <utility>
using std::make_pair;
#include <vector>
using std::vector;

//static bool ANNA_DEBUG = true;
static const bool MAX_DEBUG = false;


namespace CoCoA
{
  int GReductor::ourDefaultStatLevel = -1;

  // This function produces a new PolyRing using
  // NewPolyRing_DMPII or NewPolyRing_DMPI
  // according to the char.
  // Used internally, indets names are always x[]
  SparsePolyRing NewSparsePolyRing(const ring& CoeffRing,
                                   const vector<symbol>& IndetNames,
                                   const PPOrdering& ord)
  {
    if (IsRingFp(CoeffRing)) // special case for SmallFpImpl
      return NewPolyRing_DMPII(CoeffRing, IndetNames, ord);
    else
      return NewPolyRing_DMPI(CoeffRing, IndetNames, ord);
  }//NewSparsePolyRing


  const GRingInfo& GetGRI(const GPolyList& theGPL)
  {
    CoCoA_ASSERT(!theGPL.empty());
    return theGPL.begin()->myGRingInfo();
  }//GRI


/*  perhaps useless
bool BoolCmpLPPPoly(ConstRefRingElem f, ConstRefRingElem g)
  {
    CoCoA_ASSERT(!IsZero(f));//BUG HUNTING  ???
    CoCoA_ASSERT(!IsZero(g));//BUG HUNTING  ???
    const SparsePolyRing P = owner(f);
    return P->myCmpLPP(raw(f), raw(g))>0;
  }
  */

  bool BoolCmpLPPGPoly(const GPoly& f, const GPoly& g)
  {
    CoCoA_ASSERT(!IsZero(f));//BUG HUNTING  ???
    CoCoA_ASSERT(!IsZero(g));//BUG HUNTING  ???
    //const PPMonoid& PPM1=PPM(owner(f));
    return PPM(owner(f))->myCmp(raw(LPPForOrd(f)),raw(LPPForOrd(g)))>0;
  }

  void GReductor::myCtorAux(const BuchbergerOpTypeFlag theBuchbergerOpType,
                            const UseDynamicAlgFlag IsDynamic)
  {
    myPrepared=false;
    myAgeValue = 0;
    myWrongLPPFoundValue=false;

    myBuchbergerOpType=theBuchbergerOpType;
    IsDynamicAlgorithm=IsDynamic;

    if (!myCriteria.myCoprime) myStat.myPrintCop=false;
    if (!myCriteria.myGM) myStat.myPrintGM=false;
    if (!myCriteria.myBack) myStat.myPrintBC=false;
  }

  GReductor::GReductor(const GRingInfo& theGRI,
                       const PolyList& TheInputPolys,
                       const int StatLevel,
                       const BuchbergerOpTypeFlag theBuchbergerOpType,
                       const Reductors::UseBorelFlag UBR,
                       const UseDynamicAlgFlag IsDynamic,
                       const GBCriteria criteria):
      myGRingInfoValue(theGRI),
      myTrueReductors(theGRI, UBR),
      mySPoly(theGRI),
      myOldDeg(GradingDim(theGRI.myNewSPR())),
      myCurrentPairDeg(GradingDim(theGRI.myNewSPR())),
      myStat(len(TheInputPolys), StatLevel),
      myCriteria(criteria)
  {
    //myNReductions=0;
    myCtorAux(theBuchbergerOpType, IsDynamic);

    for (PolyList::const_iterator it=TheInputPolys.begin();it!=TheInputPolys.end();++it)
    {
      if (!IsZero(*it))
        myPolys.push_back(GPoly(*it,myGRingInfoValue));
    }
    myPolys.sort(BoolCmpLPPGPoly);

  // myNewReductors is void at start;
    // These are for fine tuning of the GReductor statistics
    //myStat.myPrintReduction=true;
    //myStat.myPrintDeg=true;
    //myStat.myPrintCop=true;
    //myStat.myPrintGM=true;
    //myStat.myPrintBC=true;
    //myStat.myPrintNumPair=true;
    //myStat.myPrintFinal=true;
    //myStat.myPrintFinalFull=true;
    //myStat.myPrintNewPairs=true;
    //myStat.myPrintPolyDeleted=true;
    //myStat.myPrintPolyDH=true;
    //myStat.myPrintPolyLen=true;
  }// GReductor ctor


  // This ctor allows to avoid transforming gpolys to polys and back if
  // you want to evaluate complex expressions. The input GPolys may contain
  // zeros and/or not be sorted. Take care.
  GReductor::GReductor(const GRingInfo& theGRI,
                       const GPolyList& TheInputGPolys,
                       const int StatLevel,
                       const BuchbergerOpTypeFlag theBuchbergerOpType,
                       const Reductors::UseBorelFlag UBR,
                       const UseDynamicAlgFlag IsDynamic,
                       const GBCriteria criteria):
      myGRingInfoValue(theGRI),
      myTrueReductors(theGRI, UBR),
      mySPoly(theGRI),
      myOldDeg(GradingDim(theGRI.myNewSPR())),
      myCurrentPairDeg(GradingDim(theGRI.myNewSPR())),
      myStat(len(TheInputGPolys), StatLevel),
      myCriteria(criteria)
  {
    myPolys=TheInputGPolys;
    //myNReductions=0;
    myCtorAux(theBuchbergerOpType, IsDynamic);

    // If you want to remove zeros and/or sort your input GPolys.
    // GPoly Zero(zero(myPolyRing),myGRingInfoValue,Component(LPP(zero(myPolyRing),myGRingInfoValue)),0);
    // myPolys.remove(Zero);
    // myPolys.sort(BoolCmpLPPGPoly);
  } // GReductor ctor


  // This ctor allows to avoid transforming gpolys to polys and back if
  // you want to evaluate complex expressions. The input GPolys may contain
  // zeros and/or not be sorted. Take care.
  GReductor::GReductor(const GRingInfo& theGRI,
                       GPolyList& TheInputGPolys,
                       const ClearMarker,// just for ctor ariety
                       const int StatLevel,
                       const BuchbergerOpTypeFlag theBuchbergerOpType,
                       const Reductors::UseBorelFlag UBR,
                       const UseDynamicAlgFlag IsDynamic,
                       const GBCriteria criteria):
      myGRingInfoValue(theGRI),
      myTrueReductors(theGRI, UBR),
      mySPoly(theGRI),
      myOldDeg(GradingDim(theGRI.myNewSPR())),
      myCurrentPairDeg(GradingDim(theGRI.myNewSPR())),
      myStat(len(TheInputGPolys), StatLevel),
      myCriteria(criteria)
  {
    for (GPolyList::iterator it=TheInputGPolys.begin();it!=TheInputGPolys.end();++it)
    {
      myPolys.push_back(GPoly(theGRI));
      myPolys.back().AssignClear(*it);
    }
    TheInputGPolys.clear();
    //myNReductions=0;
    myCtorAux(theBuchbergerOpType, IsDynamic);
    // If you want to remove zeros and/or sort your input GPolys.
    // GPoly Zero(zero(myPolyRing),myGRingInfoValue,Component(LPP(zero(myPolyRing),myGRingInfoValue)),0);
    // myPolys.remove(Zero);
    // myPolys.sort(BoolCmpLPPGPoly);
  } // GReductor ctor






//   // copy ctor - NOT WORKING, see comment below
//   GReductor::GReductor(const GReductor& theGR):
//       myGRingInfoValue(theGR.myGRingInfoValue),
//       myTrueReductors(theGR.myGRingInfoValue,theGR.myTrueReductors.IhaveBorelReductors()),
//       mySPoly(theGR.myGRingInfoValue),
//       myOldDeg(GradingDim(theGR.myGRingInfoValue.myNewSPR())),
//       myCurrentPairDeg(GradingDim(theGR.myGRingInfoValue.myNewSPR())),
//       myStat(len(theGR.myPolys),theGR.myStat.myGetLevel()),
//       myCriteria(theGR.myCriteria)
//   {
//     myPolys=theGR.myPolys;
//     // All the things that points to some element of myPolys
//     // has to be redo.
//     // myPairs,myTrueReductors,myGB
//     myPrepared=theGR.myPrepared;
//     myAgeValue = theGR.myAgeValue;
//     //myNReductions=theGR.myNReductions;
//     IsDynamicAlgorithm=theGR.IsDynamicAlgorithm;
//     myWrongLPPFoundValue=theGR.myWrongLPPFoundValue;
//     myBuchbergerOpType=theGR.myBuchbergerOpType;
//     if (!myCriteria.myCoprime) myStat.myPrintCop=false;
//     if (!myCriteria.myGM) myStat.myPrintGM=false;
//     if (!myCriteria.myBack) myStat.myPrintBC=false;
//     myStat=theGR.myStat;
//   }//GReductor copy ctor


  // To be filled out
  ostream& operator<<(ostream& out, const GReductor& GR)
  {
    out<<"The GROBNER REDUCTOR\n";
    out<<"  Reductors Len="<<GR.myReductorsLen()
        <<"  GB Len="<<GR.myGBasisLen()
        <<"  Pairs Len="<<GR.myPairsLen()
        <<"  Byte Size="<<sizeof(GR)<<"\n";
    GR.myStampaReductors(out);
    GR.myStampaGB(out);
    GR.myStampaPairs(out);
    out<<"\n";
    out<<"The Ring"<<GR.myGRingInfoValue<<"\n";
    //matrix OrdM(NewDenseMat(Z,
    //                             NumIndets(GR.myGRingInfoValue.myNewSPR())),
    //                             NumIndets(GR.myGRingInfoValue.myNewSPR())));
    //(PPM(GR.myGRingInfoValue.myNewSPR())->myOrdering())->myOrdMatCopy(OrdM);
    //out<<OrdM;
    out<<"GradingDim  is "<<GradingDim(GR.myGRingInfoValue.myNewSPR())<<endl;
    out<<"The Ring Special Index "<<ModuleVarIndex(GR.myGRingInfoValue.myNewSPR())<<"\n";
    //out<<"Reductions performed "<<GR.myNReductions<<"\n";
    out<<"Age "<< GR.myAgeValue <<"\n";
    out<<"Preparation done? "<<GR.myPrepared<<"\n";
    out<<"myOldDeg "<<GR.myOldDeg<<"\n";
    out<<"myCurrentPairDeg "<<GR.myCurrentPairDeg<<"\n";
    out<<"Is Dynamic Algorithm? "<<GR.IsDynamicAlgorithm<<"\n";
    out<<"Is Wrong LPP been Found? "<<GR.myWrongLPPFoundValue<<"\n";
    out<<"Cop Criterion " <<GR.myCriteria.myCoprime<<"\n";
    out<<"GM Criteria "   <<GR.myCriteria.myGM<<"\n";
    out<<"Back Criterion "<<GR.myCriteria.myBack<<"\n";
    out<<"Div Criterion " <<GR.myCriteria.myDiv<<"\n";
    out<<"Algorithm "<<GR.myBuchbergerOpType<<"\n";
    out<<"\n";
    return out;
  }

  // Can be used if the Greductor has been updated. If not, you are missing the
  // Spoly
  void GReductor::GetCandidateGBasis(PolyList& thePL)const
  {
    PolyList PL;
    for (GPolyPtrList::const_iterator it=myGB.begin();it!=myGB.end();++it)
      if ((**it).IsActive())
        PL.push_back((**it).myPoly());
    //if (!IsZero(mySPoly))
    //  PL.push_back(GetSPoly().myPoly());
    for (GPairList::const_iterator it=myPairs.begin();it!=myPairs.end();++it)
      if (it->IsInputPoly())
        PL.push_back(it->myFirstGPoly().myPoly());
    swap(PL,thePL);
  }//GetCandidateGBasis

  void GReductor::myStampaGB(ostream& out)const
  {
    out<<"The GBASIS\n";
    list<GPoly*>::const_iterator it=myGB.begin();
    for (;it!=myGB.end();++it) out<<**it<<endl;
    out<<endl;
  }//myStampaGB



  void GReductor::myStampaPPGB(ostream& out)const
  {
    out<<"The GBASIS\n";
    list<GPoly*>::const_iterator it=myGB.begin();
    for (;it!=myGB.end();++it) out<<LPPForDiv(**it)<<",\n";
    out<<"\n";
  }//myStampaPPGB


  void GReductor::myStampaPairs(ostream& out)const
  {
    out<<"The PAIRS\n"<<myPairs<<"\n";
  }//myStampaPairs


  void GReductor::myStampaReductors(ostream& out)const
  {
    myTrueReductors.myStampaReductors(out);
    out<<endl;
  }//myStampaReductors

// This procedure may be substituted by a transform_if
  void GReductor::myGBasis(PolyList& GBasis)
  {
    GBasis.clear();
    GPolyPtrList::const_iterator it;
    for (it=myGB.begin();it!=myGB.end();++it)
      if ((**it).IsActive()) GBasis.push_back((**it).myPoly());
  }//myGBasis


  void GReductor::myMinGens(PolyList& MinGens)
  {
    MinGens.clear();
    GPolyPtrList::const_iterator it;
    for (it=myGB.begin();it!=myGB.end();++it)
      if (IsMinimalGen(**it)) MinGens.push_back((**it).myPoly());
  }//myGBasis


 // This procedure may be substituted by a transform_if
  void GReductor::myGBasis(GPolyList& GBasis)
  {
    GBasis.clear();
    GPolyPtrList::const_iterator it;
    for (it=myGB.begin();it!=myGB.end();++it)
      if ((**it).IsActive()) GBasis.push_back((**it));
  }//myGBasis

// This procedure may be substituted by a transform_if
  void GReductor::myGBasisClear(GPolyList& theGBasis)
  {
    theGBasis.clear();
    GPolyPtrList::const_iterator it;
    for (it=myGB.begin();it!=myGB.end();++it)
      if ((**it).IsActive())
      {
        GPoly Zero(myGRingInfoValue);
        theGBasis.push_back(Zero);
        theGBasis.back().AssignClear(**it);
      }
    myGB.clear();
  }//myGBasis


  void GReductor::myGBasis(VectorList& GBasis)
  {
    GBasis.clear();
    if (myGB.empty())
      return;
    GPolyPtrList::const_iterator it;
    for (it=myGB.begin();it!=myGB.end();++it)
      if ((**it).IsActive())
        GBasis.push_back(DeEmbedPoly((**it).myPoly(),myGRingInfoValue));
  }


  void GReductor::myMinGens(VectorList& MinGens)
  {
    MinGens.clear();
    GPolyPtrList::const_iterator it;
    for (it=myGB.begin();it!=myGB.end();++it)
      if (IsMinimalGen(**it))
        MinGens.push_back(DeEmbedPoly((**it).myPoly(),myGRingInfoValue));
  }//myGBasis



//esame piu' approfondito - sia correttezza sia efficienza
//warning: il CopCriterion is false, still the coprimality is used to decide
//         which pair kill
  long GReductor::myGMInsert(GPairList& L,GPair P)
  {
    long NumPairsConsidered=0;
    bool ToBeInserted=true;
    bool erased=false;
    long P_Component=GPairComponent(P);

    GPairList::iterator it=L.begin();
    while (it!=L.end() && ToBeInserted)
    {
      ++NumPairsConsidered;
      if (P_Component==GPairComponent(*it)
          && IsDivisibleFast(LCMwMask(P), LCMwMask(*it)))
      {
        ToBeInserted=false;
        if (LCMwMask(P)==LCMwMask(*it) && (!it->IamCoprime()) && (P.IamCoprime()))
        {
          //          P.myComplete();          it->myComplete();
          *it=P;
        }
      }
      else
        if (P_Component==GPairComponent(*it)
            && IsDivisibleFast(LCMwMask(*it), LCMwMask(P)) )
        {
          it=L.erase(it);
          erased=true;
        }
      if (!erased)  ++it;
      erased=false;
    }//while
    if (ToBeInserted)
      L.push_back(P);// new pairs are sorted later
    return NumPairsConsidered;
  }


  void GReductor::myBuildNewPairs(GPairList& new_pairs)
  {
    long standing_index=len(myGB)-1;
    long walking_index=0;
    long inserted_pairs=0;//STAT

    GPolyPtrList::const_iterator last=myGB.end(); --last;
    long last_component=Component(**last);
    GPolyPtrList::const_iterator it;
    for (it=myGB.begin();it!=last;++it,++walking_index)
    {
      if ((**it).IsActive()&&last_component==Component(**it))
      {
        //std::clog<<"walking_component "<<Component(**it)<<endl;
        if (myCriteria.myDiv && IsDivisibleFast(LPPForDivwMask(**it), LPPForDivwMask(**last)) ) // speed is not necessary
        {
          if (myStat.myPrintPolyDeleted)
            clog<<"<"<<walking_index<<"> "<<LPPForDiv(**it)<<" DELETED BY NEW "
                <<"<"<<standing_index<<"> "<<LPPForDiv(**last)<<"\n";
          (*it)->Deactivate();
          ++myStat.myPolyDeleted;
        } //second if
        //else
        ++inserted_pairs;
        ++myStat.myPInserted;
        if (myCriteria.myGM)
          myStat.myGMTouched+=myGMInsert(new_pairs,
                                         GPair(**it,**last,walking_index,standing_index));
        else
        {
          myStat.myGMTouched=0;
          Ordered_Insert(new_pairs,GPair(**it,**last,walking_index,standing_index));
        }
      };//Active if
    };//for

    myStat.myGMKilled+=inserted_pairs-len(new_pairs);
    if (myStat.myPrintGM && (inserted_pairs!=len(new_pairs)))
      clog<<"[GM KILLED "<<inserted_pairs-len(new_pairs)
          <<" OUT OF "<<inserted_pairs<<"]\n";

    long pre_cop_test=len(new_pairs);
    if (myCriteria.myCoprime)
    {
      GPairList::iterator it1=new_pairs.begin();
      while (it1!=new_pairs.end())
        if (it1->IamCoprime())
          it1=new_pairs.erase(it1);
        else
          ++it1;
    }
    if (pre_cop_test!=len(new_pairs))
      if (myStat.myPrintCop)
        clog<<"[COP KILLED "<<pre_cop_test-len(new_pairs)
            <<" OUT OF "<<inserted_pairs<<"]\n";
    myStat.myCopKilled+=pre_cop_test-len(new_pairs);
    if ( myStat.myPrintNewPairs)
    {
      clog<<"NEW PAIRS "<<len(new_pairs)<<"  ";
      for (GPairList::const_iterator it=new_pairs.begin();it!=new_pairs.end();++it)
        clog<<"<"<<it->myFirstIndex()<<","<<it->mySecondIndex()<<">, ";
      clog<<"\n";
    };
  }//myBuildNewPairs

   void GReductor::myBuildNewPairsAll(GPairList& new_pairs)
  {
    new_pairs.clear();
    long standing_index=1;
    long walking_index=0;
    long inserted_pairs=0;//STAT

    GPolyPtrList::const_iterator BeginPlus=myGB.begin(); ++BeginPlus;
    for (GPolyPtrList::const_iterator last=BeginPlus;last!=myGB.end();++last,++standing_index)
    {
    long last_component=Component(**last);
    GPolyPtrList::const_iterator it;
    walking_index=0;
    for (it=myGB.begin(); it!=last; it++,walking_index++)
    {
      if ((**it).IsActive()&&last_component==Component(**it))
      {
        //clog<<"walking_component "<<Component(**it)<<endl;
        if (myCriteria.myDiv && IsDivisibleFast(LPPForDivwMask(**it), LPPForDivwMask(**last)) )
        {
          if (myStat.myPrintPolyDeleted)
            clog<<"<"<<walking_index<<"> "<<LPPForDiv(**it)<<" DELETED BY NEW "
                <<"<"<<standing_index<<"> "<<LPPForDiv(**last)<<"\n";
          (*it)->Deactivate();
          myStat.myPolyDeleted++;
        } //second if
        //else
        inserted_pairs++;
        myStat.myPInserted++;
        if (myCriteria.myGM)
          myStat.myGMTouched+=myGMInsert(new_pairs,
                                         GPair(**it,**last,walking_index,standing_index));
        else
        {
          myStat.myGMTouched=0;
          Ordered_Insert(new_pairs,GPair(**it,**last,walking_index,standing_index));
        }
      };//Active if
    };//walking for
    }// last for
    myStat.myGMKilled+=inserted_pairs-len(new_pairs);
    if (myStat.myPrintGM && (inserted_pairs!=len(new_pairs)))
      clog<<"[GM KILLED "<<inserted_pairs-len(new_pairs)
          <<" OUT OF "<<inserted_pairs<<"]\n";

    long pre_cop_test=len(new_pairs);
    if (myCriteria.myCoprime)
    {
      GPairList::iterator it1=new_pairs.begin();
      while (it1!=new_pairs.end())
        if (it1->IamCoprime())
          it1=new_pairs.erase(it1);
        else
          ++it1;
    }
    if (pre_cop_test!=len(new_pairs))
      if (myStat.myPrintCop)
        clog<<"[COP KILLED "<<pre_cop_test-len(new_pairs)
            <<" OUT OF "<<inserted_pairs<<"]\n";
    myStat.myCopKilled+=pre_cop_test-len(new_pairs);
    if ( myStat.myPrintNewPairs)
    {
      clog<<"NEW PAIRS "<<len(new_pairs)<<"  ";
      for (GPairList::const_iterator it=new_pairs.begin();it!=new_pairs.end();it++)
        clog<<"<"<<it->myFirstIndex()<<","<<it->mySecondIndex()<<">, ";
      clog<<"\n";
    };
  }//myBuildNewPairsAll

  void GReductor::myApplyBCriterion()
  {
    const PPWithMask& newPP(LPPForDivwMask(myPolys.back()));
    const long newPP_component=Component(myPolys.back());
    GPairList::iterator it=myPairs.begin();
    while (it!=myPairs.end())
    {
      ++myStat.myBTouched;
      if (GPairComponent(*it)!=newPP_component
          ||
          it->BCriterion_OK(newPP))
      {++it;}
      else
      {
        if (myStat.myPrintPolyDeleted)
          clog<<*it<<" KILLED BY NEW POLY "<<newPP<<" BCRIT\n";
        ++myStat.myBKilled;
        it=myPairs.erase(it);
      };// else
    };//While
  }

  void GReductor::myUpdateBasisAndPairs()
  {
    // HERE NO MORE - in DOGBASIS or DOAHGBASIS If non homog, probably no interred
    //if (!IsTrueGCDDomain(CoeffRing(mySPoly)))
    //myTrueReductors.OrderedInterreduce(mySPoly);  // ANNA
    if (myGRingInfo().myInputAndGrading()==AllGraded)
      myTrueReductors.interreduce(mySPoly);// this must be the first op
    ++myAgeValue;
    myPolys.push_back(mySPoly);// this must be the second op
    //myTrueReductors.SuperInterreduce(mySPoly);  // ANNA
    myTrueReductors.Insert(&myPolys.back());
    myGB.push_back(&myPolys.back());
    GPairList new_pairs;
    myBuildNewPairs(new_pairs);
    if (myCriteria.myBack)
      myApplyBCriterion();
    for_each(new_pairs.begin(), new_pairs.end(), mem_fun_ref(&GPair::myComplete));
    new_pairs.sort();
    myPairs.merge(new_pairs);
  }

  void GReductor::myUpdateBasisOnly()
  {
    // HERE NO MORE - in DOGBASIS or DOAHGBASIS If non homog, probably no interred
    //if (!IsTrueGCDDomain(CoeffRing(mySPoly)))
    //myTrueReductors.OrderedInterreduce(mySPoly);  // ANNA
    if (myGRingInfo().myInputAndGrading()==AllGraded)
      myTrueReductors.interreduce(mySPoly);// this has better to be the first op
    ++myAgeValue;
    myPolys.push_back(mySPoly);// this has better to be the second op
    myTrueReductors.Insert(&myPolys.back());
    myGB.push_back(&myPolys.back());
  }

/*
// This class exists only as a tool for the GReductor::MinPairs procedure
  template <class T1, class T2>
  struct _GPairDegBigger:public binary_function<T1,T2,bool>
  {
    bool operator()(const GPair& P, int Deg) const
    {
      return (P.IsInputPoly()||wdeg(P)>Deg);
    }
  };

  // stricly for use in the procdure below
  template <class T1, class T2>
  struct _GPairEqualLCMAndSecondIndex:public binary_function<T1,T2,bool>
  {
    bool operator()(const GPair& P1, const GPair& P2) const
      {
        return (P1.SecondIndex()==P2.SecondIndex()) && (LCM(P1)==LCM(P2));
      }
  };

  template <class T1, class T2>
  struct _GPairEqual:public binary_function<T1,T2,bool>
  {
    bool operator()(const GPair& P1, const GPair& P2) const
      {
        return (P1==P2);
      }
  };
*/

// // ANNA: apparently was used only by RingWeyl.C
//   void GReductor::myDoAFFGBasis()
//   {
//     CoCoA_ASSERT(myGetBuchbergerOpType() == AffineAlg);
//     // Input Polynomials sorted and zero polys deleted
//     if (myPolys.empty()) return;
//     if (myGRingInfoValue.AmIModule())
//       myCriteria.myCoprime = false; // CopCriterion works only for REAL ideals
//     double T=0.0;//STAT
//     GPoly Zero(myGRingInfoValue);
//     GPolyList::iterator it=myPolys.begin();
//     for (long i=0;it!=myPolys.end();++it,++i)
//       Ordered_Insert(myPairs,GPair(*it,i));
//     myOldDeg = wdeg(myPairs.front());//STAT
//     myCurrentPairDeg = myOldDeg;//STAT
//     if (myStat.myPrintDeg)
//       clog<<"\n[Log] ********** Starting_Pair_Degree="<<myOldDeg<<"\n";
//     while (!myPairs.empty()) {
//       myCurrentPairDeg=wdeg(myPairs.front());//STAT

//       if (myCurrentPairDeg!=myOldDeg)
//       {
//         if (myStat.myPrintDeg)
//           clog<<"\n[log] ********* Current_Pair_Degree_Now="<<myCurrentPairDeg<<"\n";
//         myStat.myUpgradeDegStats(myOldDeg,len(myPairs));
//         myOldDeg=myCurrentPairDeg;
//       };
//       if (!myPairs.empty())
//       {
//         if (myStat.myPrintNumPair)
//           clog<<"[log] pair="<<len(myPairs)<<" "<<flush;
//         if (myStat.myPrintReduction)
//           clog<<"doing="<<myPairs.front()<<" Len="
//               <<myStat.myPolyLens.back().first<<" reduction="<<flush;
//         if (myPairs.front().IamCoprime() && myCriteria.myCoprime)
//         {
//           mySPoly=Zero;
//           if (myStat.myPrintCop) clog<<"COP"<<flush;
//           ++myStat.myCopKilled;
//           //--myStat.myUseless;
//           T=0.0;
//         }
//         else
//         {
//           mySPoly.myAssignSPoly(myPairs.front(),myAgeValue);  // ??? SPoly computed only if not coprime
//           ++myAgeValue;
//           myStat.myPolyLens.push_back(make_pair(len(mySPoly),0));
//           T=CpuTime();
//           mySPoly.myReduce(myTrueReductors);
//           T-=CpuTime();
//           ++myStat.myNReductions;
//           myStat.myPolyLens.back().second=len(mySPoly);
//           if (IsZero(mySPoly))
//           {
//             if (myStat.myPrintReduction) clog << "0"<<flush;
//             ++myStat.myUseless;
//           }
//           else
//           {
//             if (myStat.myPrintReduction)
//               clog<<LPPForDiv(mySPoly)<<"+..<"<<len(myGB)-1<<"> Len="
//                   <<myStat.myPolyLens.back().second<<flush;
//             ++myStat.myUseful;
//             //myTrueReductors.interreduce(mySPoly);
//             if ((myStat.myPrintReduction && (len(mySPoly)<=2)))
//               clog<<"EASY REDUCTOR FOUND  LEN="<<len(mySPoly)<<" DEG="<<wdeg(mySPoly)<<endl;
//           }

//         }//else -if not coprime
//         myPairs.erase(myPairs.begin());// erase the used gpair
//         if (myStat.myPrintReduction)
//           clog << " time="<<-T<<" \n"<<flush;
//         if (!IsZero(mySPoly)) myUpdateBasisAndPairs();
//       };//end if minpairs
//     };//end main while
//     //STAT
//     if (myCurrentPairDeg!=myOldDeg)
//       myStat.myDegByDeg.push_back(DegStats(myOldDeg,0,0,0,0,0,1,0));
//     else myStat.myUpgradeDegStats(myOldDeg,0);
//     myStat.myStampa(clog);
//    }// end myDoAFFGBasis


   void GReductor::myReduceCurrentSPoly()
  {
    if (!myPrepared)
    {
      clog<<"GReductor: no preparation done ";
      return;
    }
    //if (myStat.myPrintDeg)
    //  clog<<"GReductor:myReduceCurrentSPoly performing one pair handling";
    if (myPolys.empty() || myPairs.empty())
      return;
    if (myStat.myPrintNumPair)
      clog<<"[log] pair num ="<<len(myPairs)<<" "<<flush;
    if (myStat.myPrintReduction)
        clog<<"doing="<<myPairs.front();
    if (myPairs.front().IamCoprime()&&myCriteria.myCoprime)
    {
      mySPoly = GPoly(myGRingInfoValue); // initialized as 0
      if (myStat.myPrintCop)
        clog<<"COP"<<flush;
      ++myStat.myCopKilled;
      myStat.myReductionTime=0.0;
    }
    else
    {
      mySPoly.myAssignSPoly(myPairs.front(), myAgeValue);  // ??? SPoly computed only if not coprime
      if (myPairs.front().IsInputPoly()) mySPoly.mySetMinimalGen();
      // ++myAgeValue; // moved into myUpdateBasisAndPairs/myUpdateBasisOnly
      myStat.myPolyLens.push_back(make_pair(len(mySPoly),0));
      if (myStat.myPrintReduction && IsZero(mySPoly))
        clog<<"SPOLY IS ZERO BEFORE REDUCTION !"<<endl;
      if (myStat.myPrintReduction)
        clog<<" SPolyLen="<<myStat.myPolyLens.back().first<<" reduction="<<flush;
      myStat.myReductionTime=CpuTime();
      mySPoly.myReduce(myTrueReductors);
      myStat.myReductionTime-=CpuTime();
      ++myStat.myNReductions;
      myStat.myPolyLens.back().second=len(mySPoly);
      if (IsZero(mySPoly))
      {
        if (myStat.myPrintReduction) clog << "0"<<flush;
        ++myStat.myUseless;
      }
      else
      {
        if (myStat.myPrintReduction)
          clog<<LPPForDiv(mySPoly)<<"+..<"<<len(myGB)<<"> Len="
              <<myStat.myPolyLens.back().second
              <<" Comp="<<Component(mySPoly)<<flush;
        ++myStat.myUseful;
        //if (!IsTrueGCDDomain(CoeffRing(mySPoly)) || len(mySPoly)<=2)
        //myTrueReductors.interreduce(mySPoly);
        if ((myStat.myPrintReduction && (len(mySPoly)<=2)))
          clog<<" EASY REDUCTOR FOUND  LEN="<<len(mySPoly)<<" DEG="<<wdeg(mySPoly)<<endl;
      }// else if (IsZero(mySPoly))
    }//else -if not coprime
    if (myStat.myPrintReduction)
      clog <<" time="<<-myStat.myReductionTime<<" \n"<<flush;
    myPairs.erase(myPairs.begin());// erase the used gpair
    //if (!IsZero(mySPoly)) myUpdateBasisAndPairs();
    if (!myPairs.empty())
    {
      myCurrentPairDeg=wdeg(myPairs.front());//STAT
      if (myCurrentPairDeg!=myOldDeg)
      {
        if (myStat.myPrintDeg)
        {
          clog<<
            "\n[log] ********* Current_Pair_Degree_Now="<<myCurrentPairDeg<<
            "\tSugar = " << sugar(myPairs.front()) << "\n";
        }
        myStat.myUpgradeDegStats(myOldDeg,len(myPairs));
        myOldDeg=myCurrentPairDeg;
      }//if (myCurrentPairDeg!=myOldDeg)
    }//if (!myPairs.empty())
    else myStat.myUpgradeDegStats(myCurrentPairDeg,0);
  }//myReduceCurrentSPoly


//   // ANNA: not called
//   void GReductor::myDoGBasis(const int ReductionNumber)
//   {
//     //    CoCoA_ASSERT(myGetBuchbergerOpType() == HomogeneousAlg);
//     //// now works also for Affine Algorithm
//     int RedNumber=ReductionNumber-1;
//     myPrepareGBasis();
//     while (!myPairs.empty()&&RedNumber!=0)
//     {
//       myReduceCurrentSPoly();
//       if (!IsZero(mySPoly))
//       {
//         myUpdateBasisAndPairs();
//         if (!myPairs.empty())
//           if (myCurrentPairDeg!=myOldDeg&&myTrueReductors.IhaveBorelReductors())
//             myTrueReductors.myBorelReductorsUpdateInNextDegree();
//       }
//       --RedNumber;
//     };
//   }//myDoGBasis(const int ReductionNumber)


  // ANNA: called by TmpGOperations
  void GReductor::myDoGBasis()
  {
    //    CoCoA_ASSERT(myGetBuchbergerOpType() == HomogeneousAlg);    
    myPrepareGBasis();
    while (!myPairs.empty())
    {
      myReduceCurrentSPoly();
      if (!IsZero(mySPoly))
      {
        myUpdateBasisAndPairs();
        if (!myPairs.empty())
          if (myCurrentPairDeg!=myOldDeg&&myTrueReductors.IhaveBorelReductors())
            myTrueReductors.myBorelReductorsUpdateInNextDegree();
      }
    }//while
    myFinalizeGBasis();
  }//myDoGBasis


  // ANNA: called by TmpGOperations
  void GReductor::myDoSATGBasis()
  {
    CoCoA_ASSERT(myGetBuchbergerOpType() == SaturatingAlg);
    const long HVarIndex=NumIndets(myGRingInfoValue.myNewSPR())-1;// This is OK for Ideals only
    bool satjumped=false;
    degree SPolyPredDeg(GradingDim(myGRingInfoValue.myNewSPR()));// Used for stats in the dehmog alg
    myPrepareGBasis();
    //    double T=0.0;
    while (!myPairs.empty())
    {
      myReduceCurrentSPoly();
      if (!IsZero(mySPoly))
      {  
      	SPolyPredDeg = wdeg(mySPoly);
	//mySPoly.smart_dehomog_DRL(HVarIndex);
	mySPoly.smart_dehomog(HVarIndex);
      	if (SPolyPredDeg!=wdeg(mySPoly))
      	{
          satjumped=true;
          ++myStat.myPolyDHed;
          myStat.myDegDH += ConvertTo<long>((SPolyPredDeg-wdeg(mySPoly))[0]);
          if (myStat.myPrintPolyDH & false)
            clog << "\nPOLY DEHOMOG SUCCESSFULLY    Deg "
                 <<SPolyPredDeg<< "-->"<<wdeg(mySPoly)<<" "
                 <<LPPForDiv(mySPoly)<<endl;
        }//if
        myUpdateBasisAndPairs();
 	if (!myPairs.empty())
          if (myCurrentPairDeg!=myOldDeg&&myTrueReductors.IhaveBorelReductors())
            myTrueReductors.myBorelReductorsUpdateInNextDegree();
      }
    }//while
    myFinalizeGBasis();
  }//myDoSATGBasis
  

  void GReductor::myDoGBasisByBatch()
  {
    bool NonZeroReductionPerformed=false;
    myPrepareGBasis();

    // process the special pairs
    while (!myPairs.empty())
    {
      myReduceCurrentSPoly();
      if (!IsZero(mySPoly))
        myUpdateBasisOnly();
    }//while

    do
    {
      NonZeroReductionPerformed=false;
      myCreatePairs();
//myStampaGB(clog);
    while (!myPairs.empty())
      {
        myReduceCurrentSPoly();
        if (!IsZero(mySPoly))
        {
          myUpdateBasisOnly();
	  NonZeroReductionPerformed=true;
//clog << "myDoGBasis:************************************** NON ZERO FOUND***************** "<<endl;
        }
      }//while
    }
    while (NonZeroReductionPerformed);
     myFinalizeGBasis();
   }//myDoGBasisByBatch

  void GReductor::myReduceUntilNonZeroRedSP()
  {
    bool flag=true;
    while ((!myPairs.empty())&&(IsZero(mySPoly)||flag))
    {
      flag=false;
      myReduceCurrentSPoly();
      if (!IsZero(mySPoly))
      {
        myUpdateBasisAndPairs();
 	if (!myPairs.empty())
          if (myCurrentPairDeg!=myOldDeg&&myTrueReductors.IhaveBorelReductors())
            myTrueReductors.myBorelReductorsUpdateInNextDegree();
      }
    }
  }//myReduceUntilNonZeroRedSP


  void GReductor::myCreatePairs()
  {
clog << "myCreatePairs:begin "<<endl;
      GPairList new_pairs;
      myBuildNewPairsAll(new_pairs);
      if (myCriteria.myBack)
        myApplyBCriterion();
      for_each(new_pairs.begin(), new_pairs.end(), mem_fun_ref(&GPair::myComplete));
      new_pairs.sort();
      swap(myPairs,new_pairs);
clog << "myCreatePairs:end "<<endl;
  }//myCreatePairs


  void GReductor::myCreateInputPolyPairs(GPolyList& theCandidateBasis)
  {
    GPolyList::iterator it=theCandidateBasis.begin();
    long i=0;
    for (;it!=theCandidateBasis.end();it++,i++)
      Ordered_Insert(myPairs,GPair(*it,i));
  }//myCreateInputPolyPairs


  // Prepare the first (special) pairs
  void GReductor::myPrepareGBasis()
  {
    for (GPolyList::iterator it=myPolys.begin(); it!=myPolys.end(); ++it)
      it->myInitializeSugar(myGRingInfoValue.myNewSugar(poly(*it)));
    myCreateInputPolyPairs(myPolys);
    myPrepared=true;
    if (myGRingInfoValue.AmIModule())
      myCriteria.myCoprime = false;// CopCriterion works only for REAL ideals
    CoCoA_ASSERT(len(myPairs)!=0);
    myOldDeg = wdeg(myPairs.front());//STAT
    myCurrentPairDeg = myOldDeg;//STAT
    myStat.myReductionTime=0.0;//STAT
    myStat.myTotalTime=CpuTime();//STAT
    if (myStat.myPrintDeg)
      clog<<"\n[log] ********* Starting_Pair_Degree="<<myCurrentPairDeg<<"\n";
   }//myPrepareGBasis


  // Prepare the first (special) pairs
  void GReductor::myPrepareGBasisPairsExcluded()
  {
    myPrepared=true;
    if (myGRingInfoValue.AmIModule())
      myCriteria.myCoprime = false;// CopCriterion works only for REAL ideals
    myOldDeg = wdeg(myPairs.front());//STAT
    myCurrentPairDeg = myOldDeg;//STAT
    myStat.myReductionTime=0.0;//STAT
    myStat.myTotalTime=CpuTime();//STAT
    if (myStat.myPrintDeg)
      clog<<"\n[log] ********* Starting_Pair_Degree="<<myCurrentPairDeg<<"\n";
   }//myPrepareGBasis


  void GReductor::myFinalizeGBasis()
  {
    myStat.myTotalTime-=CpuTime();//STAT
    // interreduction if the ordering is not degree compatible
    if (true)
      if (myGRingInfo().myInputAndGrading()!=AllGraded) // special case: homogeneous wrt 1st row
      {
        for (GPolyPtrList::iterator it=myGB.begin();it!=myGB.end();++it)
        {
          vector<ReductorData>::iterator it1=myTrueReductors.find(*it);
          it1->mySetIamNotToBeUsed(true);
          (**it).myReduce(myTrueReductors);
          if (IsZero(**it))    //  remove the zeros
          {
            it=myGB.erase(it);
            if (it!=myGB.begin())
              --it;
          }
          else
            it1->mySetIamNotToBeUsed(false);
        }//myGB for
      }//interreduction
    myStat.myStampa(clog);
  }//myFinalizeGBasis


//   // ANNA: not called
//   void GReductor::myDoGBasisTEST()
//   {
//     double T=0.0;
// //clog<<"myDoGBasisTEST:begin"<<endl;
//     // Input Polynomials sorted and zero polys deleted
//     if (myPolys.empty()) return;
//     if (myGRingInfoValue.AmIModule())
//       myCriteria.myCoprime = false;// CopCriterion works only for REAL ideals
//     GPoly Zero(myGRingInfoValue);
//     GPolyList::iterator it=myPolys.begin();
//     for (long i=0;it!=myPolys.end();++it,++i)
//       Ordered_Insert(myPairs,GPair(*it,i));
//     myOldDeg = wdeg(myPairs.front());
//     myCurrentPairDeg = myOldDeg;
//     if (myStat.myPrintDeg)
//       clog<<"\n[log] ********* Starting_Pair_Degree="<<myCurrentPairDeg<<"\n";
//     while (!myPairs.empty())
//     {
//       if (myStat.myPrintNumPair)
//         clog<<"[log] pair="<<len(myPairs)<<" "<<flush;
//       if (myStat.myPrintReduction)
//         clog<<"doing="<<myPairs.front();
//       if (myPairs.front().IamCoprime()&& myCriteria.myCoprime)
//       {
//         mySPoly=Zero;
//         if (myStat.myPrintCop)
//           clog<<"COP"<<flush;
//         ++myStat.myCopKilled;
//         T=0.0;
//        }
//        else
//        {
//          mySPoly.myAssignSPoly(myPairs.front(),myAgeValue);  // ??? SPoly computed only if not coprime
//          ++myAgeValue;
//          myStat.myPolyLens.push_back(make_pair(len(mySPoly),0));
//          if (myStat.myPrintReduction && IsZero(mySPoly))
//             clog<<"SPOLY IS ZERO BEFORE REDUCTION !"<<endl;
//          if (myStat.myPrintReduction)
//             clog<<" SPolyLen="<<myStat.myPolyLens.back().first<<" reduction="<<flush;
//          myPairs.erase(myPairs.begin());// erase the used gpair
//          T=CpuTime();
//          mySPoly.myReduce(myTrueReductors);
//          ++myStat.myNReductions;
//          T-=CpuTime();
//          myStat.myPolyLens.back().second=len(mySPoly);
//          if (IsZero(mySPoly))
//          {
//            if (myStat.myPrintReduction) clog << "0"<<flush;
//            ++myStat.myUseless;
//          }
//          else
//          {
//            if (myStat.myPrintReduction)
//              clog<<LPPForDiv(mySPoly)<<"+..<"<<len(myGB)+1<<"> Len="
//                  <<myStat.myPolyLens.back().second
//                  <<" Comp="<<Component(mySPoly)<<flush;
//            ++myStat.myUseful;
//            //if (!IsTrueGCDDomain(CoeffRing(mySPoly)) || len(mySPoly)<=2)
//            //myTrueReductors.interreduce(mySPoly);
//            if ((myStat.myPrintReduction && (len(mySPoly)<=2)))
//              clog<<" EASY REDUCTOR FOUND  LEN="<<len(mySPoly)<<" DEG="<<wdeg(mySPoly)<<endl;
//          }// else
//        }//else -if not coprime
//        if (myStat.myPrintReduction)
//           clog <<" time="<<-T<<" \n"<<flush;
//        if (!IsZero(mySPoly)) myUpdateBasisAndPairs();
//        if (!myPairs.empty())
//        {
//          myCurrentPairDeg=wdeg(myPairs.front());//STAT
//          if (myCurrentPairDeg!=myOldDeg)
//          {
//            if (myStat.myPrintDeg)
//              clog<<"\n[log] ********* Current_Pair_Degree_Now="<<myCurrentPairDeg<<"\n";
//            if (myTrueReductors.IhaveBorelReductors())
//              myTrueReductors.myBorelReductorsUpdateInNextDegree();
//            myStat.myUpgradeDegStats(myOldDeg,len(myPairs));
//            myOldDeg=myCurrentPairDeg;
//          }
//        }
//        else
//          myStat.myUpgradeDegStats(myCurrentPairDeg,0);
//     };//end main while
//     myStat.myTotalTime-=CpuTime();
//     myStat.myStampa(clog);
//   }// end _myDoGBasis



//At the moment, no one is calling this procedure. The caller calls myDoGBasisTEST instead
  void GReductor::myDoSATMixGBasis()
  {
    CoCoA_ERROR(ERR::SERIOUS, "myDoSATMixGBasis: should not be called!");
    // MOD the homog var and the module vars conflict
    // Input Polynomials sorted and zero polys deleted
    if (myPolys.empty()) return;

    if (myGRingInfoValue.AmIModule())
      myCriteria.myCoprime = false;// CopCriterion works only for REAL ideals

    CoCoA_ASSERT(GradingDim(myGRingInfoValue.myNewSPR()) == 1);
    const long HVarIndex=NumIndets(myGRingInfoValue.myNewSPR())-1;// This is OK for Ideals only
    bool satjumped=false;
    double T=0.0;//STAT
    degree SPolyPredDeg(GradingDim(myGRingInfoValue.myNewSPR()));// Used for stats in the dehmog alg

    GPoly Zero(myGRingInfoValue);


    GPolyList::iterator it=myPolys.begin();
    for (long i=0;it!=myPolys.end();++it,++i)
    {
      Ordered_Insert(myPairs,GPair(*it,i));
    }

    myOldDeg = wdeg(myPairs.front());//STAT
    myCurrentPairDeg = myOldDeg;//STAT

    if (myStat.myPrintDeg)
      clog<<"\n[log] ********* Starting_Pair_Degree="<<myOldDeg<<"\n";


    while (!myPairs.empty())
    {
      satjumped=false;
      myCurrentPairDeg=wdeg(myPairs.front());//STAT

      if (myCurrentPairDeg!=myOldDeg)
      {
        if (myStat.myPrintDeg)
          clog<<"\n[log] New_Degree="<<myCurrentPairDeg<<"\n";
        //      myTrueReductors.AdjustLimit(myCurrentPairDeg);
        myStat.myUpgradeDegStats(myOldDeg,len(myPairs));
        myOldDeg=myCurrentPairDeg;
      }

      if (!myPairs.empty())
      {
        // WARN - POSSIBLY DIFFERENT
        if (myStat.myPrintNumPair)
          clog<<"[log] pair="<<len(myPairs)<<" "<<flush;
        if (myStat.myPrintReduction)
          clog<<"doing="<<myPairs.front()<<" Len="
              <<myStat.myPolyLens.back().first<<"-->"<<flush;
        if (myPairs.front().IamCoprime()&& myCriteria.myCoprime)
        {
          mySPoly=Zero;
          if (myStat.myPrintCop) clog<<"COP"<<flush;
          ++myStat.myCopKilled;
          //--myStat.myUseless;
          T=0.0;
        }
        else
        {
          mySPoly.myAssignSPoly(myPairs.front(), ++myAgeValue);  // ??? SPoly computed only if not coprime
          ++myAgeValue;
          myStat.myPolyLens.push_back(make_pair(len(mySPoly),0));
          T=CpuTime();
          //mySPoly.smart_dehomog(HVarIndex);
          mySPoly.myReduce(myTrueReductors);
          ++myStat.myNReductions;
          T-=CpuTime();

          myStat.myPolyLens.back().second=len(mySPoly);
          if (IsZero(mySPoly))
          {
            if (myStat.myPrintReduction) clog << "0"<<flush;
            ++myStat.myUseless;
          }
          else
          {
            if (myStat.myPrintReduction)
              clog<<LPPForDiv(mySPoly)<<"+..<"<<len(myGB)+1<<"> Len="
                  <<myStat.myPolyLens.back().second<<flush;
            SPolyPredDeg = wdeg(mySPoly);
            mySPoly.smart_dehomog(HVarIndex);

            if (SPolyPredDeg!=wdeg(mySPoly))
            {
              satjumped=true;
              ++myStat.myPolyDHed;
              myStat.myDegDH += ConvertTo<long>((SPolyPredDeg-wdeg(mySPoly))[0]);
              if (myStat.myPrintPolyDH & false)
                clog << "\nPOLY DEHOMOG SUCCESSFULLY    Deg "
                     << SPolyPredDeg << "-->" << wdeg(mySPoly)
                     << LPPForDiv(mySPoly) << endl;
            }
            ++myStat.myUseful;
            //if (!IsTrueGCDDomain(CoeffRing(mySPoly)) && (len(mySPoly)<=2))
            //if (len(mySPoly)<=2)
            //if (!satjumped)
            //myTrueReductors.SuperInterreduce(mySPoly);
            //myTrueReductors.interreduce(mySPoly);
            if ((myStat.myPrintReduction && (len(mySPoly)<=2)))
              //if (len(mySPoly)<=2)
              clog<<"EASY REDUCTOR FOUND  LEN="<<len(mySPoly)<<" DEG="<<wdeg(mySPoly)<<endl;
          }// else
        }//else -if not coprime
        myPairs.erase(myPairs.begin());// erase the used gpair
        if (myStat.myPrintReduction)
          clog <<" time="<<-T<<" \n"<<flush;
        if (!IsZero(mySPoly)) myUpdateBasisAndPairs();
      };//end if minpairs
    };//end main while
    //STAT
    if (myCurrentPairDeg!=myOldDeg)
      myStat.myDegByDeg.push_back(DegStats(myOldDeg,0,0,0,0,0,1,0));
    else myStat.myUpgradeDegStats(myOldDeg,0);
    if (myStat.myGetLevel()>-1)
    {
      clog << "number of reductions = " << myStat.myNReductions << endl;
      clog << "# GBasis = " << myStat.myUseful-myStat.myPolyDeleted<< endl;
    }
     myStat.myStampa(clog);
  }// end myDoSATMixGBasis


  // ANNA: not called
  void GReductor::Rebuild(const PolyList& thePL)
  {
    clog << "Rebuild:begin " <<endl;
    mySPoly=GPoly(myGRingInfoValue);
    clog << "Rebuild:before for " <<endl;
    myCurrentPairDeg=degree(GradingDim(myGRingInfoValue.myNewSPR()));
    myOldDeg=degree(GradingDim(myGRingInfoValue.myNewSPR()));
    //myStat(len(thePL),myStat.myGetLevel());
    myPrepared=false;
    myAgeValue = 0;
    myStat.myNReductions=0;
    //IsDynamicAlgorithm unchanged
    myWrongLPPFoundValue=false;
    //myCriteria unchanged
    myPairs.clear();
    myTrueReductors.myClear();
    myGB.clear();
    myPolys.clear();

    for (PolyList::const_iterator it=thePL.begin();it!=thePL.end();it++)
    {
      if (!IsZero(*it))
        myPolys.push_back(GPoly(*it,myGRingInfoValue));
    }
    myPolys.sort(BoolCmpLPPGPoly);
clog << "Rebuild:end " <<endl;
  }//Rebuild


//////////////////////////// Embedding/DeEmbedding /////////////////////


ModOrdTypeForcing ModuleOrderType(const FreeModule& M)
{
  if (IsOrdPosn(ordering(M))) return WDegTOPos;
  if (IsWDegPosnOrd(ordering(M))) return WDegPosTO;
  return PosWDegTO;
}//ModOrdType



  FreeModule owner(const VectorList& theVL)
  {
    CoCoA_ASSERT(!theVL.empty());
    return owner(*theVL.begin());
  }//owner



// similar functions now moved into FreeModule.C
// FreeModule MakeNewFreeModuleForSyz2(const GPolyList& theGPL)
// {
//   if (theGPL.empty())
//     return NewFreeModule(RingQQ(),1);
//   const SparsePolyRing OldP=owner(theGPL);
//   std::vector<degree> InputShifts;
//   for (GPolyList::const_iterator it=theGPL.begin(); it!=theGPL.end(); ++it)
//      InputShifts.push_back(wdeg(*it));

//   return NewFreeModule(OldP, NewWDegPosnOrd(ordering(PPM(OldP)), InputShifts));
// }//MakeNewFreeModuleForSyz


// returns the poly ring equivalent with OldP^2, same grading
// The ordering is WDegPosnOrd if MOType==NoForcing or MOType
SparsePolyRing MakeNewPRingForSimpleEmbeddingPosFirst(const SparsePolyRing& OldP,
                                                      bool HomogInput)
{
  if (HomogInput)
    return MakeNewPRingForSimpleEmbedding(OldP, WDegPosTO);
  else
    return MakeNewPRingForSimpleEmbedding(OldP,PosWDegTO);
}					
SparsePolyRing MakeNewPRingForSimpleEmbedding(const SparsePolyRing& OldP,
                                              ModOrdTypeForcing MOType)
{
  std::vector<degree> InputShifts;
  {
    degree tmp(GradingDim(OldP));
    InputShifts.push_back(tmp);
    InputShifts.push_back(tmp);
  }
  const FreeModule FM=NewFreeModule(OldP, InputShifts, WDegPosnOrd);
  return MakeNewPRingFromModule(FM, MOType);
}//MakeNewPRingForSimpleEmbedding

SparsePolyRing MakeNewPRingForSimpleEmbedding(const SparsePolyRing& OldP)
{
  return MakeNewPRingForSimpleEmbedding(OldP,NoForcing);
}//MakeNewPRingForSimpleEmbedding

SparsePolyRing MakeNewPRingFromModule(const FreeModule& FM)
{
  return MakeNewPRingFromModule(FM,NoForcing);
}//MakeNewPRingFromModule

SparsePolyRing MakeNewPRingFromModulePosFirst(const FreeModule& FM,
					      bool HomogInput)
{
  if (HomogInput)
    return MakeNewPRingFromModule(FM, WDegPosTO);
  else
    return MakeNewPRingFromModule(FM, PosWDegTO);
}					

// This is OK for the non-homogeneous case
// For the homogenous case, PosTo this is inefficient, since
// the Deg rows in the To part are useless.
SparsePolyRing MakeNewPRingFromModule(const FreeModule& FM,
                                      ModOrdTypeForcing MOType)
  {
   const ModuleOrdering MTO = ordering(FM);
   const SparsePolyRing OldP=RingOf(FM);
   const long NumOldInds=NumIndets(OldP);
   long GrDim;
   if (MOType==PosWDegTO)
     GrDim=0;// Set simple sugar on
   else
     GrDim=GradingDim(OldP);
   const long NumNewInds=NumOldInds+GrDim+1;

   matrix OldOrdOMat(OrdMat(OldP));

   matrix NewOrdMat(NewDenseMat(RingQQ(), NumNewInds, NumNewInds));

   if (MOType == NoForcing)  MOType = ModuleOrderType(FM);

   switch (MOType)
   {
     case PosWDegTO:
       // Setting the module component ordering
       SetEntry(NewOrdMat, 0, NumNewInds-1, 1); 	
      // Part common to IsWDegPosnOrd and IsOrdPosn
       // Setting the Grading: the OldGrading		
       for (long i=0; i < GrDim+1; ++i)			
         for (long j=0; j < NumOldInds; ++j)		
            SetEntry(NewOrdMat, i+1, j, OldOrdOMat(i,j));
       // Setting the Grading: the NewGrading		
       for (long i=1; i < GrDim+1; ++i)			
         SetEntry(NewOrdMat, i, i+NumOldInds-1, 1); 	
      // Setting the TO ordering
       for (long i=GrDim; i < NumOldInds; ++i)
         for (long j=0; j < NumOldInds; ++j)
           SetEntry(NewOrdMat, i+1, j, OldOrdOMat(i,j));
        break;
	
     case WDegTOPos:
       // Part common to IsWDegPosnOrd and IsOrdPosn
//clog<<"MakeNewPRingFromModule:case OrdPosn"<<endl;
      // Setting the Grading: the OldGrading		
       for (long i=0; i < GrDim; ++i)			
         for (long j=0; j < NumOldInds; ++j)		
            SetEntry(NewOrdMat, i, j, OldOrdOMat(i,j));
       // Setting the Grading: the NewGrading		
       for (long i=0; i < GrDim; ++i)			
         SetEntry(NewOrdMat, i, i+NumOldInds, 1);
//clog<<"MakeNewPRingFromModule:the matrix graded "<<NewOrdMat<<endl;
      // Setting the TO	
       for (long i=GrDim; i < NumOldInds; ++i)
         for (long j=0; j < NumOldInds; ++j)
           SetEntry(NewOrdMat, i, j, OldOrdOMat(i,j));
//clog<<"MakeNewPRingFromModule:the matrix TO "<<NewOrdMat<<endl;
      // Setting the module component ordering
       SetEntry(NewOrdMat, NumNewInds-1-GrDim, NumNewInds-1, 1); 	
//clog<<"MakeNewPRingFromModule:the matrix TO Pos "<<NewOrdMat<<endl;
        break;

     case WDegPosTO:; // This is the default
     default:
       // Part common to IsWDegPosnOrd and IsOrdPosn
       // Setting the Grading: the OldGrading		
       for (long i=0; i < GrDim; ++i)			
         for (long j=0; j < NumOldInds; ++j)		
            SetEntry(NewOrdMat, i, j, OldOrdOMat(i,j));
       // Setting the Grading: the NewGrading		
       for (long i=0; i < GrDim; ++i)			
         SetEntry(NewOrdMat, i, i+NumOldInds, 1); 	
       // Setting the module component ordering
       SetEntry(NewOrdMat, GrDim, NumNewInds-1, 1);
       // Setting the TO ordering
       for (long i=GrDim; i < NumOldInds; ++i)
         for (long j=0; j < NumOldInds; ++j)
           SetEntry(NewOrdMat, i+1, j, OldOrdOMat(i,j));
         break;
   }
   // Filling the matrix
   for (long i=0; i < GrDim; ++i)
     for (long j=0; j < NumOldInds; ++j)
        SetEntry(NewOrdMat, NumOldInds+i+1, j, OldOrdOMat(i,j));

//clog<<"MakeNewPRingFromModule:the matrix"<<NewOrdMat<<endl;

   const PPOrdering MatNewOrd = NewMatrixOrdering(NumNewInds,
  					          GrDim,
  					          NewOrdMat);
						
  const std::vector<symbol> IndetNames = NewSymbols(NumOldInds + GrDim + 1);
//---> for DEBUGGING choose these IndetNames:
//   std::vector<symbol> IndetNames = SymbolRange("x", 0, NumOldInds-1);
//   if ( GrDim==1 ) IndetNames.push_back(symbol("s"));  // indet representing shift
//   else
//     for ( long i=0 ; i<GrDim ; ++i )
//       IndetNames.push_back(symbol("s",i));  // indet representing shift
//   IndetNames.push_back(symbol("e"));  // indet representing module component
//---> for DEBUGGING
  SparsePolyRing NewP(NewSparsePolyRing(CoeffRing(OldP),IndetNames,MatNewOrd));
  return NewP;
 }//MakeNewPRingFromModule


  GPoly EmbedPoly(ConstRefRingElem p,
    		  const GRingInfo& theGRI,
		  const long CompIndex)
  {
    //std::vector<RingElem> images;
    //for (long i=0;i!=NumIndets(theGRI.myOldSPR());++i)
    //  images.push_back(indet(theGRI.myNewSPR(), i));
    const RingHom& phi=theGRI.myOldP2NewP();
    return GPoly(phi(p)*theGRI.myE(CompIndex), theGRI);
  }//EmbedPoly

  GPoly EmbedPoly(ConstRefRingElem the_p,
                  const GRingInfo& theGRI,
                  const degree& the_d,
		  const long CompIndex)

   {
     const RingHom& phi=theGRI.myOldP2NewP();
     return GPoly(phi(the_p)*theGRI.myE(CompIndex)*theGRI.myY(the_d), theGRI);
   }//EmbedPoly


     // Embed a true polylist in a GPolyList. Component CompIndex
  GPolyList EmbedPolyList(const PolyList& thePL,
                          const GRingInfo& theGRI,
                          const long CompIndex)
  {
    GPolyList result;
    if (thePL.empty())
      return result;
    for (PolyList::const_iterator it=thePL.begin();it!=thePL.end();++it)
      result.push_back(EmbedPoly(*it, theGRI, CompIndex));
    return result;
  }//EmbedPolyList


 // Embed a true polylist in a GPolyList. Component CompIndex
  GPolyList EmbedPolyList(const PolyList& thePL,
                          const GRingInfo& theGRI,
                          const degree& the_d,
                          const long CompIndex)
  {
    GPolyList result;
    if (thePL.empty())
      return result;
    for (PolyList::const_iterator it=thePL.begin();it!=thePL.end();++it)
      result.push_back(EmbedPoly(*it, theGRI, the_d, CompIndex));
   return result;
  }//EmbedPolyList


/*
This realizes the embedding FM-->NewP of a vector v.
e_i->EY[i] and OldP2NewP gives the RingHom between BaseRing(FM) and NewP
Should be changed to avoid passing FM,NewP.
Is is here only for completeness/debug purposes.
*/
GPoly EmbedVector(const ModuleElem& v,
                  const GRingInfo& theGRI)			
{
  return EmbedVector(v, theGRI, 0);
}//EmbedVector


/*
This realizes the embedding FM-->NewP of a vector v.
e_i->EY[i] and OldP2NewP gives the RingHom between BaseRing(FM) and NewP
Should be changed to avoid passing FM,NewP.
Is is here only for completeness/debug purposes.
*/
GPoly EmbedVector(const ModuleElem& v,
                  const GRingInfo& theGRI,
                  const long StartingFromCompIndex)			
{
  RingElem p(theGRI.myNewSPR()),eMax(power(theGRI.myE(),StartingFromCompIndex));
  const RingHom& phi=theGRI.myOldP2NewP();
  for (long i=0; i<NumCompts(owner(v)); ++i)
    p+=phi(v[i])*theGRI.myEY(i)/eMax;
  return GPoly(p, theGRI);
}//EmbedVector



  GPolyList EmbedVectorList(const VectorList& InputVectorList,
                            const GRingInfo& theGRI)

  {
    return EmbedVectorList(InputVectorList, theGRI, 0);
  }//EmbedVectorList

  GPolyList EmbedVectorList(const VectorList& theVL,
                            const GRingInfo& theGRI,
			    const long StartingFromCompIndex)
 {
   GPolyList result;
   if (theVL.empty())
     return result;
   for (std::vector<ModuleElem>::const_iterator it=theVL.begin();it!=theVL.end();++it)
     result.push_back(EmbedVector(*it, theGRI, StartingFromCompIndex));
   return result;
}//EmbedVectorList		


    GPolyList SyzEmbedVectorList(const VectorList& InputVectorList,
                                 const GRingInfo& theGRI)
  {
    GPolyList result;
    if (InputVectorList.empty())
      return result;
    const SparsePolyRing NewP=theGRI.myNewSPR();
    result=EmbedVectorList(InputVectorList, theGRI);
    long k=NumCompts(theGRI.myFreeModule());
    if (theGRI.myInputAndGrading()==AffineInputWithGradedRing)
      for (GPolyList::iterator it=result.begin();it!=result.end();++it,++k)
      { // Added by JAA 2012/10/11
        RingElem Ek = theGRI.myE(k);
        (*it).myAppendClear(Ek);
      }
    else
      for (GPolyList::iterator it=result.begin();it!=result.end();++it,++k)
      { // Added by JAA 2012/10/11
        RingElem EkY = theGRI.myE(k)*theGRI.myY(wdeg(*it));
        (*it).myAppendClear(EkY);
      }

    return result;
  }//SyzEmbedVectorList

  GPolyList SyzEmbedPolyList(const PolyList& InputPolyList,
                             const GRingInfo& theGRI)
  {
    GPolyList result;
    if (InputPolyList.empty())
      return result;
    const SparsePolyRing NewP=theGRI.myNewSPR();
    result=EmbedPolyList(InputPolyList, theGRI, 0);
    RingElem SyzPP(NewP); // Gives the right degree to p+E^i
    degree d(GradingDim(NewP));
    long k=1;
    for (GPolyList::iterator it=result.begin();it!=result.end();++it,++k)
    {
      SyzPP=one(NewP);
      d=wdeg(*it);
      for (long j=0; j < GradingDim(NewP); ++j)
        SyzPP*=power(theGRI.myY(j),d[j]);
      if (theGRI.myInputAndGrading()==AffineInputWithGradedRing)
      { // Added by JAA 2012/10/11
        RingElem Ek = theGRI.myE(k);
        (*it).myAppendClear(Ek);
      }
      else
      { // Added by JAA 2012/10/11
        RingElem EkSyzPP = theGRI.myE(k)*SyzPP;
        (*it).myAppendClear(EkSyzPP);
      }
    }
    return result;
  }//SyzEmbedPolyList




  GPolyList IntEmbedPolyLists(const PolyList& thePL1,
                              const PolyList& thePL2,
                              const GRingInfo& theGRI)
  {
    GPolyList FirstPart=EmbedPolyList(thePL1, theGRI, 0);
    GPolyList SecondPart=EmbedPolyList(thePL2, theGRI, 0);
    GPolyList ThirdPart=EmbedPolyList(thePL2, theGRI, 1);
     GPolyList::iterator it1=ThirdPart.begin();
     for (GPolyList::iterator it=SecondPart.begin();it!=SecondPart.end();++it,++it1)
       (*it).myAppendClear(*it1);
     SecondPart.splice(SecondPart.begin(),FirstPart);
     return SecondPart;
  }//IntEmbedPolyLists


  GPolyList IntEmbedVectorLists(const VectorList& theVL1,
                                const VectorList& theVL2,
                                const GRingInfo& theGRI)
  {
    const long NC = NumCompts(owner(theVL1));
    //const SparsePolyRing NewP=theGRI.myNewSPR();
    GPolyList FirstPart=EmbedVectorList(theVL1, theGRI);
    GPolyList SecondPart=EmbedVectorList(theVL2, theGRI);
    GPolyList ThirdPart=EmbedVectorList(theVL2, theGRI, NC);
    GPolyList::iterator it1=ThirdPart.begin();
    for (GPolyList::iterator it=SecondPart.begin();it!=SecondPart.end();++it,++it1)
      (*it).myAppendClear(*it1);
    SecondPart.splice(SecondPart.begin(),FirstPart);
    return SecondPart;
  }//IntEmbedVectorLists

  // VL2 is a singleton
  GPolyList ColonEmbedVectorLists(const VectorList& theVL1,
                                  const VectorList& theVL2,
                                  const GRingInfo& theGRI)
  {
    GPolyList FirstPart;
    if (theVL1.empty())
      return FirstPart;
    const long NC = NumCompts(owner(theVL1));
    // const SparsePolyRing NewP=theGRI.myNewSPR();
    FirstPart=EmbedVectorList(theVL1, theGRI);
    GPoly GP1=EmbedVector(theVL2.front(), theGRI);
    degree d=wdeg(theVL2.front());
    GPoly GP2=EmbedPoly(one(theGRI.myOldSPR()), theGRI, d, NC);
    GP1.myAppendClear(GP2);
    FirstPart.push_back(GP1);
    return FirstPart;
  }//ColonEmbedVectorLists


  // PL2 is a singleton
  GPolyList ColonEmbedPolyLists(const PolyList& thePL1,
                                const PolyList& thePL2,
                                const GRingInfo& theGRI)
  {
    //std::clog<<"ComputeColon: start "<<endl;
     GPolyList FirstPart;
     if (thePL1.empty())
       return FirstPart;
     const SparsePolyRing NewP=theGRI.myNewSPR();
     FirstPart=EmbedPolyList(thePL1, theGRI, 0);
     //std::clog<<"ComputeColon: FirstPart "<<FirstPart<<endl;
     GPoly GP1=EmbedPoly(thePL2.front(), theGRI, 0);
     //std::clog<<"ComputeColon: GP1 "<<GP1<<endl;
     degree d=wdeg(thePL2.front());
     GPoly GP2=EmbedPoly(one(theGRI.myOldSPR()), theGRI, d, 1);
     GP1.myAppendClear(GP2);
     FirstPart.push_back(GP1);
     //std::clog<<"ComputeColon: end "<<endl;
     return FirstPart;
  }//ColonEmbedPolyLists


  ModuleElem DeEmbedPoly(ConstRefRingElem p,
                         const GRingInfo& theGRI)
  {
     return DeEmbedPoly(p,theGRI,0);
  }//DeEmbedPoly


  ModuleElem DeEmbedPoly(ConstRefRingElem p,
                         const GRingInfo& theGRI,
                         const long ComponentsLimit) // the component in p that goes to the 0 component of the output vector v. Lesser components of p go to higher component of v
  {
//std::clog<<"DeEmbedPoly:start"<<endl;
//std::clog<<"DeEmbedPoly:p  "<<p<<endl;
//std::clog<<"DeEmbedPoly:theGRI  "<<theGRI<<endl;
   const SparsePolyRing OldP=theGRI.myOldSPR();
   const SparsePolyRing NewP=theGRI.myNewSPR();
   const FreeModule FM=theGRI.myOutputFreeModule();
   ModuleElem v(FM);

//std::clog<<"DeEmbedPoly:v  "<<v<<endl;
//std::clog<<"DeEmbedPoly:theGRI.myNewP2OldP()  "<<theGRI.myNewP2OldP()<<endl;
//std::clog<<"DeEmbedPoly:theGRI.myOldP2NewP()  "<<theGRI.myOldP2NewP()<<endl;

   const std::vector<ModuleElem>& e = gens(FM);

   RingElem tmp(OldP);
   for (SparsePolyIter i=BeginIter(p); !IsEnded(i); ++i)
   {
//std::clog<<"DeEmbedPoly:ComponentsLimit  "<<ComponentsLimit<<endl;
//std::clog<<"DeEmbedPoly:PP(i)  "<<PP(i)<<endl;
      tmp=theGRI.myNewP2OldP()(monomial(NewP,coeff(i),PP(i)));
      CoCoA_ASSERT(theGRI.myPhonyComponent(PP(i))-ComponentsLimit < NumCompts(FM));
//std::clog<<"DeEmbedPoly:monomial(NewP,coeff(i),PP(i))  "<<monomial(NewP,coeff(i),PP(i))<<endl;
//std::clog<<"DeEmbedPoly:tmp  "<<tmp<<endl;
//std::clog<<"DeEmbedPoly:theGRI.myPhonyComponent(PP(i))-ComponentsLimit  "<<theGRI.myPhonyComponent(PP(i))-ComponentsLimit<<endl;
//std::clog<<"DeEmbedPoly:theGRI.myPhonyComponent(PP(i)) "<<theGRI.myPhonyComponent(PP(i))<<endl;
//std::clog<<"DeEmbedPoly:theGRI.myComponent(PP(i)) "<<theGRI.myComponent(PP(i))<<endl;
      v+=tmp*e[theGRI.myPhonyComponent(PP(i))-ComponentsLimit]; // reversed for coco4compatibility
   }
//std::clog<<"DeEmbedPoly:v "<<v<<endl;
//std::clog<<"DeEmbedPoly:end"<<endl;
   return v;
  }//DeEmbedPoly


  VectorList DeEmbedPolyList(const PolyList& PL,
                             const GRingInfo& theGRI)
  {
    return DeEmbedPolyList(PL,theGRI,0);
  }//DeEmbedPolyList


  VectorList DeEmbedPolyList(const PolyList& PL,
                             const GRingInfo& theGRI,
			     const long ComponentsLimit) // Poly whose LPP has last var exponent bigger than this number disappear on DeEmbedding

  {
    VectorList VL;
    if (PL.empty())
      return VL;
    PolyList::const_iterator it;
    for (it=PL.begin();it!=PL.end();++it)
      if (theGRI.myComponent(LPP(*it))<=theGRI.myComponent(ComponentsLimit))
        VL.push_back(DeEmbedPoly(*it,theGRI,ComponentsLimit));
    return VL;
  }//DeEmbedPolyList


   RingElem DeEmbedPolyToP(ConstRefRingElem p,
                           const GRingInfo& theGRI) // Poly whose LPP has last var degree bigger than this number disappear on DeEmbedding
   {
     return theGRI.myNewP2OldP()(p);
   }//DeEmbedPolyToP


  PolyList DeEmbedPolyListToPL(const PolyList& PL,
                               const GRingInfo& theGRI,
                               const long ComponentsLimit) // Poly whose LPP has phony last var exponent lesser than this number disappear on DeEmbedding
  {
    PolyList result;
    if (PL.empty())
      return result;
    PolyList::const_iterator it;
    for (it=PL.begin();it!=PL.end();++it)
      if (theGRI.myComponent(LPP(*it))<=theGRI.myComponent(ComponentsLimit)) // Copies. Copying disappears when I work with GPolys 			
        //result.push_back((theGRI.myNewP2OldP()(*it)));
        result.push_back(DeEmbedPolyToP(*it,theGRI));
    return result;
  }//DeEmbedPolyList
			
      				
   void monic(PolyList& PL)
   {
     if (PL.empty())
       return;
     const SparsePolyRing P(owner(PL));
     PPMonoidElem tmp(PPM(P));
     for (PolyList::iterator it=PL.begin();it!=PL.end();++it)
       (*it)/=monomial(P,LC(*it),tmp);
   }

    void PolyList2GPolyListClear(PolyList& thePL,
                                 GPolyList& theGPL,
                                 const GRingInfo& theGRI)
    {
      theGPL.clear();
      for (PolyList::iterator it=thePL.begin();it!=thePL.end();++it)
        theGPL.push_back(GPoly(*it, theGRI, clear));
      thePL.clear();
    }//PolyList2GPolyListClear

    void GPolyList2PolyListClear(GPolyList& theGPL,PolyList& thePL)
    {
      thePL.clear();
      if (theGPL.empty())
        return;
      const SparsePolyRing SPR=owner(theGPL);
      for (GPolyList::iterator it=theGPL.begin();it!=theGPL.end();++it)
      {
        thePL.push_back(one(SPR));
        swap(thePL.back(),it->myPolyValue);
      }
      theGPL.clear();
    }//GPolyList2PolyListClear

    void PolyList2GPolyList(const PolyList& thePL,
                             GPolyList& theGPL,
                             const GRingInfo& theGRI)
    {
      theGPL.clear();
      for (PolyList::const_iterator it=thePL.begin();it!=thePL.end();++it)
        theGPL.push_back(GPoly(*it, theGRI));
    }//PolyList2GPolyList

    void GPolyList2PolyList(const GPolyList& theGPL,PolyList& thePL)
    {
      thePL.clear();
      if (theGPL.empty())
        return;
      const SparsePolyRing SPR=owner(theGPL);
      for (GPolyList::const_iterator it=theGPL.begin();it!=theGPL.end();++it)
        thePL.push_back(it->myPolyValue);
    }//GPolyList2PolyListClear


     // Just transform the PolyList in a GPolyList
    GPolyList EmbedPolyList(PolyList& thePL,
                            const GRingInfo& theGRI)
    {
      GPolyList GPL;
      PolyList2GPolyListClear(thePL, GPL, theGRI);
      return GPL;
    }

    void DeEmbedPoly(ModuleElem& theOutputV,
                      const GPoly& the_gp)
    {
      DeEmbedPoly(theOutputV,the_gp,0l);
    }//DeEmbedPoly theOutputP

    // some copies are unavoidable when deembedding
    void DeEmbedPoly(ModuleElem& theOutputV,
                      const GPoly& the_gp,
                      const long ComponentsLimit) // the component in p that goes to the 0 component of the output vector v. Lesser components of p go to higher component of v
    {
      const FreeModule FM=owner(theOutputV);
      const GRingInfo& GRI=the_gp.myGRingInfo();
      theOutputV=zero(FM);
      const std::vector<ModuleElem>& e = gens(FM);
      const RingHom& phi=GRI.myNewP2OldP();
      for (SparsePolyIter i=BeginIter(the_gp.myPoly()); !IsEnded(i); ++i)
        theOutputV+=phi(monomial(GRI.myNewSPR(),coeff(i),PP(i)))*
                    e[GRI.myPhonyComponent(PP(i))-ComponentsLimit]; // reversed for cocoa4compatibility
    }//DeEmbedPoly  theOutputV


    void DeEmbedPolyList(VectorList& theOutputVL,
                          const GPolyList& theGPL)
    {
      DeEmbedPolyList(theOutputVL,theGPL,0);
    }//DeEmbedPolyList theOutputVL

    void DeEmbedPolyList(VectorList& theOutputVL,
                         const GPolyList& theGPL,
			 const long ComponentsLimit) // Poly whose LPP has last var degree bigger than this number disappear on DeEmbedding
    {
      theOutputVL.clear();
      if (theGPL.empty())
        return;
      const FreeModule FM=owner(theOutputVL);
      const GRingInfo& GRI(GetGRI(theGPL));
      ModuleElem tmp(FM);
      GPolyList::const_iterator it;
      for (it=theGPL.begin();it!=theGPL.end();++it)
      if (GRI.myComponent(LPPForDiv(*it))<=GRI.myComponent(ComponentsLimit))			
       {
         DeEmbedPoly(tmp,*it,ComponentsLimit);
         theOutputVL.push_back(tmp);
       }
    }//DeEmbedPolyList theOutputVL

    void DeEmbedPoly(RingElem& theOutputP,
                      const GPoly& the_gp) // Poly whose LPP has last var degree bigger than this number disappear on DeEmbedding
    {
      theOutputP=the_gp.myGRingInfo().myNewP2OldP()(the_gp.myPoly());
    }//DeEmbedPoly theOutputP

    void DeEmbedPolyList(PolyList& theOutputPL,
                         const GPolyList& theGPL,
			 const long ComponentsLimit) // Poly whose LPP has last var degree bigger than this number disappear on DeEmbedding
    {
       theOutputPL.clear();
       if (theGPL.empty())
         return;
       const GRingInfo& GRI(GetGRI(theGPL));
       RingElem tmp(GRI.myNewSPR());
       GPolyList::const_iterator it;
       for (it=theGPL.begin();it!=theGPL.end();++it)
         if (GRI.myComponent(LPPForDiv(*it))<=GRI.myComponent(ComponentsLimit)) // Copies. Copying disappears when I work with GPolys 			
         {
           DeEmbedPoly(tmp,*it);
           theOutputPL.push_back(tmp);
         }
    }//DeEmbedPolyList theOutputPL



    void SyzEmbedGPolyList(GPolyList& theGPL)
    {
      const GRingInfo& GRI(GetGRI(theGPL));
      long k=NumCompts(GRI.myFreeModule());
      for (GPolyList::iterator it=theGPL.begin() ; it!=theGPL.end() ; ++it,++k)
      { // Added by JAA 2012/10/11
        RingElem EkY = GRI.myE(k)*GRI.myY(wdeg(*it));
        (*it).myAppendClear(EkY);
      }
     }//SyzEmbedGPolyList


    void IntEmbedGPolyList(GPolyList& theGPL1, GPolyList& theGPL2)
    {
      CoCoA_ASSERT(GetGRI(theGPL1)==GetGRI(theGPL1));
      if (theGPL1.empty()||theGPL2.empty())
        return;
      const GRingInfo& GRI(GetGRI(theGPL1));
      long k = NumCompts(GRI.myFreeModule());
      
      for (GPolyList::iterator it=theGPL2.begin();it!=theGPL2.end();++it,++k)
      { // Added by JAA 2012/10/11
        RingElem EmbeddedPoly = it->myPoly()*GRI.myE(k);
        (*it).myAppendClear(EmbeddedPoly);
      }
      theGPL1.splice(theGPL1.begin(),theGPL2);
    }//IntEmbedGPolyList


    void ColonEmbedGPolyList(GPolyList& theGPL, GPoly& the_gp)
    {
      CoCoA_ASSERT(GetGRI(theGPL)==the_gp.myGRingInfo());
      const GRingInfo& GRI(the_gp.myGRingInfo());
      const long NC = NumCompts(GRI.myFreeModule());
      RingElem tmp =  GRI.myE(NC)*GRI.myY(wdeg(the_gp)); // JAA 2012-10-11
      the_gp.myAppendClear(tmp);                         // JAA 2012-10-11
      theGPL.push_back(GPoly(GRI));
      theGPL.back().AssignClear(the_gp);
    }//ColonEmbedGPolyList


    // The ordering is supposed to be Deg compatible
    RingElem homog(ConstRefRingElem the_p, const std::vector<RingElem>& theY)
    {
      const SparsePolyRing SPR=owner(the_p);
      if (GradingDim(SPR) != len(theY))
        CoCoA_ERROR("incompatible GradingDim", "homog(f, y1_yn");
      RingElem the_hp(SPR);
      RingElem tmp(SPR);
      degree MaxDeg(GradingDim(SPR));
      degree TmpDeg(GradingDim(SPR));
      // Compute the maximum degree
      for (SparsePolyIter it=BeginIter(the_p);!IsEnded(it);++it)
        MaxDeg=top(MaxDeg,wdeg(PP(it)));
      // Homogenizing
      for (SparsePolyIter it=BeginIter(the_p);!IsEnded(it);++it)
        {
          //std::clog<<"Homogenized:tmp"<<tmp<<endl;
          tmp = monomial(SPR,coeff(it),PP(it));
          TmpDeg=wdeg(PP(it));
          for (long i=0; i != len(theY); ++i)
            tmp*=power(theY[i],(MaxDeg[i]-TmpDeg[i]));
          SPR->myAddClear(raw(the_hp), raw(tmp));
        }
      return the_hp;
     }//homog (multihomog)


    void homogenized(ModuleElem& /*the_hv*/,
                     const ModuleElem& /*the_v*/,
                     const GRingInfo& /*theGRI*/)
    {
      CoCoA_ERROR(ERR::NYI, "homogenized(ModuleElem&, ModuleElem&, const GRingInfo");
    }//homogenized
/*
      const FreeModule NewFM=owner(the_hv);
      the_hv=zero(NewFM);
      const std::vector<ModuleElem>& e = gens(NewFM);
      degree maxdegree=wdeg(the_v);
      for (long i=0;i<NumCompts(owner(the_v));++i)
        for (SparsePolyIter it=BeginIter(the_v[i]); !IsEnded(it); ++it)
          the_hv=the_hv+
                   e[i]*(theGRI.myOldP2NewP())(monomial(theGRI.myNewSPR(),coeff(it),PP(it)))
                   *theGRI.myY(maxdegree-wdeg(PP(it)));
    }//Homogeneize
*/

    std::vector<long> PolyList2IndexList(const PolyList& thePL)
    {
      std::vector<long> result;
      std::vector<long> tmp;
      long i=0;

      for (PolyList::const_iterator it=thePL.begin();it!=thePL.end();++it)
      {
        SparsePolyIter it_p=BeginIter(*it);
        exponents(tmp,PP(it_p));
        i=0;
        while (i!=len(tmp)&&tmp[i]==0)
          ++i;
        CoCoA_ASSERT(i!=len(tmp));
        result.push_back(i);
      }
      return result;
    }


  PPMonoidElem IndexList2PPMonoidElem(const PPMonoid& thePPM,
                                      const std::vector<long>& theIL)
  {
    PPMonoidElem t(thePPM);
    for (std::vector<long>::const_iterator it=theIL.begin();it!=theIL.end();++it)
      t*=indet(thePPM,*it);
    return t;
  }//IndexList2PPMonoidElem


  std::vector<long> PPMonoidElem2IndexList(ConstRefPPMonoidElem thePP)
  {
    std::vector<long> IndexList;
    std::vector<long> tmp;
    exponents(tmp,thePP);
    for (long i=0; i < len(tmp);++i)
      if (tmp[i]!=0)
        IndexList.push_back(i);
    return IndexList;
  }//PPMonoidElem2IndexList


  SparsePolyRing MakeElimRingFromOld(const SparsePolyRing& theOldP,
                                     const std::vector<long>& IndetsToElim,
                                     const bool IsHomog)
  {
    matrix OldOrdOMat = OrdMat(theOldP);
    vector<long> rows = LongRange(0, GradingDim(theOldP)-1);
    vector<long> cols = LongRange(0, NumIndets(theOldP)-1);
    
    MatrixView Grading = submat(OldOrdOMat, rows, cols);

    matrix NewOrdMat = ElimMat(Grading, IndetsToElim);
    long NewGrDim = 0;
    if (GradingDim(theOldP) != 0 && IsHomog)
    {
      NewOrdMat = HomogElimMat(Grading, IndetsToElim);
      NewGrDim = GradingDim(theOldP);
    }

    const PPOrdering& NewOrd = NewMatrixOrdering(NumIndets(theOldP),
                                                 NewGrDim,
                                                 NewOrdMat);
    vector<symbol> IndetNames = NewSymbols(NumIndets(theOldP));
    SparsePolyRing NewP=NewSparsePolyRing(CoeffRing(theOldP),IndetNames,NewOrd);
    return NewP;
  }//MakeElimRingFromOld


    bool IsHomog(const PolyList& v)
    {
      // morally:  return find_if(v.begin(), v.end(), not1(IsHomog)) == v.end();
      const long n = len(v);
      for (long i=0; i < n; ++i)
	if (!IsHomog(v[i])) return false;
      return true;
//    We *DO NOT USE* STL algorithm because ptr_fun fails when fun has arg of reference type.
//       return find_if(v.begin(), v.end(),
//                      not1(ptr_fun(static_cast<bool(*)(ConstRefRingElem)>
//                                   (CoCoA::IsHomog))))
//         == v.end();
    }



    bool IsHomog(const VectorList& theVL)
    {
//       // morally:  return find_if(v.begin(), v.end(), not1(IsHomog)) == v.end();
//       return find_if(v.begin(), v.end(),
//                      not1(ptr_fun(static_cast<bool(*)(const ModuleElem&)>
//                                   (CoCoA::IsHomog))))
//         == v.end();
      if (theVL.empty())
        return true;
      for (VectorList::const_iterator it=theVL.begin();it!=theVL.end();++it)
        if (!IsHomog(*it))  return false;
      return true;
    }


/////////////////
/*
    std::vector<degree> DegStructure(ConstRefRingElem the_p)
    {
      std::vector<degree> L;
      for (SparsePolyIter it=BeginIter(the_p);!IsEnded(it); ++it)
        L.push_back(wdeg(PP(it)));
      return L;
    }//DegStructure

   std::vector<std::vector<degree> > DegStructure(const ModuleElem& the_v)
    {
      std::vector<std::vector<degree> > L;
      std::vector<degree> tmp;
      const FreeModule FM=owner(the_v);
      const std::vector<degree> S=shifts(FM);
      for (long i=0;i<NumCompts(owner(the_v));++i)
      {
        tmp.clear();
        for (SparsePolyIter it=BeginIter(the_v[i]);!IsEnded(it); ++it)
            tmp.push_back(wdeg(PP(it))+S[i]);
        L.push_back(tmp);
      }
      return L;								
    }//DegStructure

    std::vector<std::vector<degree> > DegStructure(const PolyList& thePL)
    {
      std::vector<std::vector<degree> > L;
      for (PolyList::const_iterator it=thePL.begin();it!=thePL.end();++it)
          L.push_back(DegStructure(*it));
      return L;
    }//DegStructure



    std::vector<std::vector<std::vector<degree> > > DegStructure(const VectorList& theVL)
    {
      std::vector<std::vector<std::vector<degree> > > L;
      for (VectorList::const_iterator it=theVL.begin();it!=theVL.end();++it)
          L.push_back(DegStructure(*it));
      return L;
    }//DegStructure

*/

    PolyList MakePolyList(ConstRefRingElem the_p)
    {
      PolyList L;
      L.push_back(the_p);
      return L;
    }//MakePolyList

    VectorList MakeVectorList(const ModuleElem& the_v)
    {
     VectorList L;
      L.push_back(the_v);
      return L;
    }//MakeVectorList


  PolyList WithoutDenominators(const PolyList& PL, SparsePolyRing Rx)
  {
    if (PL.empty()) return PL;
    CoCoA_ASSERT(HasUniqueOwner(PL));
    const SparsePolyRing Kx(owner(PL));
    CoCoA_ASSERT(IsFractionFieldOfGCDDomain(CoeffRing(Kx)));
    CoCoA_ASSERT(BaseRing(CoeffRing(Kx)) == CoeffRing(Rx));
    CoCoA_ASSERT(ordering(PPM(Kx)) == ordering(PPM(Rx)));
    PolyList result;
    std::vector<long> e(NumIndets(Rx));
    RingElem f(Rx);
    for (PolyList::const_iterator it=PL.begin() ; it!=PL.end() ; ++it )
    {
      RingElem d(CommonDenom(*it));
      f = 0;
      for (SparsePolyIter i=BeginIter(*it); !IsEnded(i); ++i)
      {
        exponents(e, PP(i));
        PushBack(f, (d/den(coeff(i))*num(coeff(i))), e); // same PPO
        //        f += monomial(Rx, (d/den(coeff(i))*num(coeff(i))), e);
      }
      result.push_back(f);
    }
    return result;
  }//WithoutDenominators


  PolyList WithDenominator1Hom(const PolyList& PL, SparsePolyRing Kx)
  {
    CoCoA_ASSERT(IsFractionField(CoeffRing(Kx)));
    if (PL.empty()) return PL;
    CoCoA_ASSERT(HasUniqueOwner(PL));
    const SparsePolyRing Rx(owner(PL));
    CoCoA_ASSERT(BaseRing(CoeffRing(Kx)) == CoeffRing(Rx));
    CoCoA_ASSERT(ordering(PPM(Kx)) == ordering(PPM(Rx)));
    PolyList result;
    // result.reserve(len(PL)); // if PL is a vector
    const RingHom EmbedIntoFrF = EmbeddingHom(CoeffRing(Kx));
    const RingHom EmbedCoeff   = CoeffEmbeddingHom(Kx);
    const RingHom phi = PolyRingHom(Rx, Kx, EmbedCoeff(EmbedIntoFrF), indets(Kx));
    //    transform(PL.begin(), PL.end(), back_inserter(result), phi);
    return apply(phi, PL);
  }//WithDenominator1Hom

// void ReadInt(std::istream& in, int& the_int,SkipTagType ST)
// {
//   if (ST == GetTag) SkipTag(in, "<int>");
//   in >> the_int;
//   SkipTag(in, "</int>");
// }//ReadInt

//////////////////

}// end namespace cocoa


// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/TmpGReductor.C,v 1.67 2014/07/31 14:45:18 abbott Exp $
// $Log: TmpGReductor.C,v $
// Revision 1.67  2014/07/31 14:45:18  abbott
// Summary: Merged io.H and UtilsTemplate.H into new header VectorOperations.H
// Author: JAA
//
// Revision 1.66  2014/07/31 13:10:45  bigatti
// -- GetMatrix(PPO) --> OrdMat(PPO)
// -- added OrdMat and GradingMat to PPOrdering, PPMonoid, SparsePolyRing
//
// Revision 1.65  2014/07/30 14:10:32  abbott
// Summary: Changed BaseRing into RingOf
// Author: JAA
//
// Revision 1.64  2014/07/22 10:47:33  bigatti
// -- using NewSymbols for elim
//
// Revision 1.63  2014/07/14 15:08:16  abbott
// Summary: Changed include of tmp.H into UtilsTemplate.H
// Author: JAA
//
// Revision 1.62  2014/07/09 14:27:53  abbott
// Summary: Removed AsFreeModule and AsFGModule
// Author: JAA
//
// Revision 1.61  2014/07/08 08:39:48  abbott
// Summary: Removed AsFractionField
// Author: JAA
//
// Revision 1.60  2014/07/07 13:06:32  abbott
// Summary: Removed AsSparsePolyRing
// Author: JAA
//
// Revision 1.59  2014/04/30 16:16:04  abbott
// Summary: Replaced X.size() by len(X)
// Author: JAA
//
// Revision 1.58  2014/04/17 13:38:57  bigatti
// -- MatrixViews --> MatrixView
//
// Revision 1.57  2014/04/11 15:44:27  abbott
// Summary: Renamed MatrixArith to MatrixOperations (in includes)
// Author: JAA
//
// Revision 1.56  2014/03/26 15:23:44  bigatti
// -- added MinGens for submodules
//
// Revision 1.55  2013/12/03 12:09:11  bigatti
// -- now incrementing myAgeValue only if added to GB
//
// Revision 1.54  2013/10/28 13:18:05  bigatti
// -- myCreateSpecialPairs --> myCreateInputPolyPairs
//
// Revision 1.53  2013/06/12 08:52:19  bigatti
// -- added mySetMinimalGen
// -- commented out unused code
//
// Revision 1.52  2013/06/03 09:11:50  bigatti
// renamed ModuleTermOrdering into ModuleOrdering
//
// Revision 1.51  2013/05/27 17:00:21  bigatti
// -- new names for ModuleOrdering
//
// Revision 1.50  2013/03/26 14:56:06  abbott
// Updated the conversion fns (in ptic removed procedure "convert");
// numerous consequential changes.
//
// Revision 1.49  2013/02/21 17:17:13  bigatti
// -- moved NewFreeModuleForSyz to FreeModule
//
// Revision 1.48  2013/02/14 17:34:35  bigatti
// -- cleaned up code for elimination matrices
//
// Revision 1.47  2012/10/24 12:23:27  abbott
// Removed use of std::ptr_fun.
//
// Revision 1.46  2012/10/16 09:54:27  abbott
// Replaced  RefRingElem  by  RingElem&
// Modified many calls to  myAppendClear  (because a temporary value will
// not match a non-const reference in a fn call).
//
// Revision 1.45  2012/10/03 12:26:08  bigatti
// -- changed homogenized --> homog: now returns a RingElem (instead of void)
//
// Revision 1.44  2012/05/28 16:47:23  bigatti
// -- added useful comment for debugging on anonymous symbols
//
// Revision 1.43  2012/05/28 09:18:20  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.42  2012/05/22 10:02:37  abbott
// Removed IsGCDDomain; substituted by IsTrueGCDDomain.
// Added IsFractionFieldOfGCDDomain.
//
// Revision 1.41  2012/05/10 14:44:53  abbott
// Employed anonymous symbols in MakeNewPRingFromModule (fixing issue #155).
//
// Revision 1.40  2012/03/30 17:31:04  bigatti
// -- making ordering matrices over QQ
//
// Revision 1.39  2012/02/10 10:29:07  bigatti
// -- changed RingZ.H, RingQ.H --> RingZZ.H, RingQQ.H
//
// Revision 1.38  2012/02/08 16:14:49  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.37  2012/01/26 16:49:51  bigatti
// -- using apply instead of transform
//
// Revision 1.36  2011/12/07 16:00:07  bigatti
// -- removed double "#include <functional>"
//
// Revision 1.35  2011/08/14 15:52:16  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.34  2011/05/27 12:13:14  bigatti
// -- added CoCoA_ASSERT for empty list of generators
//
// Revision 1.33  2011/05/24 08:16:00  bigatti
// -- only spaces
//
// Revision 1.32  2011/03/22 11:26:53  bigatti
// -- unsigned long --> long
// -- myAge --> myAgeValue
//
// Revision 1.31  2011/03/11 17:39:10  bigatti
// -- changed  unsigned int --> long
//
// Revision 1.30  2011/03/11 15:39:05  bigatti
// -- changes size_t --> long
// -- changes size --> len
//
// Revision 1.29  2011/03/11 14:50:39  abbott
// Changed some size_t into long.
//
// Revision 1.28  2011/03/10 16:39:33  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.27  2011/03/08 18:06:33  bigatti
// -- removed useless vector<long> copy in MakeElimRingFromOld
//
// Revision 1.26  2010/12/26 13:04:37  abbott
// Changed "GlobalXXXput" into corresponding std C++ stream
// (even in commented out code).
//
// Revision 1.25  2010/08/06 17:11:02  bigatti
// -- restored owner(const VectorList&)  // should be moved into TmpGTypes
//
// Revision 1.24  2010/08/06 10:28:15  bigatti
// -- commented out myDoAFFGBasis and owner(VL): unused
//
// Revision 1.23  2010/08/05 16:05:36  bigatti
// -- simpler code after new class criteria
//
// Revision 1.22  2010/07/16 09:29:58  bigatti
// -- minor cleaning and coding conventions
//
// Revision 1.21  2010/05/14 09:53:09  bigatti
// -- removed empty ctor for SugarDegree
// -- added marker for SugarDegree(uninitialized)
// -- SugarDegree for GBasis input is initialized by myPrepareGBasis
//
// Revision 1.20  2010/03/23 14:43:07  bigatti
// -- class GRingInfo estracted from TmpGPoly
//
// Revision 1.19  2009/11/26 16:59:19  bigatti
// -- minor cleaning of WithoutDenominators and WithDenominator1Hom
// -- sorted includes
//
// Revision 1.18  2009/10/26 16:29:57  bigatti
// -- minor change in creting a zero GPoly
// -- added printing of sugar when starting a new degree
//
// Revision 1.17  2009/09/22 13:35:55  bigatti
// -- following coding conventions in function names Matrix --> Mat
// -- forced all matrices to be over RingZ
//
// Revision 1.16  2009/05/14 12:30:33  abbott
// Added missing colon.
//
// Revision 1.15  2009/04/27 14:11:22  bigatti
// -- added BuchbergerOpTypeMarker flag
// -- added cocoa_assert
//
// Revision 1.14  2009/03/20 15:33:26  bigatti
// -- just comments and spaces
//
// Revision 1.13  2009/01/30 13:41:50  bigatti
// -- enum instead of bool arguments
//
// Revision 1.12  2008/09/22 15:58:58  bigatti
// -- commented out unused "double T=0.0"
//
// Revision 1.11  2008/09/19 14:08:16  bigatti
// -- modified GRStats (M.Caboara)
//
// Revision 1.10  2008/09/16 15:03:42  bigatti
// -- added LPPForDiv
// -- changed LPP into LPPForOrd
//
// Revision 1.9  2008/07/09 16:10:11  abbott
// Added missing & (for C++ const reference).
//
// Revision 1.8  2008/04/21 11:23:11  abbott
// Separated functions dealing with matrices and PPOrderings into a new file.
// Added matrix norms, and completed adjoint.
//
// Revision 1.7  2008/04/18 15:35:57  abbott
// (long overdue) Major revision to matrices
//
// Revision 1.6  2008/03/12 16:35:18  bigatti
// -- changed: IsHomogeneous --> IsHomog
// -- changed: ERR:ZeroPoly --> ERR::ZeroRingElem
//
// Revision 1.5  2007/11/09 10:45:52  bigatti
// -- [caboara] preparation for self-saturating algorithm
//
// Revision 1.4  2007/10/30 17:14:07  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.3  2007/03/23 18:38:42  abbott
// Separated the "convert" functions (and CheckedCast) into their own files.
// Many consequential changes.  Also corrected conversion to doubles.
//
// Revision 1.2  2007/03/12 16:31:10  bigatti
// -- fixed iterators issues (M.Abshoff)
//
// Revision 1.1  2007/03/09 18:56:56  bigatti
// -- added Tmp prefix to Groebner related files
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.30  2007/03/08 18:22:29  cocoa
// Just whitespace cleaning.
//
// Revision 1.29  2007/03/08 16:55:06  cocoa
// Changed name of "range" function to "SymbolRange".
//
// Revision 1.28  2007/03/08 11:07:12  cocoa
// Made pseudo ctors for polynomial rings more uniform.  This allowed me to
// remove an include of CoCoA/symbol.H  from the RingDistrM*.H files, but then
// I had to put the include in several .C files.
//
// Revision 1.27  2007/03/07 17:04:31  cocoa
// -- several changes by M.Caboara: more operations on ideals,
//    exception cleaner, coding conventions, WSugar, dynamic
//
// Revision 1.26  2007/03/05 21:06:07  cocoa
// New names for homomorphism pseudo-ctors: removed the "New" prefix.
//
// Revision 1.25  2007/03/02 10:47:53  cocoa
// First stage of RingZ modifications -- tests do not compile currently, Anna will fix this.
//
// Revision 1.24  2007/02/12 16:44:05  bigatti
// -- restore "GlobalLogput()" in logging info -- to be discussed with Max
//
// Revision 1.23  2007/02/10 18:44:03  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.22  2006/12/21 13:48:33  cocoa
// Made all increment/decrement calls prefix (except where the must be postfix).
//
// Revision 1.21  2006/12/04 13:55:54  cocoa
// -- added: sugar for GradingDim > 0  (called wsugar)
//
// Revision 1.20  2006/11/27 14:24:49  cocoa
// -- reorganised #include files
//
// Revision 1.19  2006/11/24 17:39:14  cocoa
// -- reorganized includes of header files
//
// Revision 1.18  2006/11/22 15:33:39  cocoa
// -- changed handling of timings and number of reductions (M.Caboara)
//
// Revision 1.17  2006/11/20 14:57:17  cocoa
// -- added: (standard) sugar for modules
// -- fixed: non-homogeneous sysygies
// -- minor fixes     [M.Caboara]
//
// Revision 1.16  2006/11/14 17:20:23  cocoa
// -- uncommented myDoAFFGBasis for RingWeyl (should be updated soon)
//
// Revision 1.15  2006/11/09 17:36:06  cocoa
// -- bug fix for WDegPos ordering (by M.Caboara)
//
// Revision 1.14  2006/10/06 17:32:26  cocoa
// -- wip: evolution of Groebner Framework (Max)
//
// Revision 1.13  2006/10/06 16:32:06  cocoa
// -- changed: GPoly::SPoly --> GPoly::myAssignSPoly
// -- changed: Len(const GPoly&) --> len(const GPoly&)
// -- added:   poly(const GPoly&)
// -- added:   GPoly::myUpdateLenLPPDegComp()
// -- in reduce.C added functions for computing sugar during reduction
//
// Revision 1.12  2006/10/06 14:04:15  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.11  2006/10/06 10:48:34  cocoa
// Removed #include references to GradedFreeModule.H
//
// Revision 1.10  2006/10/06 10:16:26  cocoa
// Replaced two dynamic_casts by calls to the new function IsRingFp.
//
// Revision 1.9  2006/09/27 14:33:50  cocoa
// -- fixed check for RingFp in NewSparsePolyRing
//
// Revision 1.8  2006/08/17 09:32:26  cocoa
// -- added: flags for homogeneous input
//
// Revision 1.7  2006/08/07 21:23:25  cocoa
// Removed almost all publicly visible references to SmallExponent_t;
// changed to long in all PPMonoid functions and SparsePolyRing functions.
// DivMask remains to sorted out.
//
// Revision 1.6  2006/07/20 14:43:56  cocoa
// -- just spaces and a CoCoA_ASSERT
//
// Revision 1.5  2006/07/18 11:01:00  cocoa
// -- changed symbol for the names of the "module" indeterminates
// -- added variables for calling ring homomorphisms
// -- moved IsHomogeneous to SparsePolyRing/FreeModule
//
// Revision 1.4  2006/06/20 17:17:45  cocoa
// -- reinserted TotalTime printing: this is to compare timings for
//    copying the result
//
// Revision 1.3  2006/06/13 16:40:55  cocoa
// -- simpler code in CoCoAServer to print length of the result and time
//    (easier for benchmarcks)
//
// Revision 1.2  2006/06/09 13:24:53  cocoa
// -- just a "\t" before printing "TotalTime"
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.24  2006/05/16 09:01:26  cocoa
// -- added "const" to arg of homogenized
// -- changed ordering matrix for modules (compatibility with GradedFreeModule)
//
// Revision 1.23  2006/05/12 17:03:16  cocoa
// -- swapped arguments in homogenized
//
// Revision 1.22  2006/05/12 10:46:53  cocoa
// -- commented out printouts
//
// Revision 1.21  2006/05/11 16:00:22  cocoa
// -- fixed spelling of "homogenize"
//
// Revision 1.20  2006/05/04 14:25:16  cocoa
// -- major cleaning of FreeModule: created GradedFreeModule and moved
//    some code around
//
// Revision 1.19  2006/05/02 14:37:12  cocoa
// -- only changes to logging info (by M.Abshoff)
//
// Revision 1.18  2006/04/26 16:59:27  cocoa
// -- restored printing of "time" (for benchmarks)
//
// Revision 1.17  2006/04/26 09:50:40  cocoa
// -- added GReductor::ourDefaultStatLevel variable to allow CoCoAServer
//    to set statistics level
//
// Revision 1.16  2006/04/21 16:45:17  cocoa
// -- new functions by Max
//
// Revision 1.15  2006/04/11 16:59:52  cocoa
// -- commented out "sort" in myUpdateBasisAndPairs: probably unnecessary
//
// Revision 1.14  2006/04/11 16:44:19  cocoa
// -- added functions for Elim
// -- minor cleanup
//
// Revision 1.13  2006/04/10 17:02:47  cocoa
// -- changed: BCriterion_OK now uses PPWithMask instead of ConstRefPPMonoidElem
//
// Revision 1.12  2006/03/17 18:17:16  cocoa
// -- changed: use of ReductionCog for reduction (major cleanup)
//
// Revision 1.11  2006/03/02 13:45:57  cocoa
// -- changes by Max
//
// Revision 1.10  2006/02/13 14:46:45  cocoa
// -- changes by Max
//
// Revision 1.9  2006/02/13 13:53:04  cocoa
// -- changes by Max (GRingInfo)
//
// Revision 1.8
// date: 2006/01/20 15:43:30;  author: cocoa;  state: Exp;  lines: +6 -3
// -- fixed: use of RefPPMonoidElem and ConstRefPPMonoidElem

// Revision 1.7
// date: 2006/01/18 15:58:20;  author: cocoa;  state: Exp;  lines: +133 -131
// -- new changes by Max

// Revision 1.6  2006/01/17 15:44:56  cocoa
// -- chamges by Max for operations with modules
//
// Revision 1.5  2006/01/17 10:23:08  cocoa
// Updated DivMask; many consequential changes.
// A few other minor fixes.
//
// Revision 1.4  2005/12/31 12:22:17  cocoa
// Several minor tweaks to silence the Microsoft compiler:
//  - added some missing #includes and using directives
//  - moved some function defns into the right namespace
//  - etc.
//
// Revision 1.3  2005/11/17 13:21:54  cocoa
// -- used "zero(myPolyRing)"
//
// Revision 1.2  2005/11/15 13:05:42  cocoa
// Forgot to check these before checking in the new timing functions.
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
// Revision 1.3  2005/04/19 14:06:04  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.2  2005/02/11 14:15:20  cocoa
// New style ring elements and references to ring elements;
// I hope I have finally got it right!
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.11  2004/11/19 15:44:27  cocoa
// Changed names of "casting" functions which convert a ring into
// one with a more special structure (e.g. FractionField).  These
// functions now have names starting with "As".  There were several
// consequential changes.
//
// Revision 1.10  2004/11/02 15:58:33  cocoa
// -- added some (commented) lines for new DivMask::SingleBitWrapImpl
//
// Revision 1.9  2004/10/29 16:02:52  cocoa
// -- added field myDivMaskImplPtr for creating GPolys with the same
//    DivMask::base
// -- calls to LCMwMask (instead of LCM) for fast divisibility tests
// -- function IsDivisible had wrong semantics --> swapped arguments everywhere
//
// Revision 1.8  2004/06/29 17:10:22  cocoa
// Partially tidied use of "protected" and "private" in various
// base classes.  Checking in at the end of the day -- it works,
// and I wouldn't want it to be lost next time point's disk
// misbehaves.
//
// Revision 1.7  2004/06/16 16:13:41  cocoa
// Improved I/O facilities with knock-on changes
//
// Revision 1.6  2004/05/27 16:52:12  cocoa
// -- removed ";" at the end of function bodies (g++ gives error on them)
//
// Revision 1.5  2004/03/20 17:46:10  cocoa
// Check in prior to departure to RWCA
//
// Revision 1.4  2004/03/04 11:37:18  cocoa
// -- updated code for Borel reductors:
//    ReductorData fields: myGPoly->myGPolyPtr;  new: myCount, IamBorelUpdated
//    myBorelReductors is now in Reductors (was in GReductor)
//    Reductors: field: IhaveBorelReductors;  type: enum UseBorelMarker
//
// Revision 1.3  2003/10/09 12:48:17  cocoa
// New coding convention for rings.
//
// Revision 1.2  2003/10/01 10:35:31  cocoa
// - applied "my" coding convention to PPMonoid and PPOrdering
//
// Revision 1.1.1.1  2003/09/24 12:55:43  cocoa
// Imported files
//
// Revision 1.15  2003/09/22 17:30:10  bigatti
// - new field myOrdPoly to order the GPairs
// - ==> not ordered list of pairs which are completed and sorted later
//
// Revision 1.14  2003/06/23 17:10:00  abbott
// Minor cleaning prior to public release.
// Improved the include directives,
//
// Revision 1.13  2003/05/29 16:56:56  bigatti
// - change: myRingSpecialIndex is now an int
// - added:  flag for criterion (to be disables for non-comm Groebner)
// - fix:    all zero polynomials are now removed from input list
// - change: each GReductor may now set its own statistics level
//
// Revision 1.12  2003/05/28 14:29:43  bigatti
// - new code for modules
//
// Revision 1.11  2003/05/14 17:09:28  bigatti
// - new functions for "BorelReductorsPolys" and saturating algorithm
// - optimization of PPMonoidElem function calls
//
// Revision 1.10  2002/11/15 16:46:40  bigatti
// - clear code after new PP.H (with DivMask)
//
// Revision 1.9  2002/09/19 17:23:32  bigatti
// - Cleaner code based on PolyRing
//
// Revision 1.8  2002/05/13 11:48:36  bigatti
// - new data structure for "Reductors"
//
// Revision 1.7  2002/04/29 17:06:29  bigatti
// - removed "true" from "BuildNewPairs" [if (true||myStat.print_poly_deleted)..]
//
// Revision 1.6  2002/04/15 17:18:59  bigatti
// - Max's new code
//
// Revision 1.5  2002/04/09 14:20:27  bigatti
// - SPoly now takes a GPair as argument
// - CKR criterion
//
// Revision 1.4  2002/01/31 17:23:39  bigatti
// - interreduce
//
// Revision 1.3  2002/01/11 10:34:34  bigatti
// - removed all "--myGB.end()"
//
// Revision 1.2  2001/12/12 18:16:39  bigatti
// - new structure of reduction
//
// Revision 1.1  2001/12/05 13:25:45  bigatti
// Initial revision
//

// Ordinamento delle coppie - cambiare ?
// inline e controllo efficienza
// Sugar
// liste-vettori   poly con rc
// risistemazioni delle stampa
// GPairs con i cofattori
// poly const * f - per theGPoly e per Pairs, se non c'e' il rc
// current paur? SPoly? or the actual structure using the last of .. .. vedere se c'e' rc
// ordinamento dei riduttori
// le statistiche deg by deg sono appiccicate con lo sputo, in quanto non
//   matchano con la struttura.
