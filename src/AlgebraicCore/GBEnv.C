//   Copyright (c)  2005-2010 Anna Bigatti, Massimo Caboara

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


#include "CoCoA/GBEnv.H"
#include "CoCoA/DenseMatrix.H" // for DetermineIfMyGradingIsPosPlus
#include "CoCoA/FreeModule.H"
#include "CoCoA/IntOperations.H"
#include "CoCoA/PPMonoidEv.H"
#include "CoCoA/RingQQ.H" // for DetermineIfMyGradingIsPosPlus
#include "CoCoA/assert.H"
#include "CoCoA/VectorOperations.H" // just for debugging and statistics
#include "CoCoA/matrix.H" // for DetermineIfMyGradingIsPosPlus
#include "CoCoA/symbol.H"

using std::vector;
#include <limits>
using std::numeric_limits;
#include <algorithm>
using std::min;
using std::max;
#include <iostream>
using std::ostream;
using std::endl;
//using std::swap;
#include <iterator>



namespace CoCoA
{  

// Two utilities for GRingInfo ctors

  RingHom CreateNewP2OldPRingHom(const SparsePolyRing& theNewSPR,
                                 const SparsePolyRing& theOldSPR)
  {
    if (theNewSPR==theOldSPR)
      return  IdentityHom(theNewSPR);
    std::vector<RingElem> images;
    for (long i=0; i!=NumIndets(theOldSPR); ++i)
      images.push_back(indet(theOldSPR, i)); // x[i] |-> x[i]
    for (long i=0; i!=GradingDim(theNewSPR); ++i)
      images.push_back(one(theOldSPR)); // y[i] |-> 1
    images.push_back(one(theOldSPR)); // e |-> 1
    return  PolyRingHom(theNewSPR,theOldSPR,CoeffEmbeddingHom(theOldSPR),images);
  }//CreateNewP2OldPRingHom

  RingHom CreateOldP2NewPRingHom(const SparsePolyRing& theNewSPR,
                                 const SparsePolyRing& theOldSPR)
  {
     if (theNewSPR==theOldSPR)
       return  IdentityHom(theNewSPR);
     std::vector<RingElem> images;
     for (long i=0; i < NumIndets(theOldSPR); ++i)
       images.push_back(indet(theNewSPR, i));       // x[i] |-> x[i]
     return PolyRingHom(theOldSPR,theNewSPR,CoeffEmbeddingHom(theNewSPR),images);
  }//CreateOldP2NewPRingHom





/*
---------class GRingInfo-------------------------------
*/

ComputationInputAndGradingType  DetermineComputationType(long GrDim,
                                                         const bool IsHomog,
							 const bool IsSatAlg)
{
// This check is temporarily disable to experimento with the bigrading idea.
//clog<<"DetermineComputationType:IsSatAlg="<<IsSatAlg<<std::endl;
//clog<<"DetermineComputationType:GrDim="<<GrDim<<std::endl;
  if (IsSatAlg)
  { 
    //clog<<"DetermineComputationType:=SaturatingAlgNoDRL"<<std::endl;
    if (!IsHomog) CoCoA_ERROR(ERR::SERIOUS, "DetermineComputationType");
    if (GrDim==0) return SaturatingAlgNoDRL;
    return SaturatingAlg;
  }
  if (GrDim==0) return AllAffine;
  if (IsHomog) return AllGraded;
  return AffineInputWithGradedRing;
}//DetermineComputationType

  // ----------------------------------------------------------------------
  // GRingInfo ctors

  void GRingInfo::myCtorAux(const SparsePolyRing& theNewSPR,
                            const bool HomogeneousInput,
                            const bool SaturatingAlgorithm)
  {
    myInputAndGradingValue=DetermineComputationType(GradingDim(theNewSPR),
                                                    HomogeneousInput,
                                                    SaturatingAlgorithm);
    myGradingPosPlusValue=DetermineIfMyGradingIsPosPlus(theNewSPR);
    mySetCoeffRingType(CoeffEncoding::Field);
  }
  

  GRingInfo::GRingInfo(const SparsePolyRing& theNewSPR,
                       const bool HomogeneousInput,
	               const bool SaturatingAlgorithm,
		       const DivMaskRule& DivMaskR):
    myNewSPRValue(theNewSPR),
    myOldSPRValue(theNewSPR),
    myPPMValue(NewPPMonoidEv(SymbolRange("x",0,NumIndets(theNewSPR)-1), ordering(PPM(theNewSPR)))),
    myFreeModuleValue(NewFreeModule(myNewSPRValue,1)),
    myOutputFreeModuleValue(NewFreeModule(myNewSPRValue,1)),
    myNewP2OldPValue(IdentityHom(myNewSPRValue)),
    myOldP2NewPValue(IdentityHom(myNewSPRValue)),
    myDivMaskRuleValue(DivMaskR),
    AmIModuleValue(false)
  {
    myCtorAux(theNewSPR, HomogeneousInput, SaturatingAlgorithm);
  }// ctor GRingInfo

  GRingInfo::GRingInfo(const SparsePolyRing& theNewSPR,
                       const SparsePolyRing& theOldSPR,
                       const FreeModule& theFM,
                       const FreeModule& theOutputFM,
                       const bool HomogeneousInput,
	               const bool SaturatingAlgorithm,
                       const DivMaskRule& DivMaskR):
    myNewSPRValue(theNewSPR),
    myOldSPRValue(theOldSPR),
    myPPMValue(NewPPMonoidEv(SymbolRange("x",0,NumIndets(theNewSPR)-1), ordering(PPM(theNewSPR)))),
    myFreeModuleValue(theFM),
    myOutputFreeModuleValue(theOutputFM),
    myNewP2OldPValue(CreateNewP2OldPRingHom(myNewSPRValue,myOldSPRValue)),
    myOldP2NewPValue(CreateOldP2NewPRingHom(myNewSPRValue,myOldSPRValue)),
    myDivMaskRuleValue(DivMaskR),
    AmIModuleValue(theNewSPR!=theOldSPR)
  {
    //    if (!IsField(CoeffRing(theOldSPR)))
    //      CoCoA_ERROR("coefficients are not in a field", "ComputeGBasis");
    std::vector<RingElem> Y; // The grading vars
    const std::vector<RingElem>& x = indets(myNewSPRValue);
    // Fill Y
    for (long i=0; i < GradingDim(myNewSPRValue); ++i)
       Y.push_back(x[i+NumIndets(theOldSPR)]);

    const std::vector<degree> S=shifts(myFreeModuleValue);
   RingElem tmp(myNewSPRValue);
   for (long i=0; i < NumCompts(myFreeModuleValue); ++i)
   {
      tmp=power(myE(),this->myComponent(i));
      for (long j=0; j < GradingDim(myNewSPRValue); ++j)
        tmp*=power(Y[j],S[i][j]);
      myEYValue.push_back(tmp);
    }
    myCtorAux(theNewSPR, HomogeneousInput, SaturatingAlgorithm);
  }// ctor GRingInfo

  GRingInfo::GRingInfo(const SparsePolyRing& theNewSPR,
                       const SparsePolyRing& theOldSPR,
                       const FreeModule& theOutputFM,
                       const bool HomogeneousInput,
	               const bool SaturatingAlgorithm,
                       const DivMaskRule& DivMaskR):
    myNewSPRValue(theNewSPR),
    myOldSPRValue(theOldSPR),
    myPPMValue(NewPPMonoidEv(SymbolRange("x",0,NumIndets(theNewSPR)-1), ordering(PPM(theNewSPR)))),
    myFreeModuleValue(NewFreeModule(theNewSPR,1)),
    myOutputFreeModuleValue(theOutputFM),
    myNewP2OldPValue(CreateNewP2OldPRingHom(myNewSPRValue,myOldSPRValue)),
    myOldP2NewPValue(CreateOldP2NewPRingHom(myNewSPRValue,myOldSPRValue)),
    myDivMaskRuleValue(DivMaskR),
    AmIModuleValue(theNewSPR!=theOldSPR)
  {
    //    if (!IsField(CoeffRing(theOldSPR)))
    //      CoCoA_ERROR("coefficients are not in a field", "ComputeGBasis");
    myCtorAux(theNewSPR, HomogeneousInput, SaturatingAlgorithm);
  }// ctor GRingInfo

  GRingInfo::GRingInfo(const SparsePolyRing& theNewSPR,
                       const SparsePolyRing& theOldSPR,
                       const bool HomogeneousInput,
	               const bool SaturatingAlgorithm,
                       const DivMaskRule& DivMaskR):
    myNewSPRValue(theNewSPR),
    myOldSPRValue(theOldSPR),
    myPPMValue(NewPPMonoidEv(SymbolRange("x",0,NumIndets(theNewSPR)-1), ordering(PPM(theNewSPR)))),
    myFreeModuleValue(NewFreeModule(theNewSPR,1)),
    myOutputFreeModuleValue(NewFreeModule(theNewSPR,1)),
    myNewP2OldPValue(CreateNewP2OldPRingHom(myNewSPRValue,myOldSPRValue)),
    myOldP2NewPValue(CreateOldP2NewPRingHom(myNewSPRValue,myOldSPRValue)),
    myDivMaskRuleValue(DivMaskR),
    AmIModuleValue(theNewSPR!=theOldSPR)
  {
    //    if (!IsField(CoeffRing(theOldSPR)))
    //      CoCoA_ERROR("coefficients are not in a field", "ComputeGBasis");
    myCtorAux(theNewSPR, HomogeneousInput, SaturatingAlgorithm);
  }// ctor GRingInfo

  // GRingInfo ctors
  // ----------------------------------------------------------------------

  void GRingInfo::mySetCoeffRingType(CoeffEncoding::type CT)
  { myCoeffRingTypeValue = CT; }


  bool GRingInfo::operator==(const GRingInfo& theGRI)const
  {
    if (myNewSPRValue==theGRI.myNewSPRValue
         &&
	myOldSPRValue==theGRI.myOldSPRValue
         &&
	myPPMValue==theGRI.myPPMValue
         &&
        myOutputFreeModuleValue==theGRI.myOutputFreeModuleValue
         &&
	myEYValue==theGRI.myEYValue
         //&& // I want to do this, the == operator is not there
         //myDivMaskRuleValue==theGRI.myDivMaskRuleValue
         )
      return true;
    else
      return false;
  }//operator==



long GRingInfo::myComponent(ConstRefPPMonoidElem T)const
{
  if (!AmIModule()) return 0;// True Ring
  return exponent(T,ModuleVarIndex(myNewSPRValue));
}

long GRingInfo::myPhonyComponent(ConstRefPPMonoidElem T)const
{
  if (!AmIModule()) return 0;// True Ring
  return myComponent(exponent(T,ModuleVarIndex(myNewSPRValue)));
}

RingElem GRingInfo::myY(const degree& the_d)const
{
   RingElem result(one(myNewSPR()));
   for (long j=0; j < GradingDim(myNewSPR()); ++j)
      result*=power(myY(j),the_d[j]);
   return result;
}//myY


  SugarDegree GRingInfo::myNewSugar(ConstRefRingElem f) const
  {
    switch (myInputAndGrading())
    {
    case AllGraded: // ANNA: (w)graded + homogeneous
      return NewWSugarConst(f);
    case SaturatingAlg:
      return NewWSugarSat(f);
    case AffineInputWithGradedRing: // ANNA: (w)graded + non-homogeneous
      {
        if (/*module && */ IsMyGradingPosPlus())
        { // ANNA: should be implemented with proper weights
          int idx = ModuleVarIndex(myNewSPR());
          return NewStdSugarNoIdx(f, idx);
        }
        return NewWDeg1CompTmp(f);
      }
    case AllAffine: // ANNA: GradingDim = 0 --> StandardSugarAlgorithm
      //      if (/*module && */ IsMyGradingPosPlus())
      if (AmIModule())
      {
        int idx = ModuleVarIndex(myNewSPR());
        return NewStdSugarNoIdx(f, idx);
      }
      return NewStdSugar(f);
    case SaturatingAlgNoDRL: // GradingDim = 0
      if (/*module && */ IsMyGradingPosPlus())
      {
        int idx = ModuleVarIndex(myNewSPR());
        return NewStdSugarNoIdxSat(f, idx);
      }
      return NewStdSugarSat(f);
    default: CoCoA_ERROR(ERR::SERIOUS, "GRingInfo::mySugar");
    }//switch
    return NewStdSugar(f); // just to keep the compiler quiet
  }


ostream& operator<<(ostream& out, const GRingInfo& theGRI)
{
  out<<"the ring is "<<theGRI.myNewSPR()<<endl
     <<" the old ring is "<<theGRI.myOldSPR()<<endl
     //<<" Input Free Module "<<theGRI.myFreeModule()<<endl
     //<<" Output Free Module "<<theGRI.myOutputFreeModule()<<endl
     <<" AmIModule "<<theGRI.AmIModule()<<endl
     <<" myInputAndGrading = "<<theGRI.myInputAndGrading()<<endl
     <<" myGradingPosPlusValue = "<<theGRI.IsMyGradingPosPlus()<<endl
     <<" embedding grading "
     <<" EY=\n";
     for (std::vector<RingElem>::const_iterator it=theGRI.myEYValue.begin();
                                                it!=theGRI.myEYValue.end();++it)
       {out<<*it<<endl;}						
     out<<endl;;
  return out;
}


long ModuleVarIndex(const SparsePolyRing& P)
{
  long tmp = NumIndets(P);
  if (tmp!=0)
    return tmp-1;
  else
    return tmp;
}//ModuleVarIndex


bool AreCompatible(const GRingInfo& theGRI1,const GRingInfo& theGRI2)
{
  if (theGRI1.myNewSPRValue==theGRI2.myNewSPRValue
        &&
      theGRI1.myOldSPRValue==theGRI2.myOldSPRValue
        &&
      theGRI1.myPPMValue==theGRI2.myPPMValue
      )
       //&& // I want to do this, the == operator is not there
         //theGRI1.myDivMaskRuleValue==theGRI2.myDivMaskRuleValue
    return true;
  else
    return false;
}//AreCompatible


// A member field?
std::vector<RingElem> GRingInfo::myY()const
{
  vector<RingElem> Y;
  for (long i=0; i < GradingDim(myNewSPRValue); ++i)
    Y.push_back(indet(myNewSPRValue,i+NumIndets(myOldSPRValue)));
  return Y;
}//myY()

 // Grdim>=2, order matrix first row is 0,..,0,1
 bool GRingInfo::DetermineIfMyGradingIsPosPlus(const SparsePolyRing& theSPR)
 {
   // This checks if indeed the order is a PosPlus.
   // Another option is to SET this field at the right time.
   // Slightly more efficient, but more risky.
   return false;
   if (GradingDim(theSPR)<1)
     return false;
   matrix OrdM(NewDenseMat(RingQQ(),NumIndets(theSPR),NumIndets(theSPR)));
   PPM(theSPR)->myOrdering()->myOrdMatCopy(OrdM);
   for (long i=0; i<NumIndets(theSPR)-1; ++i)
     if (OrdM(0,i)!=0)
     {
       return false;
     }
   if (OrdM(0,NumIndets(theSPR)-1)!=1)
   {
     return false;
   }
   return true;
 }

}// end namespace cocoa

// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/GBEnv.C,v 1.10 2014/07/31 14:45:17 abbott Exp $
// $Log: GBEnv.C,v $
// Revision 1.10  2014/07/31 14:45:17  abbott
// Summary: Merged io.H and UtilsTemplate.H into new header VectorOperations.H
// Author: JAA
//
// Revision 1.9  2014/07/31 13:10:46  bigatti
// -- GetMatrix(PPO) --> OrdMat(PPO)
// -- added OrdMat and GradingMat to PPOrdering, PPMonoid, SparsePolyRing
//
// Revision 1.8  2012/05/28 09:18:21  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.7  2012/04/02 15:28:23  bigatti
// -- ZZ --> QQ
//
// Revision 1.6  2012/02/10 10:26:40  bigatti
// -- changed RingZ.H, RingQ.H --> RingZZ.H, RingQQ.H
//
// Revision 1.5  2012/02/08 17:10:10  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.4  2011/03/11 16:50:04  bigatti
// -- changes  unsigned int  --> long
//
// Revision 1.3  2011/03/10 16:39:34  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.2  2010/12/26 13:04:37  abbott
// Changed "GlobalXXXput" into corresponding std C++ stream
// (even in commented out code).
//
// Revision 1.1  2010/03/23 14:40:55  bigatti
// -- first import
//
