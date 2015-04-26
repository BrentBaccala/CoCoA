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

#include "CoCoA/FractionField.H"
#include "CoCoA/FreeModule.H"
#include "CoCoA/IntOperations.H"
#include "CoCoA/RingDistrMPolyInlFpPP.H"
#include "CoCoA/RingDistrMPolyInlPP.H"
#include "CoCoA/TmpGOperations.H"
// #include "CoCoA/TmpGReductor.H"  // included in TmpGOperations.H 
#include "CoCoA/symbol.H"
#include "CoCoA/VectorOperations.H"


#include <algorithm>
using std::stable_sort;
#include <list>
using std::list;
#include <functional>
using std::less;
using std::binary_function;
using std::mem_fun_ref; // for calling GPair::complete on GPairList
#include <utility>
using std::make_pair;
#include <iostream>
using std::ostream;
using std::endl;
using std::flush;

namespace CoCoA
{


  namespace // anonymous
  {
    bool IsHomogGrD0(const PolyList& L)
    {
      if (L.empty()) return true;
      if (GradingDim(owner(L[0]))==0) return true;
      return IsHomog(L);
    }

    bool IsHomogGrD0(const VectorList& L)
    {
      if (L.empty()) return true;
      if (GradingDim(RingOf(owner(L[0])))==0) return true;
      return IsHomog(L);
    }

    bool IsHomogGrD0(ConstRefRingElem f)
    {
      if (GradingDim(owner(f))==0) return true;
      return IsHomog(f);
    }

  }
      

// GBasis

  void ComputeGBasis(VectorList& outGB, VectorList& outMinGens, const VectorList& inGens)
  {
    if (inGens.empty())
    {
      outGB.clear();
      outMinGens.clear();
      return;
    }
    const FreeModule FM=owner(inGens);
    const SparsePolyRing NewP(MakeNewPRingFromModule(FM));
    const SparsePolyRing OldP(RingOf(FM));
    if (!IsField(CoeffRing(OldP)))
      CoCoA_ERROR("coefficients are not in a field", "ComputeGBasis");
   bool IsSatAlg=false;
   const GRingInfo GRI(NewP,OldP,FM,FM,IsHomogGrD0(inGens),
                       IsSatAlg,NewDivMaskEvenPowers());
   GPolyList EmbeddedPolys=EmbedVectorList(inGens,GRI);
   GReductor GBR(GRI, EmbeddedPolys, GReductor::ourDefaultStatLevel);
   GBR.myDoGBasis();// homog input standard alg interred
   VectorList TmpGB;
   VectorList TmpMinGens;
   GBR.myGBasis(TmpGB);
   if (GradingDim(OldP)>0 && IsHomog(inGens)) GBR.myMinGens(TmpMinGens);
   outGB.clear(); // just to remember to clean this up
   outMinGens.clear(); // just to remember to clean this up
//    outGB =      DeEmbedPolyList(TmpGB, GRI);
//    outMinGens = DeEmbedPolyList(TmpMinGens, GRI);
   outGB =      TmpGB;
   outMinGens = TmpMinGens;
 }//ComputeGBasis


  void ComputeGBasis(PolyList& outGB, PolyList& outMinGens, const PolyList& inPL, const int StatLevel)
  {
    Stats S(len(inPL), StatLevel);  // AMB totally useless....
    ComputeGBasis(outGB, outMinGens, S, inPL, StatLevel);
  }
  
  void ComputeSATGBasis(PolyList& outGB, const PolyList& inGens, const int StatLevel)
  {
    Stats S(len(inGens), StatLevel);  // AMB totally useless....
    ComputeSATGBasis(outGB, S, inGens, StatLevel);
  }
  
  void ComputeSATMixGBasis(PolyList& theGB, const PolyList& inGens, const int StatLevel)
  {
    Stats S(len(inGens), StatLevel);  // AMB totally useless....
    ComputeSATMixGBasis(theGB, S, inGens, StatLevel);
  }
  

  void ComputeGBasis(PolyList& outGB, PolyList& outMinGens, Stats& outStats, const PolyList& inGens, const int StatLevel)
  {
    if (inGens.empty())
    {
      outGB.clear();
      outMinGens.clear();
      return;
    }
    SparsePolyRing SPR(owner(inGens));
    // MAX: Adding to ComputeGBasis for modules should suffice.
    // The new ring should be created by MakeNewPRingFromModule
    // and the embeddings:  WithoutDenominators by EmbedVectorList
    //                   :  WithDenominator1Hom by DeEmbedPolyList
    if (!IsField(CoeffRing(SPR)))
      CoCoA_ERROR("coefficients are not in a field", "ComputeGBasis");
    bool IsSatAlg=false;
    if (IsFractionFieldOfGCDDomain(CoeffRing(SPR)))
    {
      const ring R = BaseRing(CoeffRing(SPR));
      SparsePolyRing Rx = NewPolyRing_DMPI(R, symbols(PPM(SPR)), ordering(PPM(SPR)));
      GRingInfo GRI(Rx, IsHomogGrD0(inGens),IsSatAlg,NewDivMaskEvenPowers());
      GRI.mySetCoeffRingType(CoeffEncoding::FrFldOfGCDDomain);
      //GReductor GBR(GRI, WithoutDenominators(inGens, Rx), GReductor::ourDefaultStatLevel);
      GReductor GBR(GRI, WithoutDenominators(inGens, Rx), StatLevel);
      GBR.myDoGBasis();// homog input standard alg interred
      PolyList TmpGB;
      PolyList TmpMinGens;
      GBR.myGBasis(TmpGB);
      if (GradingDim(SPR)>0 && IsHomog(inGens)) GBR.myMinGens(TmpMinGens);
      outGB.clear();//just to remember to clean this up
      outMinGens.clear();//just to remember to clean this up
      outGB = WithDenominator1Hom(TmpGB, SPR);
      outMinGens = WithDenominator1Hom(TmpMinGens, SPR);
      outStats = GBR.myStats();
    }
    else
    {
      GRingInfo GRI(SPR,IsHomogGrD0(inGens),IsSatAlg,NewDivMaskEvenPowers());
      GReductor GBR(GRI, inGens, StatLevel);
      GBR.myDoGBasis();// homog input standard alg interred
      GBR.myGBasis(outGB);
      outMinGens.clear();//just to remember to clean this up
      if (GradingDim(SPR)>0 && IsHomog(inGens)) GBR.myMinGens(outMinGens);
      outStats = GBR.myStats();
    }          
  }//ComputeGBasis  
  
  
  void ComputeSATMixGBasis(PolyList& theGB, Stats& theStats, const PolyList& thePL,const int StatLevel)
  {      
    if (thePL.empty())
    {
      theGB.clear();
      return;
    }
    bool IsSatAlg=false;
    SparsePolyRing SPR(owner(thePL));
    // MAX: Adding to ComputeGBasis for modules should suffice.
    // The new ring should be created by MakeNewPRingFromModule
    // and the embeddings:  WithoutDenominators by EmbedVectorList
    //                   :  WithDenominator1Hom by DeEmbedPolyList
    if (!IsField(CoeffRing(SPR)))
      CoCoA_ERROR("coefficients are not in a field", "ComputeGBasis");
    if (IsFractionFieldOfGCDDomain(CoeffRing(SPR)))
    {
      const ring R = BaseRing(CoeffRing(SPR));
      SparsePolyRing Rx = NewPolyRing_DMPI(R, symbols(PPM(SPR)), ordering(PPM(SPR)));
      GRingInfo GRI(Rx, IsHomogGrD0(thePL),IsSatAlg,NewDivMaskEvenPowers());
      GRI.mySetCoeffRingType(CoeffEncoding::FrFldOfGCDDomain);
      //GReductor GBR(GRI, WithoutDenominators(thePL, Rx), GReductor::ourDefaultStatLevel);
      GReductor GBR(GRI, WithoutDenominators(thePL, Rx), StatLevel);
      GBR.myDoGBasis();// homog input standard alg interred
      PolyList TmpGB;
      GBR.myGBasis(TmpGB);
      theGB.clear();//just to remember to clean this up
      theGB=WithDenominator1Hom(TmpGB, SPR);
      theStats = GBR.myStats();
    }
    else
    {
      GRingInfo GRI(SPR,IsHomogGrD0(thePL),IsSatAlg,NewDivMaskEvenPowers());
      GReductor GBR(GRI, thePL, StatLevel);
      GBR.myDoGBasis();// homog input standard alg interred
      GBR.myGBasis(theGB);
      theStats = GBR.myStats();
    }
  }//ComputeSATMixGBasis


  void ComputeSATGBasis(PolyList& theGB, Stats& theStats, const PolyList& thePL,const int StatLevel)
  {
    if (thePL.empty())
    {
      theGB.clear();
      return;
    }
    bool IsSatAlg=true;
    SparsePolyRing SPR(owner(thePL));
    // MAX: Adding to ComputeGBasis for modules should suffice.
    // The new ring should be created by MakeNewPRingFromModule
    // and the embeddings:  WithoutDenominators by EmbedVectorList
    //                   :  WithDenominator1Hom by DeEmbedPolyList
    if (!IsField(CoeffRing(SPR)))
      CoCoA_ERROR("coefficients are not in a field", "ComputeSATGBasis");
    if (IsFractionFieldOfGCDDomain(CoeffRing(SPR)))
    {
      const ring R = BaseRing(CoeffRing(SPR));
      SparsePolyRing Rx = NewPolyRing_DMPI(R, symbols(PPM(SPR)), ordering(PPM(SPR)));
      GRingInfo GRI(Rx, IsHomogGrD0(thePL),IsSatAlg,NewDivMaskEvenPowers());
      GRI.mySetCoeffRingType(CoeffEncoding::FrFldOfGCDDomain);
      //GReductor GBR(GRI, WithoutDenominators(thePL, Rx), GReductor::ourDefaultStatLevel);
      GReductor GBR(GRI, WithoutDenominators(thePL, Rx), StatLevel,
                    GReductor::SaturatingAlg);
      GBR.myDoSATGBasis();// homog input standard alg interred
      PolyList TmpGB;
      GBR.myGBasis(TmpGB);
      theGB.clear();//just to remember to clean this up
      theGB=WithDenominator1Hom(TmpGB, SPR);
      theStats = GBR.myStats();
    }
    else
    {
      GRingInfo GRI(SPR,IsHomogGrD0(thePL),IsSatAlg,NewDivMaskEvenPowers());
      GReductor GBR(GRI, thePL, StatLevel, GReductor::SaturatingAlg);
      GBR.myDoSATGBasis();// homog input standard alg interred
      GBR.myGBasis(theGB);
      theStats = GBR.myStats();
    }
  }//ComputeSATGBasis

/*
void ComputeGBasisFrameWork(PolyList& theGB, const PolyList& theInputPolyList)

 {
   theGB.clear();
   if (theInputPolyList.empty())
     return;
   const SparsePolyRing SPR(owner(theInputPolyList));
   GRingInfo GRI(SPR,IsHomogGrD0(theInputPolyList),NewDivMaskEvenPowers());
   GReductor GBR(GRI, theInputPolyList,GReductor::ourDefaultStatLevel);


   if (false)
   {
     GBR.myDoGBasis();// homog input standard alg interred
     GBR.myGBasis(theGB);
   }


   if (false)
   {
     GBR.myPrepareGBasis();
     while (GBR.myPairsLen()!=0)
     {
       GBR.myDoGBasis(1);
     }
     GBR.myGBasis(theGB);
   }

  if (false)
  {
   //non null Spoly by non null Spoly
   GBR.myPrepareGBasis();
   //GPoly SP(GRI);
   while (GBR.myPairsLen()!=0)
   {
     GBR.myReduceUntilNonZeroRedSP();
     //SP=GBR.GetSPoly();
     //.......Things.......
     //PolyList PL;
     //GBR.GetCandidateGBasis(PL);
     //GBR.SetSPoly(SP);
     if (GBR.myPairsLen()!=0)
     {
       myTrueReductors.interreduce(mySPoly);
       myUpdateBasisAndPairs();
       if (myCurrentPairDeg!=myOldDeg&&myTrueReductors.IhaveBorelReductors())
         myTrueReductors.myBorelReductorsUpdateInNextDegree(myCurrentPairDeg);
     }
    }
   GBR.myGBasis(theGB);
  }

 }//ComputeGBasisFrameWork
*/



  void ComputeElim(PolyList& theElimResult, const PolyList& thePL, ConstRefPPMonoidElem inds)
  {
    if (thePL.empty())
    {
      theElimResult.clear();
      return;
    }
    PolyList ElimResult;
    PolyList NewPL;
    const SparsePolyRing OldSPR=owner(thePL);
    std::vector<long> IndexList=PPMonoidElem2IndexList(inds);
    //    bool IsHomogGrD0PL = (GradingDim(OldSPR)==0 ? false : IsHomogGrD0(thePL));
    bool IsHomogGrD0PL = false;
    if  (GradingDim(OldSPR)!=0 && IsHomogGrD0(thePL)) IsHomogGrD0PL = true;
    SparsePolyRing ElimSPR=MakeElimRingFromOld(OldSPR,IndexList, IsHomogGrD0PL);
    std::vector<RingElem> NewImages=indets(ElimSPR);
    std::vector<RingElem> OldImages=indets(OldSPR);
    RingHom Old2New=PolyRingHom(OldSPR,ElimSPR,CoeffEmbeddingHom(ElimSPR),NewImages);
    RingHom New2Old=PolyRingHom(ElimSPR,OldSPR,CoeffEmbeddingHom(OldSPR),OldImages);
    for (PolyList::const_iterator it=thePL.begin();it!=thePL.end();++it)
      NewPL.push_back(Old2New(*it));
    PPMonoidElem NewInds(LPP(Old2New(monomial(OldSPR,1,inds))));
    PolyList GB;
    PolyList MinGens; // useless
    ComputeGBasis(GB, MinGens, NewPL, GReductor::ourDefaultStatLevel);
    for (PolyList::iterator it=GB.begin();it!=GB.end();++it)
      if (IsCoprime(LPP(*it),NewInds))
        ElimResult.push_back(New2Old(*it));
    swap(theElimResult,ElimResult);
    return;
  }//ComputeElim

  void ComputeElim(VectorList& /*theElimResult*/,
                   const VectorList& /*theVL*/,
                   ConstRefPPMonoidElem /*inds*/)
  {
    // Remember to check if ordering is TOPOS ordering here
    // if not, use hom
    CoCoA_ERROR(ERR::NYI, "ComputeElim(const VectorList& theVL,ConstRefPPMonoidElem)");
    return;
  }//ComputeElim


  void ComputeSyz(VectorList& theSyzResult, const FreeModule& SyzFM, const VectorList& theVL)
  {
    if (!IsHomogGrD0(theVL))
      CoCoA_ERROR("Not Yet Tested for non-homogeneous input", "ComputeSyz(const VectorList&)");
    if (theVL.empty())
    {
      theSyzResult.clear();
      return;
    }
    bool IsSatAlg=false;
    VectorList SyzResult;
    ModOrdTypeForcing MOType;
    if (IsHomogGrD0(theVL))
      MOType=WDegPosTO;
    else
      MOType=PosWDegTO;
    const FreeModule FM=owner(theVL);
    const SparsePolyRing NewP(MakeNewPRingFromModule(FM,MOType));
    const SparsePolyRing OldP(RingOf(FM));
    // Note: the GRI should build itself SyzFM and NewP from the data and deduce FM and OldP.
    //       All the embedding/deembedding functions should be memebers of GRI.
    GRingInfo GRI(NewP,OldP,FM,SyzFM,IsHomogGrD0(theVL),IsSatAlg,NewDivMaskEvenPowers());
    GPolyList EmbeddedPolys=SyzEmbedVectorList(theVL,GRI);
    GReductor GBR(GRI, EmbeddedPolys, GReductor::ourDefaultStatLevel);
    GBR.myDoGBasis();
    PolyList OutputPolyList;
    GBR.myGBasis(OutputPolyList);
    SyzResult=DeEmbedPolyList(OutputPolyList,GRI,NumCompts(FM));
    swap(theSyzResult,SyzResult);
    return;
  }//ComputeSyz


  void ComputeSyz(VectorList& theSyzResult, const FreeModule& SyzFM, const PolyList& thePL)
  {
    //if (!IsHomogGrD0(thePL))
    //  CoCoA_ERROR("Not Yet Tested for non-homogeneous input", "ComputeSyz(const VectorList&)");
    if (thePL.empty())
    {
      theSyzResult.clear();
      return;
    }
    bool IsSatAlg=false;
    VectorList SyzResult;
    ModOrdTypeForcing MOType;
    if (IsHomogGrD0(thePL))
      MOType=WDegPosTO;
    else
      MOType=PosWDegTO;
    const SparsePolyRing OldP(owner(thePL));
    const SparsePolyRing NewP(MakeNewPRingForSimpleEmbedding(OldP,MOType));
    GRingInfo GRI(NewP,OldP,SyzFM,IsHomogGrD0(thePL),IsSatAlg,NewDivMaskEvenPowers());
    GPolyList EmbeddedPolys=SyzEmbedPolyList(thePL,GRI);
    GReductor GBR(GRI, EmbeddedPolys, GReductor::ourDefaultStatLevel);
    GBR.myDoGBasis();
    PolyList OutputPolyList;
    GBR.myGBasis(OutputPolyList);
    SyzResult=DeEmbedPolyList(OutputPolyList,GRI,1);
    swap(theSyzResult,SyzResult);
    return;
  }//ComputeSyz


  void ComputeIntersection(VectorList& theIntersectionResult,
                           const VectorList& theVL1,
                           const VectorList& theVL2)
  {
    bool IsSatAlg=false;
    VectorList IntersectionResult;
    const FreeModule FM=owner(theVL1);
    //const SparsePolyRing NewP(MakeNewPRingFromModule(FM,WDegPosTO));
    const SparsePolyRing NewP(MakeNewPRingFromModulePosFirst(FM,IsHomogGrD0(theVL1)&&IsHomogGrD0(theVL2)));
    const SparsePolyRing OldP(RingOf(FM));
    GRingInfo GRI(NewP,OldP,FM,FM,IsHomogGrD0(theVL1)&&IsHomogGrD0(theVL2),
                  IsSatAlg,NewDivMaskEvenPowers());
    GPolyList EmbeddedPolys=IntEmbedVectorLists(theVL1,theVL2,GRI);
    GReductor GBR(GRI, EmbeddedPolys, GReductor::ourDefaultStatLevel);
    GBR.myDoGBasis();// homog input standard alg interred
    PolyList OutputPolyList;
    GBR.myGBasis(OutputPolyList);
    IntersectionResult=DeEmbedPolyList(OutputPolyList,GRI,NumCompts(FM));
    swap(IntersectionResult,theIntersectionResult);
    return;
  }//ComputeIntersection


  void ComputeIntersection(PolyList& theIntersectionResult,
                           const PolyList& thePL1,
                           const PolyList& thePL2)
  {
    bool IsSatAlg=false;
    PolyList IntersectionResult;
    const SparsePolyRing OldP(owner(thePL1));
    const SparsePolyRing NewP(MakeNewPRingForSimpleEmbeddingPosFirst(OldP,IsHomogGrD0(thePL1)&&IsHomogGrD0(thePL2)));
    GRingInfo GRI(NewP,OldP,IsHomogGrD0(thePL1)&&IsHomogGrD0(thePL2),IsSatAlg,
                  NewDivMaskEvenPowers());
    GPolyList EmbeddedPolys=IntEmbedPolyLists(thePL1,thePL2,GRI);
    GReductor GBR(GRI, EmbeddedPolys, GReductor::ourDefaultStatLevel);
    GBR.myDoGBasis();// homog input standard alg interred
    PolyList OutputPolyList;
    GBR.myGBasis(OutputPolyList);
    IntersectionResult=DeEmbedPolyListToPL(OutputPolyList,GRI,1);
    swap(theIntersectionResult,IntersectionResult);
    return;
  }//ComputeIntersection


// Colon by a single vector
  void ComputeColonByPrincipal(PolyList& theColonResult, const VectorList& theVL1, const VectorList& theVL2)
  {
    bool IsSatAlg=false;
    //    for (VectorList::const_iterator it=theVL1.begin();it!=theVL1.end();++it)
    //      if (!IsHomogGrD0(*it))
      if (!IsHomogGrD0(theVL1))
        CoCoA_ERROR(ERR::NYI, "ComputeColonByPrincipal(const VectorList&, const VectorList&) Not Yet Implemented for non-homogeneous input");
      //    for (VectorList::const_iterator it=theVL2.begin();it!=theVL2.end();++it)
      //      if (!IsHomogGrD0(*it))
      if (!IsHomogGrD0(theVL2))
        CoCoA_ERROR(ERR::NYI, "ComputeColonByPrincipal(const VectorList&, const VectorList&) Not Yet Implemented for non-homogeneous input");
    PolyList ColonResult;
    const FreeModule FM=owner(theVL1);
    const SparsePolyRing NewP(MakeNewPRingFromModule(FM,WDegPosTO));
    const SparsePolyRing OldP(RingOf(FM));
    GRingInfo GRI(NewP,OldP,FM,FM,IsHomogGrD0(theVL1)&&IsHomogGrD0(theVL2),
                  IsSatAlg,NewDivMaskEvenPowers());
    GPolyList EmbeddedPolys=ColonEmbedVectorLists(theVL1,theVL2,GRI);
    GReductor GBR(GRI, EmbeddedPolys, GReductor::ourDefaultStatLevel);
    GBR.myDoGBasis();// homog input standard alg interred
    PolyList OutputPolyList;
    GBR.myGBasis(OutputPolyList);
    ColonResult=DeEmbedPolyListToPL(OutputPolyList,GRI,NumCompts(FM));
    swap(ColonResult,theColonResult);
    return;
  }//ComputeColon

// Colon
  void ComputeCColon(PolyList& theCColonResult,
                     const VectorList& theVL1,
                     const VectorList& theVL2)
  {
    if (theVL1.empty() && theVL2.empty())
      CoCoA_ERROR("One argument must be non empty","ComputeCColon");
    PolyList CColonResult;
    if (theVL2.empty())
    {
      CColonResult.push_back(one(RingOf(owner(theVL1))));
      swap(CColonResult,theCColonResult);
      return;
    }
    if (theVL1.empty())
    {
      swap(CColonResult,theCColonResult);
      return;
    }
    PolyList tmp;
    VectorList::const_iterator it=theVL2.begin();
    ComputeColonByPrincipal(CColonResult,theVL1,MakeVectorList(*it));
    ++it;
    for (;it!=theVL2.end();++it)
    {
      ComputeColonByPrincipal(tmp,theVL1,MakeVectorList(*it));
      ComputeIntersection(CColonResult,CColonResult,tmp);
    }
    swap(CColonResult,theCColonResult);
    return;
  }//ComputeCColon


  void ComputeColonByPrincipal(PolyList& theColonResult,
                    const PolyList& thePL1,
                    ConstRefRingElem f)
  {
    PolyList TmpColonResult;
    if (IsZero(f))
    {
      TmpColonResult.push_back(one(owner(thePL1)));
      swap(theColonResult, TmpColonResult);
      return;    
    }
    bool IsSatAlg=false;
    const SparsePolyRing OldP(owner(thePL1));
    const SparsePolyRing NewP(MakeNewPRingForSimpleEmbeddingPosFirst(OldP,IsHomogGrD0(thePL1)&&IsHomogGrD0(f)));
    GRingInfo GRI(NewP,OldP,IsHomogGrD0(thePL1)&&IsHomogGrD0(f),
                  IsSatAlg,NewDivMaskEvenPowers());
    GPolyList EmbeddedPolys=ColonEmbedPolyLists(thePL1, std::vector<RingElem>(1,f), GRI);
    GReductor GBR(GRI, EmbeddedPolys, GReductor::ourDefaultStatLevel);
    GBR.myDoGBasis();// homog input standard alg interred
    PolyList OutputPolyList;
    GBR.myGBasis(OutputPolyList);
    TmpColonResult = DeEmbedPolyListToPL(OutputPolyList,GRI,1);
    swap(theColonResult, TmpColonResult);
    return;
  }//ComputeColonByPrincipal


// Colon
  void ComputeCColon(PolyList& theCColonResult,
                     const PolyList& thePL1,
                     const PolyList& thePL2)
  {
    if (thePL1.empty() && thePL2.empty())
      CoCoA_ERROR("One argument must be non empty","ComputeCColon");
    PolyList CColonResult;
    if (thePL2.empty())
    {
      CColonResult.push_back(one(owner(thePL1)));
      swap(theCColonResult,CColonResult);
      return;
    }
    if (thePL1.empty())
    {
      theCColonResult.clear();
      return;
    }
    PolyList tmp;
    PolyList::const_iterator it=thePL2.begin();
    ComputeColonByPrincipal(CColonResult, thePL1, *it);
    ++it;
    for (;it!=thePL2.end();++it)
    {
      ComputeColonByPrincipal(tmp, thePL1, *it);
      ComputeIntersection(CColonResult,CColonResult,tmp);
    }
    swap(theCColonResult,CColonResult);
    return;
  }//ComputeCColon


  void ComputeColonByPrincipal(VectorList& /*theColonResult*/,
                    const VectorList& /*theVL*/,
                    const PolyList& /*thePL*/)
  {
    CoCoA_ERROR(ERR::NYI,"ComputeColonByPrincipal(VectorList, PolyList)");
    return;
  }//ComputeColonByPrincipal
  

  namespace // anonymous namespace for file local functions and definitions
  {
  // These procedures are for EqualLPPs used by Saturation
  
  // PolyList::const_iterators are ordered according to LPP of their polys
  bool ByLPP(const PolyList::const_iterator& it,
             const PolyList::const_iterator& it1)
  {
    return (SparsePolyRingPtr(owner(*it))->myCmpLPP(raw(*it),raw(*it1)) == -1);
  }


  bool AreEqualLPPs(std::vector<PolyList::const_iterator>& theV1,
                    std::vector<PolyList::const_iterator>& theV2)
  {
    const long lenV1 = len(theV1);
    const long lenV2 = len(theV2);
    if (lenV1 != lenV2)  return false;
    stable_sort(theV1.begin(), theV1.end(), ByLPP);
    stable_sort(theV2.begin(), theV2.end(), ByLPP);
    for (int i=0; i!=lenV1; ++i)
      if (LPP(*(theV1[i])) != LPP(*(theV2[i])))  return false;
    return true;
  }
    
    
//     bool AreEqualLPPs(const std::vector<VectorList::const_iterator>& theV1,
//                       const std::vector<VectorList::const_iterator>& theV2)
//     {
    // const long lenV1 = len(V1);
    // const long lenV2 = len(V2);
//       if (lenV1 != lenV2)  return false;
//       stable_sort(theV1.begin(), theV1.end(), ByLPP);
//       stable_sort(theV2.begin(), theV2.end(), ByLPP);
//       for (int i=0; i!=lenV1; ++i)
//         //if (*(theV1[i])!=*(theV2[i]))
//         // when LT exists if (LT(*(theV1[i]))!=LT(*(theV2[i])))
//         return false;
//       return true;
//     }

    // Useful when you have I\subset J and you want to check I==J
    bool AreEqualLPPs(const PolyList& thePL1, const PolyList& thePL2)
    {
      std::vector<PolyList::const_iterator> V1,V2;
      for (PolyList::const_iterator it=thePL1.begin(); it!=thePL1.end(); ++it)
        V1.push_back(it);
      for (PolyList::const_iterator it=thePL2.begin(); it!=thePL2.end(); ++it)
        V2.push_back(it);
      return AreEqualLPPs(V1,V2);
    }
    
//     // Useful when you have M\subset N and you want to check M==N
//     bool AreEqualLPPs(const VectorList& theVL1, const VectorList& theVL2)
//     {
//       std::vector<VectorList::const_iterator> V1,V2;
//       for (VectorList::const_iterator it=theVL1.begin();it!=theVL1.end();++it)
//         V1.push_back(it);
//       for (VectorList::const_iterator it=theVL2.begin();it!=theVL2.end();++it)
//         V2.push_back(it);
//       return AreEqualLPPs(V1,V2);
//     }
  
  } // anonymous namespace

   
//---------------------------------------------------------

  void ComputeSaturationByPrincipal(PolyList& theSaturationResult,
                         const PolyList& thePL1,
                         ConstRefRingElem f)
  {
    PolyList TmpA;
    ComputeColonByPrincipal(TmpA, thePL1, f);
    PolyList TmpB;
    ComputeColonByPrincipal(TmpB, TmpA, f);
    while (!AreEqualLPPs(TmpA,TmpB))
    {
      swap(TmpA, TmpB);
      ComputeColonByPrincipal(TmpB, TmpA, f);
    }
    swap(theSaturationResult, TmpA);
    return;
  }//ComputeSaturationByPrincipal


  void ComputeSSaturation(PolyList& theSSaturationResult,
                          const PolyList& thePL1,
                          const PolyList& thePL2)
  {
    if (thePL1.empty() && thePL2.empty())
      CoCoA_ERROR("One argument must be non empty","ComputeSSaturation");
    if (thePL2.empty())
    {
      theSSaturationResult.clear();// this or sawp? this look better
      theSSaturationResult.push_back(one(owner(thePL1)));
      return;
    }
    if (thePL1.empty())
    {
      theSSaturationResult.clear();
      return;
    }

    PolyList TmpPL1;
    ComputeCColon(TmpPL1,thePL1,thePL2);
    PolyList TmpPL2;
    ComputeCColon(TmpPL2,TmpPL1,thePL2);
    while (!AreEqualLPPs(TmpPL1,TmpPL2))
    {
      swap(TmpPL1,TmpPL2);
      ComputeCColon(TmpPL2,TmpPL1,thePL2);
    }
    swap(theSSaturationResult,TmpPL2);
    return;
  }//ComputeSSaturation


  void ComputeSaturationByPrincipal(VectorList& /*theSaturation*/,
                         const VectorList& /*theVL*/,
                         const PolyList& /*thePL*/)
  {
    CoCoA_ERROR(ERR::NYI,"ComputeSaturationByPrincipal(VectorList, PolyList)");
    return;
  }//ComputeSaturationByPrincipal


  void ComputeSSaturation(VectorList& /*theSSaturation*/,
                          const VectorList& /*theVL*/,
                          const PolyList& /*thePL*/)
  {
    CoCoA_ERROR(ERR::NYI,"ComputeSSaturation(VectorList, PolyList)");
    return;
  }//ComputeSSaturation


  void ComputeHomogenization(VectorList& /*theHomogResult*/,
                             const VectorList& /*theVL*/,
                             const PolyList& /*thePL*/)
  {
    CoCoA_ERROR(ERR::NYI,"ComputeHomogenization(VectorList, PolyList)");
    return;
  }//ComputeHomogenization


  void ComputeHomogenization(PolyList& theHomogResult,
                             const PolyList& thePL1,
                             const PolyList& theIndets)
  {
    if (thePL1.empty())
    {
      theHomogResult.clear();
      return;
    }
    PolyList HomogResult;
    const SparsePolyRing SPR=owner(thePL1);
    RingElem IndetProduct(product(theIndets));
    RingElem tmp(SPR);
    for (PolyList::const_iterator it=thePL1.begin();it!=thePL1.end();++it)
      HomogResult.push_back(homog(*it, theIndets));
    ComputeSaturationByPrincipal(HomogResult,HomogResult,IndetProduct);
    swap(theHomogResult,HomogResult);
    return;
  }//ComputeHomogenization


// WARN: it supposes ComputeSaturationByPrincipal returns a GB
  bool RadicalMembership(const PolyList& PL, ConstRefRingElem the_f)
  {
    PolyList PL1;
    ComputeSaturationByPrincipal(PL1, PL, the_f);
    if (len(PL) != 1) return false;
    //    monic(PL1);
    //    if (IsOne(PL1.front()))
    //      return true;
    //    else
    //      return false;
    return IsInvertible(PL1.front());
  }//RadicalMembership

  void ComputeLT(VectorList& /*theLTResult*/,const VectorList& /*theVL*/)
  {
    // Waiting for LT of a ModuleElem
    CoCoA_ERROR(ERR::NYI,"ComputeLT(VectorList)");
    return;
  }//ComputeLT


  void ComputeLT(PolyList& theLTResult, const PolyList& thePL)
  {
    PolyList GB;
    PolyList MinGens; // useless
    ComputeGBasis(GB, MinGens, thePL, GReductor::ourDefaultStatLevel);
    if (GB.empty())
    {
      swap(theLTResult,GB);
      return;
    }
    PolyList L;
    SparsePolyRing P(owner(*(GB.begin())));
    for (PolyList::const_iterator it=GB.begin();it!=GB.end();++it)
    {
      SparsePolyIter it_p=BeginIter(*it);
      L.push_back(monomial(P, one(CoeffRing(P)), PP(it_p)));
    }
    swap(theLTResult,L);
  }//ComputeLT

}// end namespace cocoa


// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/TmpGOperations.C,v 1.31 2014/07/31 14:45:18 abbott Exp $
// $Log: TmpGOperations.C,v $
// Revision 1.31  2014/07/31 14:45:18  abbott
// Summary: Merged io.H and UtilsTemplate.H into new header VectorOperations.H
// Author: JAA
//
// Revision 1.30  2014/07/30 14:10:18  abbott
// Summary: Changed BaseRing into RingOf
// Author: JAA
//
// Revision 1.29  2014/07/14 15:07:57  abbott
// Summary: Changed include of tmp.H into UtilsTemplate.H
// Author: JAA
//
// Revision 1.28  2014/07/09 14:27:53  abbott
// Summary: Removed AsFreeModule and AsFGModule
// Author: JAA
//
// Revision 1.27  2014/07/08 08:39:14  abbott
// Summary: Removed AsFractionField
// Author: JAA
//
// Revision 1.26  2014/07/07 13:02:17  abbott
// Summary: Removed AsSparsePolyRing
// Author: JAA
//
// Revision 1.25  2014/04/30 16:14:38  abbott
// Summary: Replaced X.size() by len(X)
// Author: JAA
//
// Revision 1.24  2014/03/26 15:23:44  bigatti
// -- added MinGens for submodules
//
// Revision 1.23  2014/03/21 13:08:09  bigatti
// -- cosmetics
//
// Revision 1.22  2013/10/30 09:55:56  bigatti
// -- dealing with 0 generator in ComputeColonByPrincipal
//
// Revision 1.21  2013/06/12 08:54:07  bigatti
// -- added computation of MinGens (in ComputeGBasis)
// -- changed some "the" in "out"/"in" in argument names
//
// Revision 1.20  2013/02/21 17:16:45  bigatti
// -- changed syntax for ComputeSyz
//
// Revision 1.19  2013/02/12 16:36:49  bigatti
// -- added temporary (ehm..) function IsHomogGrD0 for GradingDim==0
//
// Revision 1.18  2013/01/31 11:43:22  bigatti
// -- added Stats field to ComputeXXXGBasis for returning statistics
//
// Revision 1.17  2012/10/03 12:24:46  bigatti
// -- some (..)ByPrincipal now take a RingElem instead of PolyList
// -- some cleaning
//
// Revision 1.16  2012/10/02 16:44:41  bigatti
// -- just aesthetics
//
// Revision 1.15  2012/05/28 09:18:20  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.14  2012/05/22 10:02:37  abbott
// Removed IsGCDDomain; substituted by IsTrueGCDDomain.
// Added IsFractionFieldOfGCDDomain.
//
// Revision 1.13  2011/12/07 15:54:34  bigatti
// -- renamed ambiguous "operator<" and hidden into anonymous namespace
// -- renamed ambiguous "operator==" into AreEqualLPPs (used by Saturation)
//
// Revision 1.12  2011/12/05 16:32:10  bigatti
// -- fixed bug about saturation (by non-principal ideal)
//
// Revision 1.11  2010/07/27 07:37:13  bigatti
// -- new class GBCriteria, simplified GReductor ctor
//
// Revision 1.10  2009/04/27 13:18:21  bigatti
// -- added SaturatingAlg flag
//
// Revision 1.9  2009/01/30 13:41:50  bigatti
// -- enum instead of bool arguments
//
// Revision 1.8  2008/09/19 13:33:42  bigatti
// -- added: Sat algorithm (M.Caboara)
//
// Revision 1.7  2008/09/19 11:34:15  bigatti
// -- new mechanism for passing verbosity level (or StatLevel)
//    [only partially tested]
//
// Revision 1.6  2008/07/09 16:09:47  abbott
// Cosmetic tidying.  Removed bogus comments. Added missing & (for C++ const reference).
//
// Revision 1.5  2008/03/12 16:35:18  bigatti
// -- changed: IsHomogeneous --> IsHomog
// -- changed: ERR:ZeroPoly --> ERR::ZeroRingElem
//
// Revision 1.4  2007/11/09 10:45:52  bigatti
// -- [caboara] preparation for self-saturating algorithm
//
// Revision 1.3  2007/10/30 17:14:07  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.2  2007/03/27 13:19:21  bigatti
// -- removed printout
//
// Revision 1.1  2007/03/09 18:56:56  bigatti
// -- added Tmp prefix to Groebner related files
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.17  2007/03/08 18:22:29  cocoa
// Just whitespace cleaning.
//
// Revision 1.16  2007/03/07 17:04:31  cocoa
// -- several changes by M.Caboara: more operations on ideals,
//    exception cleaner, coding conventions, WSugar, dynamic
//
// Revision 1.15  2007/03/05 21:06:07  cocoa
// New names for homomorphism pseudo-ctors: removed the "New" prefix.
//
// Revision 1.14  2007/02/10 18:44:03  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.13  2007/01/25 12:05:20  bigatti
// -- fixed: removed coefficients in LT output
//
// Revision 1.12  2006/12/21 13:48:33  cocoa
// Made all increment/decrement calls prefix (except where the must be postfix).
//
// Revision 1.11  2006/12/04 13:55:54  cocoa
// -- added: sugar for GradingDim > 0  (called wsugar)
//
// Revision 1.10  2006/11/24 17:01:42  cocoa
// -- reorganized includes of header files
//
// Revision 1.9  2006/11/22 14:43:32  cocoa
// -- minor cleaning (indicated by Intel compiler)
//
// Revision 1.8  2006/11/20 14:57:17  cocoa
// -- added: (standard) sugar for modules
// -- fixed: non-homogeneous sysygies
// -- minor fixes     [M.Caboara]
//
// Revision 1.7  2006/10/06 16:46:17  cocoa
// -- syzygies for non-homogenous polynomials (Max)
// -- wip: evolution of Groebner Framework (Max)
//
// Revision 1.6  2006/10/06 10:48:34  cocoa
// Removed #include references to GradedFreeModule.H
//
// Revision 1.5  2006/08/17 09:35:49  cocoa
// -- changed: checks for homogeneity of input
//
// Revision 1.4  2006/08/07 21:23:25  cocoa
// Removed almost all publicly visible references to SmallExponent_t;
// changed to long in all PPMonoid functions and SparsePolyRing functions.
// DivMask remains to sorted out.
//
// Revision 1.3  2006/07/18 10:52:24  cocoa
// -- added check for non-homogeneous input in ComputeElim
//
// Revision 1.2  2006/06/12 13:36:35  cocoa
// -- fixed embarrassing bug in WithDenominator1Hom clearing denominators
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.17  2006/05/16 08:59:16  cocoa
// -- added function for interactive Groebner
//
// Revision 1.16  2006/05/12 17:03:16  cocoa
// -- swapped arguments in homogenized
//
// Revision 1.15  2006/05/11 16:00:22  cocoa
// -- fixed spelling of "homogenize"
//
// Revision 1.14  2006/04/27 14:01:11  cocoa
// -- tidied up include files (using GTypes.H)
//
// Revision 1.13  2006/04/26 10:02:12  cocoa
// -- added GReductor::ourDefaultStatLevel variable to allow CoCoAServer
//    to set statistics level
//
// Revision 1.12  2006/04/21 16:47:06  cocoa
// -- new syntax for ComputeGBasis by Max
//
// Revision 1.11  2006/04/11 16:43:16  cocoa
// -- changed call of LPPNoCopy --> LPP (does not make a copy)
//
// Revision 1.10  2006/04/11 16:22:40  cocoa
// -- added: Elim, LT
//
