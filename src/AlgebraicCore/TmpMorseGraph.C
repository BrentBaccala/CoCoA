//   Copyright (c)  2013  Mario Albert

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


#include "CoCoA/TmpMorseGraph.H"

#include "CoCoA/DenseMatrix.H"
#include "CoCoA/MatrixOperations.H"
#include "CoCoA/MatrixView.H"
#include "CoCoA/TmpResolutionMinimization.H"

using std::make_pair;

namespace CoCoA
{

  std::vector<MorseElement> MorseGraph::myComputeGeneralBasis() const
  {
    std::vector<MorseElement> basis;
    for (PommBasisElem it = myBasis.begin(); it != myBasis.end(); ++it)
    {
      const std::vector<long> NonMultVars(myDynamicBitsetToLong(flip(it->second)));
      const std::vector<DynamicBitset> bitsets(myPossibleWedges(NonMultVars, NumIndets(myRing)));
      for (std::vector<DynamicBitset>::const_iterator BitsetIter = bitsets.begin(); BitsetIter != bitsets.end(); ++BitsetIter)
      {
        basis.push_back(MorseElement(*BitsetIter, it));
      }
    }
    return basis;
  }


  std::vector<DynamicBitset> MorseGraph::myPossibleWedges(const std::vector<long>& NonMultVars, long length) const
  {
    std::vector<std::vector<long> > ResultAsVector;
    const long NumNonMultVars = len(NonMultVars);
    for (long i = 0; i <= NumNonMultVars; ++i)
    {
      myVariationWithoutRepetition(ResultAsVector, std::vector<long>(), NonMultVars, i);
    }
    std::vector<DynamicBitset> result; result.reserve(len(ResultAsVector));
    for (std::vector<std::vector<long> >::iterator it = ResultAsVector.begin(); it != ResultAsVector.end(); ++it)
    {
      result.push_back(myVectorLongToDynamicBitset(*it, length));
    }
    return result;
  }


  void MorseGraph::myVariationWithoutRepetition(std::vector<std::vector<long> >& result,
                                                const std::vector<long>& CurrentResult,
                                                const std::vector<long>& InputSet,
                                                long length)  const
  {
    if (length == 0)
    {
      result.push_back(CurrentResult);
      return;
    }
    for (std::vector<long>::const_iterator it = InputSet.begin(); it != InputSet.end(); ++it)
    {
      std::vector<long> TmpCurrentResult(CurrentResult);
      TmpCurrentResult.push_back(*it);
      std::vector<long>::const_iterator iter(it);
      ++iter;
      std::vector<long> TmpInputSet(iter, InputSet.end());
      myVariationWithoutRepetition(result, TmpCurrentResult, TmpInputSet, length - 1);
    }
  }


  void MorseGraph::myComputeBasicGraph(const std::vector<MorseElement>& elems, StandardRepresentationContainer& container)
  {
    for (std::vector<MorseElement>::const_iterator it = elems.begin(); it != elems.end(); ++it)
    {
      ConstResIter IterToRes(myResolution.insert(make_pair(*it, MorsePaths())).first);
      std::vector<std::pair<MorseElement, RingElem> > maps(it->myComputeBasicMaps(myGetBasisRange(), container));
      myAddMapsToResolution(maps, IterToRes);
    }
  }


  void MorseGraph::myComputeBasicConstantGraph(const std::vector<MorseElement>& elems, StandardRepresentationContainer& container)
  {
    for (std::vector<MorseElement>::const_iterator it = elems.begin(); it != elems.end(); ++it)
    {
      ConstResIter IterToRes(myResolution.insert(make_pair(*it, MorsePaths())).first);
      std::vector<std::pair<MorseElement, RingElem> > maps(it->myComputeBasicConstantMaps(myGetBasisRange(), container));
      myAddMapsToResolution(maps, IterToRes);
    }
  }


  void MorseGraph::myAddMapsToResolution(const std::vector<std::pair<MorseElement, RingElem> >& maps, const ConstResIter& origin)
  {
    for (std::vector<std::pair<MorseElement, RingElem> >::const_iterator MapIter = maps.begin(); MapIter != maps.end(); ++MapIter)
    {
      myResolution[MapIter->first].myAddPath(origin, MapIter->second);
    }
  }



  void MorseGraph::myDirectMorseReduction(StandardRepresentationContainer& container)
  {
    std::map<MorseElement, MorsePaths>::iterator ResolutionIter(myResolution.end());
    --ResolutionIter;
    while (ResolutionIter != myResolution.begin())
    {
      if ((ResolutionIter->first).IAmBasisElement())
      {
        --ResolutionIter;
        continue;
      }
      MorseElement morse(ResolutionIter->first);
      //assert !(ResolutionIter->second).IAmEmpty() 
      const long maximum(morse.myMaxTypeOne());
      if ((maximum == -1) || (maximum != morse.myMaxTypeTwo()) || (ResolutionIter->second).IamEmpty())
      {
        ResolutionIter = myDeleteAndJumpToPrevious(ResolutionIter);
        continue;
      }
      const PathMap& paths((ResolutionIter->second).myGetPaths());
      std::vector<long> longs(morse.myGetWedgeProductAsLongs());
      // because of myMaxTypeOne i is not in longs
      longs.push_back(maximum);
      for (std::vector<long>::iterator LongIter = longs.begin(); LongIter != longs.end(); ++LongIter)
      {
        DynamicBitset NewWedgeProduct(morse.myGetWedgeProduct());
        NewWedgeProduct.mySet(maximum, true);
        NewWedgeProduct.mySet(*LongIter, false);
        myLeftMinimization(paths, NewWedgeProduct, morse, maximum, *LongIter);
        myRightMinimization(paths, NewWedgeProduct, morse, maximum, *LongIter, container);
      }
      ResolutionIter = myDeleteAndJumpToPrevious(ResolutionIter);
    }
  }


  void MorseGraph::myConstantDirectMorseReduction(StandardRepresentationContainer& container)
  {
    std::map<MorseElement, MorsePaths>::iterator ResolutionIter(myResolution.end());
    --ResolutionIter;
    while (ResolutionIter != myResolution.begin())
    {
      if ((ResolutionIter->first).IAmBasisElement())
      {
        --ResolutionIter;
        continue;
      }
      MorseElement morse(ResolutionIter->first);
      //assert !(ResolutionIter->second).IAmEmpty() 
      const long maximum(morse.myMaxTypeOne());
      if ((maximum == -1) || (maximum != morse.myMaxTypeTwo()) || (ResolutionIter->second).IamEmpty())
      {
        ResolutionIter = myDeleteAndJumpToPrevious(ResolutionIter);
        continue;
      }
      const PathMap& paths((ResolutionIter->second).myGetPaths());
      std::vector<long> longs(morse.myGetWedgeProductAsLongs());
      // because of myMaxTypeOne i is not in longs
      longs.push_back(maximum);
      for (std::vector<long>::iterator LongIter = longs.begin(); LongIter != longs.end(); ++LongIter)
      {
        DynamicBitset NewWedgeProduct(morse.myGetWedgeProduct());
        NewWedgeProduct.mySet(maximum, true);
        NewWedgeProduct.mySet(*LongIter, false);
        myRightMinimization(paths, NewWedgeProduct, morse, maximum, *LongIter, container);
      }
      ResolutionIter = myDeleteAndJumpToPrevious(ResolutionIter);
    }
  }


  std::map<MorseElement, MorsePaths>::iterator MorseGraph::myDeleteAndJumpToPrevious(std::map<MorseElement, MorsePaths>::iterator iter)
  {
    std::map<MorseElement, MorsePaths>::iterator TmpIter(iter);
    --iter;
    myResolution.erase(TmpIter);
    return iter;
  }

  RingElem MorseGraph::myCreateNewBasisElement(const MorseElement& m, long IndexMult, long IndexDiv) const
  {
    RingElem result(indet(myRing, IndexMult));
    result *= monomial(myRing, one(CoeffRing(myRing)), m.myGetRightFactor());
    result = result / indet(myRing, IndexDiv);
    result *= m.myGetBasisElement();
    return result;
  }

  void MorseGraph::myLeftMinimization(const PathMap& paths, DynamicBitset NewWedgeProduct, const MorseElement& origin, long maximum, long LongIter)
  {
    MorseElement NewMorseElement(origin);
    NewMorseElement.mySetWedgeProduct(NewWedgeProduct);
    NewMorseElement.myDivideRightProductWith(maximum);

    for (PathMap::const_iterator PathIter = paths.begin(); PathIter != paths.end(); ++PathIter)
    {
      myResolution[NewMorseElement].myAddPath(PathIter->first, PathIter->second * origin.myEpsilon(maximum, maximum) * origin.myEpsilon(LongIter, maximum) * indet(myRing, LongIter));
    }
  }


  void MorseGraph::myRightMinimization(const PathMap& paths,  DynamicBitset NewWedgeProduct, const MorseElement& origin, long maximum, long LongIter, StandardRepresentationContainer& container)
  {

    if (LongIter != maximum && (!origin.IamMultiplicativeIndex(LongIter)))
    {
      const bool InPolyRing(myMapRing == myRing);
      RingElem NewBasis(myCreateNewBasisElement(origin, LongIter, maximum));
      std::vector<RingElem> NormalForm(container.myComputeStandardRepresentation(NewBasis));
      std::vector<std::pair<RingElem, DynamicBitset> >::iterator BasisIter(myBasis.begin());
      for (std::vector<RingElem>::iterator NormalFormIter = NormalForm.begin(); NormalFormIter != NormalForm.end(); ++NormalFormIter)
      {
        for (SparsePolyIter SPI = myRing->myBeginIter(raw(*NormalFormIter)); SPI != myRing->myEndIter(raw(*NormalFormIter)); ++SPI)
        {
          MorseElement NewMorseElement(NewWedgeProduct, PP(SPI), BasisIter);
          for (PathMap::const_iterator PathIter = paths.begin(); PathIter != paths.end(); ++PathIter)
          {
            RingElem path(PathIter->second);
            if (InPolyRing)
            {
              SparsePolyRingPtr(myMapRing)->myMulByCoeff(raw(path), raw(coeff(SPI)));
            } else {
              myMapRing->myMul(raw(path), raw(path), raw(coeff(SPI)));
            }
            if (origin.myEpsilon(maximum, maximum) * origin.myEpsilon(LongIter, maximum) == 1)
            {
              myMapRing->myNegate(raw(path), raw(path));
            }
            myResolution[NewMorseElement].myAddPath(PathIter->first, path);
          }
        }
        ++BasisIter;
      }
    }
  }


  void MorseGraph::myComputeResolution()
  {
    myResolution.clear();
    StandardRepresentationContainer container(myMill);
    myComputeBasicGraph(myComputeGeneralBasis(), container);
    myDirectMorseReduction(container);
  }


  matrix MorseGraph::myComputeBettiNumbers()
  {
    myComputeConstantResolution();
    std::vector< std::vector<long> > ranks(myComputeRanks());
    std::vector< std::vector<long> > WasteRanks(myComputeWasteRanks(ranks));
    ranks = myMatrixMinus(ranks, WasteRanks);
    return myTransformRanksToBettis(ranks);
  }


  std::pair<matrix, matrix> MorseGraph::myComputePseudoBettiNumbers()
  {
    myComputeConstantResolution();
    std::vector< std::vector<long> > ranks(myComputeRanks());
    std::vector< std::vector<long> > WasteRanks(myComputeWasteRanks(ranks));
    return make_pair(myTransformRanksToBettis(ranks), myTransformRanksToBettis(myMatrixMinus(ranks,WasteRanks)));
  }


  std::vector< std::vector<long> > MorseGraph::myComputeWasteRanks(const std::vector< std::vector<long> >& ranks) const
  {
    std::map<MorseElement, MorsePaths>::const_iterator iter(myResolution.end());
    --iter;
    const long size((iter->first).myCountWedgeBasis());
    std::vector< std::vector<long> > WasteRanks (len(ranks), std::vector<long>(len(ranks[0]), 0));
    std::map<MorseElement, MorsePaths>::const_iterator ResIter(myResolution.begin());
    for (long i = 0; i < size; ++i)
    {
      myComputeWasteRanksPerMap(WasteRanks, ResIter, ranks, i);
    }
    return WasteRanks;
  }


  void MorseGraph::myComputeWasteRanksPerMap(std::vector< std::vector<long> >& WasteRanks, std::map<MorseElement, MorsePaths>::const_iterator& ResIter, const std::vector< std::vector<long> >& ranks, long PosInRes) const
  {
    std::vector<std::pair<long, long> > RanksOfDegreeMatrices(myComputeWasteRanksPerDegree(ResIter, ranks[PosInRes], ranks[PosInRes + 1], PosInRes));
    for (std::vector<std::pair<long, long> >::iterator RMIter = RanksOfDegreeMatrices.begin(); RMIter != RanksOfDegreeMatrices.end(); ++RMIter)
    {
      WasteRanks[PosInRes][RMIter->first] = WasteRanks[PosInRes][RMIter->first] + RMIter->second;
      long TooMuch(0);
      if (WasteRanks[PosInRes][RMIter->first] > ranks[PosInRes][RMIter->first])
      {
        TooMuch = WasteRanks[PosInRes][RMIter->first] - ranks[PosInRes][RMIter->first];
        WasteRanks[PosInRes][RMIter->first] = ranks[PosInRes][RMIter->first];
      }
      WasteRanks[PosInRes + 1][RMIter->first] = RMIter->second - TooMuch;
    }
  }


  std::vector<std::pair<long, long> > MorseGraph::myComputeWasteRanksPerDegree(std::map<MorseElement, MorsePaths>::const_iterator& ResIter, const std::vector<long>& RowRanks, const std::vector<long>& ColRanks, long PosInRes) const
  {
    std::vector<std::pair<long, long> > RanksOfDegreeMatrices;

    std::map<std::pair<long, MorseElement>, long> identifier(myGradedPositionMorseElements(PosInRes));
    for (long degree = 0; degree < len(RowRanks); ++degree)
    {
      if (RowRanks[degree] == 0 || ColRanks[degree] == 0)
      {
        if (RowRanks[degree] != 0)
        {
          while((ResIter->first).myDegree() == degree && (ResIter->first).myCountWedgeBasis() == PosInRes)
          {
            ++ResIter;
          }
        }
        continue;
      }
      matrix SubMat(ConstructDegreeMatrix(ResIter, RowRanks[degree], ColRanks[degree], degree, PosInRes, identifier));
      RanksOfDegreeMatrices.push_back(make_pair(degree, rank(SubMat)));
    }
    return RanksOfDegreeMatrices;
  }


  matrix MorseGraph::ConstructDegreeMatrix(std::map<MorseElement, MorsePaths>::const_iterator& ResIter, long rows, long cols, long degree, long PosInRes, const std::map<std::pair<long, MorseElement>, long>& identifier) const
  {
    matrix SubMat(NewDenseMat(myMapRing, rows, cols));
    while((ResIter->first).myDegree() == degree && (ResIter->first).myCountWedgeBasis() == PosInRes)
    {
      std::map<std::pair<long, MorseElement>, long>::const_iterator RowIter(identifier.find(make_pair(degree, ResIter->first)));
      if (RowIter != identifier.end())
      {
        const long PositionRow(RowIter->second);
        PathMap maps((ResIter->second).myGetPaths());
        for (PathMap::const_iterator MapIter = maps.begin(); MapIter != maps.end(); ++ MapIter)
        {
          std::map<std::pair<long, MorseElement>, long>::const_iterator ColIter(identifier.find(make_pair(degree, (MapIter->first)->first)));
          if (ColIter != identifier.end())
          {
            const long PositionCol(ColIter->second);
            SetEntry(SubMat, PositionRow, PositionCol, MapIter->second);
          }
        }
      }
      ++ResIter;
    }
    return SubMat;      
  }


  std::vector< std::vector<long> > MorseGraph::myComputeRanks() const
  {
    const long LengthRows = ((myResolution.rbegin())->first).myCountWedgeBasis() + 1;
    long LengthCols(0);
    std::vector< std::pair<long, long> > RanksVector;
    for (std::map<MorseElement, MorsePaths>::const_iterator it = myResolution.begin(); it != myResolution.end(); ++it)
    {
      RanksVector.push_back(make_pair((it->first).myCountWedgeBasis(), (it->first).myDegree()));
      if ((it->first).myDegree() > LengthCols)
      {
        LengthCols = (it->first).myDegree();
      }
    }
    std::vector< std::vector<long> > ranks (LengthRows, std::vector<long>(LengthCols + 1, 0));
    for (std::vector<std::pair<long, long> >::iterator it = RanksVector.begin(); it != RanksVector.end(); ++it)
    {
      ++ranks[it->first][it->second];
    }
    return ranks;
  }


  matrix MorseGraph::myTransformRanksToBettis(const std::vector<std::vector<long> >& ranks) const
  {
    const long RowsBetti(myNumRowsBettiDiagram(ranks));
    const long LenRanks = len(ranks);
    if (LenRanks == 1 && len(ranks[0]) == 1)
    {
      return NewDenseMat(RingZZ(), 1, 1);
    }
    matrix bettis(NewDenseMat(RingZZ(), RowsBetti, LenRanks + 1));
    for (long i = 0; i < LenRanks; ++i)
    {
      for (long j = i + 1; j < RowsBetti + i + 1; ++j)
      {
        SetEntry(bettis, j - i - 1, i + 1, ranks[i][j]);
      }
    }
    SetEntry(bettis, 0, 0, 1);
    for (long i = 1; i < NumRows(bettis); ++i)
    {
      SetEntry(bettis, i, 0, 0);
    }
    return bettis;
  }


  long MorseGraph::myNumRowsBettiDiagram(const std::vector<std::vector<long> >& ranks) const
  {
    CoCoA_ASSERT(!ranks.empty());
    const long res(len(ranks[0]) - len(ranks)); //skipping 0 degree
    if (res <= 0)
    {
      return 1;
    }
    return res;
  }


  std::vector<matrix> MorseGraph::myMapsAsMatrices() const
  {
    std::vector<matrix> res;
    std::map<MorseElement, MorsePaths>::const_iterator iter(myResolution.end());
    --iter;
    long size((iter->first).myCountWedgeBasis());
    for (long i = 0; i < size; ++i)
    {
      res.push_back(myMapsAsMatrix(i));
    }
    return res;
  }


  matrix MorseGraph::myMapsAsMatrix(long pos) const
  {
    matrix res = myInitialMapsAsMatrix(pos);
    std::map<MorseElement, long> identifier(myPositionMorseElements(pos));
    for (std::map<MorseElement, MorsePaths>::const_iterator it = myResolution.begin(); it != myResolution.end(); ++it)
    {
      if ((it->first).myCountWedgeBasis() == pos)
      {
        const long PositionRow(identifier[it->first]);
        PathMap maps((it->second).myGetPaths());
        for (PathMap::iterator MapIter = maps.begin(); MapIter != maps.end(); ++ MapIter)
        {
          const long PositionCol(identifier[(MapIter->first)->first]);
          SetEntry(res, PositionRow, PositionCol, MapIter->second);
        }
      }
    }
    return res;
  }


  matrix MorseGraph::myInitialMapsAsMatrix(long pos) const
  {
    long row(0);
    long col(0);
    for (std::map<MorseElement, MorsePaths>::const_iterator it = myResolution.begin(); it != myResolution.end(); ++it)
    {
      const long NumWedgeBasis((it->first).myCountWedgeBasis());
      if (pos > NumWedgeBasis)
      {
        continue;
      }
      if (pos == NumWedgeBasis)
      {
        ++row;
      }
      if ((pos + 1) == NumWedgeBasis)
      {
        ++col;
      }
      if ((pos + 1) < NumWedgeBasis)
      {
        break;
      }
    }
    return NewDenseMat(myMapRing, row, col);
  }


  std::vector<matrix> MorseGraph::myInitialMapsAsMatrices() const
  {
    std::vector<matrix> res;
    long SizePrecedingWedge(0);
    long row(-1);
    long col(0);
    for (std::map<MorseElement, MorsePaths>::const_iterator it = myResolution.begin(); it != myResolution.end(); ++it)
    {
      if (SizePrecedingWedge != (it->first).myCountWedgeBasis())
      {
        SizePrecedingWedge = (it->first).myCountWedgeBasis();
        if (row > 0)
        {
          res.push_back(NewDenseMat(myMapRing, row, col));
        }
        row = col;
        col = 0;
      }
      ++col;
    }
    if (row > 0)
    {
      res.push_back(NewDenseMat(myMapRing, row, col));
    }
    return res;
  }


// we can now access an element via <NumberWedges, Position In NumberWedges>
  std::map<MorseElement, long> MorseGraph::myPositionMorseElements(long NumberWedges) const
  {
    std::map<MorseElement, long> res;
    long pos1(0);
    long pos2(0);
    for (std::map<MorseElement, MorsePaths>::const_iterator it = myResolution.begin(); it != myResolution.end(); ++it)
    {
      if ((it->first).myCountWedgeBasis() == NumberWedges)
      {
        res.insert(make_pair(it->first, pos1));
        ++pos1;
      }
      if ((it->first).myCountWedgeBasis() == NumberWedges + 1)
      {
        res.insert(make_pair(it->first, pos2));
        ++pos2;
      }
    }
    return res;
  }


// we can now access an element via <NumberWedges, Position In NumberWedges>
// orderd by degree
  std::map<std::pair<long, MorseElement>, long> MorseGraph::myGradedPositionMorseElements(long NumberWedges) const
  {
    std::map<std::pair<long, MorseElement>, long> res;
    long pos(0);
    long degree(0);
    long OldNumberWedges(0);
    for (std::map<MorseElement, MorsePaths>::const_iterator it = myResolution.begin(); it != myResolution.end(); ++it)
    {
      const long CurrentNumberWedges((it->first).myCountWedgeBasis());
      if (CurrentNumberWedges < NumberWedges)
      {
        continue;
      }
      if (degree != (it->first).myDegree() || CurrentNumberWedges != OldNumberWedges)
      {
        pos = 0;
        degree = (it->first).myDegree();
        if (CurrentNumberWedges > NumberWedges + 1)
        {
          break;
        }
        OldNumberWedges = CurrentNumberWedges;
      }
      res.insert(make_pair(make_pair(degree, it->first), pos));
      ++pos;
    }
    return res;
  }


  std::vector<std::vector<long> > MorseGraph::myMatrixMinus(std::vector<std::vector<long> > m1, const std::vector<std::vector<long> >& m2) const
  {
    const long nrows = len(m1);
    const long ncols = len(m1[0]);
    for (long i = 0; i < nrows; ++i)
    {
      for (long j = 0; j < ncols; ++j)
      {
        m1[i][j] -= m2[i][j];
      }
    }
    return m1;
  }


  std::vector<matrix> MorseGraph::myComputeMinimalResolution()
  {
    //initialization
    myResolution.clear();
    StandardRepresentationContainer container(myMill);
    //compute resolution
    myComputeBasicGraph(myComputeGeneralBasis(), container);
    myDirectMorseReduction(container);
    //minimalization
    std::vector<matrix> maps(1, myZerothMatrix());
    std::vector<matrix> OtherMaps(myMapsAsMatrices());
    maps.insert(maps.end(), OtherMaps.begin(), OtherMaps.end());
    ResolutionMinimization minim(myRing, maps);
    minim.myMinimization();
    return minim.myGetResolution();
  }


  matrix MorseGraph::myZerothMatrix() const
  {
    std::vector<RingElem> vec;
    for (std::map<MorseElement, MorsePaths>::const_iterator it = myResolution.begin(); it != myResolution.end(); ++it)
    {
      if (0 != (it->first).myCountWedgeBasis())
      {
        break;
      }
      vec.push_back((it->first).myGetBasisElement());
    }
    CoCoA_ASSERT(!vec.empty());
    return NewDenseMat(RowMat(vec));
  }


  std::vector<matrix> MorseGraph::myGetResolution() const
  {
    std::vector<matrix> res(1, myZerothMatrix());
    std::vector<matrix> OtherMaps(myMapsAsMatrices());
    res.insert(res.end(), OtherMaps.begin(), OtherMaps.end());
    return res;
  }


  std::vector<std::pair<RingElem, DynamicBitset> > MorseGraph::myTransform(const std::vector<std::pair<RingElem, std::vector<bool> > >& VecWithBool) const
  {
    std::vector<std::pair<RingElem, DynamicBitset> > res; 
    for (std::vector<std::pair<RingElem, std::vector<bool> > >::const_iterator it = VecWithBool.begin(); it != VecWithBool.end(); ++it)
    {
      res.push_back(make_pair(it->first, myVectorBoolToDynamicBitset(it->second)));
    }
    return res;
  }


  void MorseGraph::myComputeConstantResolution()
  {
    myResolution.clear();
    myMapRing = CoeffRing(myRing);
    StandardRepresentationContainer container(myMill);
    myComputeBasicConstantGraph(myComputeGeneralBasis(), container);
  
    myConstantDirectMorseReduction(container);
  }


  matrix JBBettiDiagram(JBMill mill)
  {
    if (!(IsStdDegRevLex(ordering(owner(LPP(*(JBReturnJB(mill).begin())))))))
    {
      CoCoA_ERROR(ERR::PPOrder, "It must be the degrevlex ordering!!!");
    }
    if (!JBIsPommaretBasis(mill) || !JBIsHomogenous(mill))
    {
      CoCoA_ERROR("The ideal isn't homogenous or pommaret", "JBBettiDiagram");
    }
    MorseGraph mg(mill);
    return mg.myComputeBettiNumbers();
  }

  /*
   * computes the minimal free Resolution
   */
  std::vector<matrix> JBMinimalResolution(JBMill mill)
  {
    if (!(IsStdDegRevLex(ordering(owner(LPP(JBReturnJB(mill).front()))))))
    {
      CoCoA_ERROR(ERR::PPOrder, "It must be the degrevlex ordering!!!");
    }
    if (!JBIsPommaretBasis(mill) || !JBIsHomogenous(mill))
    {
      CoCoA_ERROR("The ideal isn't homogenous or pommaret", "JBMinimalResolution");
    }
    MorseGraph mg(mill);
    return mg.myComputeMinimalResolution();

  }


  /*
   * computes a free Resolution
   */
  std::vector<matrix> JBResolution(JBMill mill)
  {
    if (!(IsStdDegRevLex(ordering(owner(LPP(*(JBReturnJB(mill).begin())))))))
    {
      CoCoA_ERROR(ERR::PPOrder, "It must be the degrevlex ordering!!!");
    }
    if (!JBIsPommaretBasis(mill))
    {
      CoCoA_ERROR("The ideal isn't pommaret", "JBResolution");
    }
    MorseGraph mg(mill);
    mg.myComputeResolution();
    std::vector<matrix> maps(1, mg.myZerothMatrix());
    std::vector<matrix> OtherMaps(mg.myMapsAsMatrices());
    maps.insert(maps.end(), OtherMaps.begin(), OtherMaps.end());
    return maps;
  }


} // end of namespace CoCoA
