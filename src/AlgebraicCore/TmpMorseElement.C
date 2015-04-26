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


#include "CoCoA/TmpMorseElement.H"

using std::make_pair;

namespace CoCoA
{

  long myMaxTrueOfBitset(DynamicBitset bitset)
  {
    const long length(bitset.myLen());
    for (long i = length - 1; i > -1; --i)
    {
      if (bitset.Iam1At(i))
      {
        return i;
      }
    }
    return -1;
  }


  DynamicBitset myVectorLongToDynamicBitset(const std::vector<long>& longs, const long& length)
  {
    DynamicBitset result(length);
    for (std::vector<long>::const_iterator it = longs.begin(); it != longs.end(); ++it)
    {
      result.mySet(*it);
    }
    return result;
  }

  DynamicBitset myVectorBoolToDynamicBitset(const std::vector<bool>& BoolVector)
  {
    const long n = len(BoolVector);
    DynamicBitset result(n);
    for (long i=0; i < n; ++i)
      result.mySet(i, BoolVector[i]);
    return result;
  }

  std::vector<long> myDynamicBitsetToLong(const DynamicBitset& DynBitset)
  {
    std::vector<long> result;
    long length(len(DynBitset));
    for (long i = 0; i < length; ++i)
    {
      if (DynBitset.Iam1At(i))
      {
        result.push_back(i);
      }
    }
    return result;
  }

  const std::vector<RingElem>& StandardRepresentationContainer::myComputeStandardRepresentation(ConstRefRingElem r)
  {
    typedef std::pair<RingElem, std::vector<RingElem> > StdRepr;
    ++myOriginNormalForms;
    std::multimap<PPMonoidElem, StdRepr>::iterator res(myContainer.end());
    std::pair<std::multimap<PPMonoidElem, StdRepr>::iterator,
      std::multimap<PPMonoidElem, StdRepr>::iterator> iters(myContainer.equal_range(LPP(r)));
  for (std::multimap<PPMonoidElem, StdRepr>::iterator it = iters.first; it != iters.second; ++it)
  {
    if ((it->second).first == r)
    {
      res = it;
      break;
    } 
  }
  if (res == myContainer.end())
  {
    ++myReallyNormalForms;
    res = myContainer.insert(make_pair(LPP(r), make_pair(r, myMill.myStandardRepresentationWithoutRestShort(r))));
  }
  return (res->second).second;
}


  //------------------------------------------------------------------
  // MorseElement

  // 2 ctors
  MorseElement::MorseElement(const DynamicBitset& WedgeProduct, const PommBasisElem basis)
      : myWedgeProduct(WedgeProduct)
      , myRightFactor(LPP(one(owner(basis->first))))
      , myBasis(basis)
      , myRightProduct(LPP(myBasis->first))
      , myProduct(myRightProduct * NewPP(owner(myRightFactor), myWedgeProduct))
  {
    myWedgeProductAsLongs = myDynamicBitsetToLong(myWedgeProduct);
  }

  MorseElement::MorseElement(const DynamicBitset& WedgeProduct, const PPMonoidElem& myRightFactor, const PommBasisElem basis)
      : myWedgeProduct(WedgeProduct)
      , myRightFactor(myRightFactor)
      , myBasis(basis)
      , myRightProduct(LPP(myBasis->first) * myRightFactor)
      , myProduct(myRightProduct * NewPP(owner(myRightFactor), myWedgeProduct))
  {
    myWedgeProductAsLongs = myDynamicBitsetToLong(myWedgeProduct);
  }


  void MorseElement::mySetWedgeProduct(const DynamicBitset& elem)
  {
    myProduct = (myProduct / NewPP(owner(myRightFactor), myWedgeProduct)) * NewPP(owner(myRightFactor), elem);
    myWedgeProduct = elem;
    myWedgeProductAsLongs = myDynamicBitsetToLong(myWedgeProduct);
  }


  void MorseElement::mySetRightFactor(const PPMonoidElem& elem)
  {
    myProduct = (myProduct / myRightFactor) * elem;
    myRightFactor = (myRightFactor / myRightFactor) * elem;
    myRightFactor = elem;
  }


  void MorseElement::myDivideRightProductWith(long i)
  {
    const PPMonoidElem t = indet(owner(myProduct), i);
    myProduct = myProduct / t;
    myRightProduct = myRightProduct / t;
    myRightFactor = myRightFactor / t;
  }


  bool operator <(const MorseElement& m1, const MorseElement& m2)
  {
    const long l1 = len(m1.myWedgeProductAsLongs);
    const long l2 = len(m2.myWedgeProductAsLongs);
    if (l1 != l2) return (l1 < l2);
    // Same length...
    const int CmpProduct(owner(m1.myProduct)->myCmp(raw(m1.myProduct), raw(m2.myProduct)));
    if (CmpProduct != 0) return (CmpProduct < 0);

    // products are equal...
    if (m1.myBasis->second != m2.myBasis->second)
      return IsSubset(m2.myBasis->second, m1.myBasis->second);
    // Same mult vars...
    for (long i=0; i < l1; ++i)
      if (m1.myWedgeProductAsLongs[i] != m2.myWedgeProductAsLongs[i])
        return (m2.myWedgeProductAsLongs[i] < m1.myWedgeProductAsLongs[i]);

//    return (m1.myRightProduct < m2.myRightProduct);
    return owner(m1.myProduct)->myCmp(raw(m1.myRightProduct), raw(m2.myRightProduct)) < 0;
  }

  bool operator <=(const MorseElement& m1, const MorseElement& m2) { return !(m2 < m1); }
  bool operator > (const MorseElement& m1, const MorseElement& m2) { return  (m2 < m1); }
  bool operator >=(const MorseElement& m1, const MorseElement& m2) { return !(m1 < m2); }


  int MorseElement::myEpsilon(long test, long add) const
  {
    bool even = true;
    if ((add < test) && !myWedgeProduct.Iam1At(test))
      even = !even;
    for (long i = 0; i < test; ++i)
    {
      even ^= myWedgeProduct.Iam1At(i);
    }
    if (even) return 1;
    else return -1;
  }


  long MorseElement::myMaxTypeOne() const
  {
    const DynamicBitset support(myRightFactor);
    return myMaxTrueOfBitset((support) - myWedgeProduct);
  }


  long MorseElement::myMaxTypeTwo() const
  {
    const DynamicBitset support(myRightFactor);
    return myMaxTrueOfBitset((support | myWedgeProduct) & myBasis->second);
  }


  std::vector<std::pair<MorseElement, RingElem> > MorseElement::myComputeBasicMaps(const std::pair<PommBasisElem, PommBasisElem>& BasisIters, StandardRepresentationContainer& container) const
  {
    std::vector< std::pair<MorseElement, RingElem> > result;
    const long length(len(myWedgeProduct));
    for (std::vector<long>::const_iterator it = myWedgeProductAsLongs.begin(); it != myWedgeProductAsLongs.end(); ++it)
    {
      DynamicBitset NewWedge(myWedgeWithOneRemoved(myWedgeProductAsLongs, length, it));
      myComputeLeftMaps(result, *it, NewWedge);
      myComputeRightMaps(result, *it, NewWedge, BasisIters, container, myGetPolyRing());
    }
    return result; 
  }


  std::vector<std::pair<MorseElement, RingElem> > MorseElement::myComputeBasicConstantMaps(const std::pair<PommBasisElem, PommBasisElem>& BasisIters, StandardRepresentationContainer& container) const
  {
    std::vector<std::pair<MorseElement, RingElem> > result;
    const long length(len(myWedgeProduct));
    for (std::vector<long>::const_iterator it = myWedgeProductAsLongs.begin(); it != myWedgeProductAsLongs.end(); ++it)
    {
      DynamicBitset NewWedge(myWedgeWithOneRemoved(myWedgeProductAsLongs, length, it));
      myComputeRightMaps(result, *it, NewWedge, BasisIters, container, CoeffRing(myGetPolyRing()));
    }
    return result; 
  }


  DynamicBitset MorseElement::myWedgeWithOneRemoved(const std::vector<long>& WedgeAsLong, long LengthWedge, std::vector<long>::const_iterator it) const
  {
    DynamicBitset ans = myVectorLongToDynamicBitset(WedgeAsLong, LengthWedge);
    ans.mySet(*it, false);
    return ans;
  }


  void MorseElement::myComputeLeftMaps(std::vector<std::pair<MorseElement, RingElem> >& maps, long i, DynamicBitset NewWedge) const
  {
    MorseElement m1(NewWedge, myBasis);
    const SparsePolyRing P = owner(myBasis->first);
    maps.push_back(make_pair(m1, myEpsilon(i,i)*indet(P,i)));
  }


  void MorseElement::myComputeRightMaps(std::vector<std::pair<MorseElement, RingElem> >& maps,
                                        long i,
                                        DynamicBitset NewWedge,
                                        const std::pair<PommBasisElem, PommBasisElem>& BasisIters,
                                        StandardRepresentationContainer& container,
                                        const ring& MapRing)  const
  {
    const SparsePolyRing P = owner(myBasis->first);

    const std::vector<RingElem>& SecondPart(container.myComputeStandardRepresentation(myBasis->first * indet(P,i)));
    PommBasisElem BasisIter(BasisIters.first);
    for (std::vector<RingElem>::const_iterator RepIter = SecondPart.begin(); RepIter != SecondPart.end(); ++RepIter)
    {
      for (SparsePolyIter SPI = BeginIter(*RepIter); !IsEnded(SPI); ++SPI)
      {
        MorseElement m(NewWedge, PP(SPI), BasisIter);
        RingElem map(MapRing, (-myEpsilon(i, i)) * coeff(SPI));
        maps.push_back(make_pair(m, map));
      }
      ++BasisIter;
    }
    CoCoA_ASSERT(BasisIter == BasisIters.second);
  }

} // end of namespace CoCoA
