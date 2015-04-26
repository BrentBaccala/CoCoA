//   Copyright (c)  2006-2010  Anna Bigatti, Massimo Caboara

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


#include "CoCoA/DynamicBitset.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/VectorOperations.H" // for printing vectors
//#include "CoCoA/PPMonoid.H"

#include <algorithm>
using std::transform;
//#include <bitset>
using std::bitset;
#include <functional>
//??
#include <iterator>
using std::back_inserter;  // for transform
#include <iostream>
//using std::ostream;
///??? UNUSED!!! //#include <list>
///??? UNUSED!!! using std::list;
//#include <vector>
using std::vector;

//static const bool MAX_DEBUG = false;

namespace CoCoA
{
  DynamicBitset::OutputStyle DynamicBitset::ourOutputStyle; // clean;

// class DynamicBitset


  void DynamicBitset::myResize(long n) // only for ctors
  {
    CoCoA_ASSERT(n >= 0);
    myLenValue = n;
    if (n == 0) return;
    myVec.resize((n-1)/ourNumBitsInBlock + 1); // all 0s
  }
  

  DynamicBitset::DynamicBitset(ConstRefPPMonoidElem pp)
  {
    vector<long> v(NumIndets(owner(pp)));
    exponents(v, pp);
    // same as DynamicBitset(v)
    const long n = len(v);
    myResize(n);
    for (long i=0; i!=n; ++i)
      if (v[i] != 0)  mySet(i);
  }


//   DynamicBitset::DynamicBitset(const vector<long>& v)
//   {
//     myResize(len(v));
//     for (long i=0; i!=len(v); ++i)
//       if (v[i] != 0)  mySet(i);
//   }


  DynamicBitset::DynamicBitset(long n)
  {
    if (n < 0) CoCoA_ERROR(ERR::NotNonNegative, "DynamicBitset(long n)");
    myResize(n);
  }


  DynamicBitset::DynamicBitset(const DynamicBitset& rhs)
  {
    myLenValue = rhs.myLenValue;
    myVec = rhs.myVec;
  }


  // Anna 15Apr2010:
  // if ">>=" appears to be too slow make vector of precomputed masks (static)
  bool DynamicBitset::IamAll1s() const
  {
    const BitBlock mask = ~0; // 111111111
    for (vector<BitBlock>::const_iterator it=myVec.begin(); it!=myVec.end(); ++it)
      if ((std::operator^(*it, mask)).any())  // xor for bitset
      {
        if (it+1 != myVec.end() || myLenValue%ourNumBitsInBlock==0) return false;
        const int shift = (ourNumBitsInBlock - myLenValue%ourNumBitsInBlock);
//JAA        mask >>= (ourNumBitsInBlock - myLenValue%ourNumBitsInBlock); // 00000111111
        if (std::operator^(*it, mask>>shift).any()) return false;
      }
    return true;
  }
  

  void DynamicBitset::myOutputSelf(std::ostream& out) const
  {
    if (myLenValue==0)  return;
    vector<BitBlock>::const_reverse_iterator rit=myVec.rbegin();
    // first block without spurious 0s
    for (long n=(myLenValue-1)%ourNumBitsInBlock; n>=0; --n)
      out << rit->test(n);
    //    for (++rit; rit!=myVec.rend(); ++rit)  out << '-' << *rit;
    for (++rit; rit!=myVec.rend(); ++rit) out << *rit;
  }


  void DynamicBitset::myOutputSelf8(std::ostream& out) const
  {
    if (myLenValue==0)  return;
    vector<BitBlock>::const_reverse_iterator rit=myVec.rbegin();
    // first block without spurious 0s
    for (long n=(myLenValue-1)%ourNumBitsInBlock; n>=0; --n)
    {
      out << rit->test(n);
      if (n%8==0 && n!=0)  out << '.';
    }
    for (++rit; rit!=myVec.rend(); ++rit)
    {
      out << '-';
      for (long n=ourNumBitsInBlock-1; n>=0; --n)
      {
        out << rit->test(n);
        if (n%8==0 && n!=0)  out << '.';
      }
    }
  }


  namespace // Anna 15Apr2010: maybe not the best way.  Just practising STL.
  {
    inline unsigned long ToULong(const DynamicBitset::BitBlock& b)
    { return b.to_ulong(); }
  }
  
  void DynamicBitset::myOutputSelfLong(std::ostream& out) const
  {
    vector<unsigned long> v;
    v.reserve(len(myVec));
    transform(myVec.rbegin(), myVec.rend(), back_inserter(v), ToULong);
    out << v;
  }


  long count(const DynamicBitset& b)
  {
    long NumOnes = 0;
    for (vector<DynamicBitset::BitBlock>::const_iterator it=b.myVec.begin(); it!=b.myVec.end(); ++it)
      NumOnes += it->count();
    return NumOnes;
  }


  DynamicBitset flip(DynamicBitset DB)
  {
    const DynamicBitset::BitBlock All1s = ~0;
    const long w = DynamicBitset::ourNumBitsInBlock;
    const long lenDB = DB.myLenValue;
    const long n = lenDB/w; // integer division!
    for (long i=0; i < n; ++i)
      DB.myVec[i] ^= All1s;

    const int shift = (w - lenDB%w);
    if (shift != 0)
      DB.myVec.back() ^= (All1s >> shift);
    return DB;
  }


  std::ostream& operator<<(std::ostream& out, const DynamicBitset& DB)
  {
    switch (DynamicBitset::ourOutputStyle)
    {
    case DynamicBitset::clean:          DB.myOutputSelf(out);  break;
    case DynamicBitset::AsRevVecOfLong: DB.myOutputSelfLong(out);  break;
    case DynamicBitset::WithSeparators: DB.myOutputSelf8(out);  break;
    default: CoCoA_ERROR(ERR::SERIOUS, "operator<<(ostream, DynamicBitset)");
    }
    return out;
  }


  // can be improved making vector of exponents (for dense PPMonoid)
  PPMonoidElem NewPP(const PPMonoid& PPM, const DynamicBitset& b)
  {
    if (NumIndets(PPM) != len(b))  CoCoA_ERROR(ERR::MixedSizes, "NewPP");
    PPMonoidElem pp(PPM);
    for (long i=0; i!=len(b); ++i)
      if (b.Iam1At(i))  pp *= indet(PPM, i);
    return pp;
  }


}// end namespace cocoa





/*

Some future optimization:

for IsFLesser: proceed in the test word by word:
compute the first word of g1-f, the first word of g2-f, the check

ConnectionBlock: use ptr and not iterators.

sparse representation for facets? Use vector, (sort?).
Only reasonable if density is much much lower than #VARS!

*/

// RCS header/log on the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/DynamicBitset.C,v 1.16 2014/07/31 14:45:17 abbott Exp $
// $Log: DynamicBitset.C,v $
// Revision 1.16  2014/07/31 14:45:17  abbott
// Summary: Merged io.H and UtilsTemplate.H into new header VectorOperations.H
// Author: JAA
//
// Revision 1.15  2014/04/30 16:06:14  abbott
// Summary: Replaced X.size() by len(X)
// Author: JAA
//
// Revision 1.14  2013/06/27 16:47:48  abbott
// Added flip fn.
//
// Revision 1.13  2013/04/16 14:45:25  abbott
// Added new fn count for DynamicBitset.
//
// Revision 1.12  2012/09/28 14:05:52  abbott
// Commented out apparently useless include of STL list stuff.
//
// Revision 1.11  2012/01/26 16:53:04  bigatti
// -- added #include <iterator>
//
// Revision 1.10  2011/03/11 16:50:31  bigatti
// -- changed  unsigned int --> long
//
// Revision 1.9  2011/03/11 12:36:21  bigatti
// -- changed size_t --> long
// -- changed size --> len
//
// Revision 1.8  2010/05/28 15:50:42  bigatti
// -- cleaning
// -- commented out ctor taking a list of long
// -- moved some "facet" functions to TmpIsTree
//
// Revision 1.7  2010/04/27 16:08:58  bigatti
// -- changed error messages
//
// Revision 1.6  2010/04/21 11:35:55  bigatti
// -- changed Iam1InPos --> Iam1At
//
// Revision 1.5  2010/04/21 09:43:49  bigatti
// -- added static member field for choosing printing style
//
// Revision 1.4  2010/04/17 16:47:29  bigatti
// -- added ctor with ConstRefPPMonoidElem
// -- added include <algorithm>
// -- more consistent myOutputSelf* functions
//
// Revision 1.3  2010/04/15 16:01:39  bigatti
// -- DynamicBitset going towards final version
//
// Revision 1.2  2010/04/13 15:30:25  bigatti
// -- reorganized into (almost) final shape
//
// Revision 1.1  2010/03/30 15:20:43  bigatti
// -- was "facet" in TmpIsTree
//
