//   Copyright (c)  2006  John Abbott

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


#include "CoCoA/QBGenerator.H"
#include "CoCoA/assert.H"
#include "CoCoA/error.H"
#include "CoCoA/VectorOperations.H"   // for template fn to output vector/list.

//#include <list>
using std::list;
#include <iostream>
using std::ostream;
#include <algorithm>
using std::find_if;

namespace CoCoA
{

  namespace // anonymous namespace for file local definitions
  {

    class IsFactorOf
    {
    public:
      explicit IsFactorOf(ConstRefPPMonoidElem pp);
      bool operator()(ConstRefPPMonoidElem candidate);
    private:
      PPMonoidElem myPP;
    };

    inline IsFactorOf::IsFactorOf(ConstRefPPMonoidElem pp):
        myPP(pp)
    {}

    inline bool IsFactorOf::operator()(ConstRefPPMonoidElem candidate)
    {
      return IsDivisible(myPP, candidate);
    }

  } // end of anonymous namespace


  QBGenerator::QBGenerator(const PPMonoid& PPM):
      myPPMValue(PPM),
      myCornerList(),
      myNewCornerList(),
      myAvoidList(),
      myQBList()
  {
    myCornerList.push_back(one(PPM));
    myNewCornerList.push_back(one(PPM));
  }


  void QBGenerator::myCornerPPIntoQB(PPMonoidElem pp) // NB I do want a private copy of the arg!!!
  {
    CoCoA_ASSERT(owner(pp) == myPPMValue);
    list<PPMonoidElem> CornerList(myCornerList); // local copy to allow exception safety
    list<PPMonoidElem>::iterator it = find(CornerList.begin(), CornerList.end(), pp);
    if (it == CornerList.end())
      CoCoA_ERROR("PP does not belong to the corner set", "myCornerPPIntoQB");
    CornerList.erase(it); // remove pp from CornerList

    // To ensure EXCEPTION SAFETY, we wait until later to put pp into myQBList!

    // Now find out which multiples of pp (if any) are to be put into the corner lists.
    // We consider pp*x[i] in turn for each indeterminate x[i].
    list<PPMonoidElem> NewCornerList;
    PPMonoidElem ppx(myPPMValue);
    for (long xi=0; xi < NumIndets(myPPMValue); ++xi)
    {
      ppx = pp * indet(myPPMValue, xi);
      if (find_if(myAvoidList.begin(), myAvoidList.end(), IsFactorOf(ppx)) != myAvoidList.end()) continue;
      if (find_if(CornerList.begin(), CornerList.end(), IsFactorOf(ppx)) != CornerList.end()) continue;
      NewCornerList.push_back(ppx);
    }

    // Sort the new corner elements, and merge a copy of them into the full set of corner elements.
    NewCornerList.sort();
    list<PPMonoidElem> CopyOfNewCornerElems(NewCornerList);
    myQBList.push_back(pp);
    CornerList.merge(CopyOfNewCornerElems); // merge moves the PPs into CornerList
    myNewCornerList.swap(NewCornerList);
    myCornerList.swap(CornerList);
  }


  void QBGenerator::myCornerPPIntoAvoidSet(ConstRefPPMonoidElem pp)
  {
    CoCoA_ASSERT(owner(pp) == myPPMValue);
    const list<PPMonoidElem>::iterator it = find(myCornerList.begin(), myCornerList.end(), pp);
    if (it == myCornerList.end())
      CoCoA_ERROR("PP does not belong to the corner set", "myCornerPPIntoAvoidSet");
    myAvoidList.push_back(pp); // Do this *before* deleting pp from myCornerList as pp may alias the element we delete!
    myCornerList.erase(it);
    myNewCornerList.clear();
  }


  const std::list<PPMonoidElem>& QBGenerator::myNewCorners() const
  {
    return myNewCornerList;
  }


  const std::list<PPMonoidElem>& QBGenerator::myCorners() const
  {
    return myCornerList;
  }


  const std::vector<PPMonoidElem>& QBGenerator::myQB() const
  {
    return myQBList;
  }


  const PPMonoid& QBGenerator::myPPM() const
  {
    return myPPMValue;
  }


  void QBGenerator::myOutputSelf(std::ostream& out) const
  {
    out << "QBGenerator("
        << "QB=" << myQBList
        << ", corners=" << myCornerList
        << ", NewCorners=" << myNewCornerList
        << ", avoid=" << myAvoidList
        << ")";
  }


  const PPMonoid& PPM(const QBGenerator& QBG)
  {
    return QBG.myPPM();
  }


  std::ostream& operator<<(std::ostream& out, const QBGenerator& QBG)
  {
    QBG.myOutputSelf(out);
    return out;
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/QBGenerator.C,v 1.5 2014/07/31 14:45:18 abbott Exp $
// $Log: QBGenerator.C,v $
// Revision 1.5  2014/07/31 14:45:18  abbott
// Summary: Merged io.H and UtilsTemplate.H into new header VectorOperations.H
// Author: JAA
//
// Revision 1.4  2013/01/21 13:29:53  abbott
// Changed impl so that it is exception safe.
//
// Revision 1.3  2011/03/10 16:39:34  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.2  2007/10/30 17:14:07  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.6  2007/02/10 18:44:03  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.5  2007/01/09 15:52:08  cocoa
// Changed QBGenerator to use std::vector instead of std::list for the result.
// Minor mod to configure script.
//
// Revision 1.4  2006/11/24 17:06:10  cocoa
// -- reorganized includes of header files
//
// Revision 1.3  2006/10/06 14:04:14  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.2  2006/08/07 21:23:25  cocoa
// Removed almost all publicly visible references to SmallExponent_t;
// changed to long in all PPMonoid functions and SparsePolyRing functions.
// DivMask remains to sorted out.
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.2  2006/04/27 13:45:30  cocoa
// Changed name of NewIdentityRingHom to NewIdentityHom.
// Changed name of member functions which print out their own object
// into myOutputSelf (to distinguish from "transitive" myOutput fns).
//
// Revision 1.1  2006/04/21 15:03:23  cocoa
// New code for Buchberger-Moeller and variants.
//
//
