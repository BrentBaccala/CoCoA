//   Copyright (c)  2005  John Abbott

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


#include "CoCoA/RingHom.H"
#include "CoCoA/error.H"
#include "CoCoA/ideal.H"

#include <iostream>
using std::ostream;
#include<vector>
using std::vector; // used only in IsInKer

namespace CoCoA
{

  //-------------------------------------------------------
  // CLASS PartialRingHom

  RingElem PartialRingHom::operator()(ConstRefRingElem x) const
  {
    if (owner(x) != domain(*this))
      CoCoA_ERROR(ERR::BadRingHomArg, "Applying PartialRingHom to RingElem");
    RingElem ans(codomain(*this));
    mySmartPtr->myApply(raw(ans), raw(x));
    return ans;
  }


  RingElem PartialRingHom::operator()(const MachineInt& n) const
  {
    RingElem arg(domain(*this), n);
    RingElem ans(codomain(*this));
    mySmartPtr->myApply(raw(ans), raw(arg));
    return ans;
  }

  RingElem PartialRingHom::operator()(const BigInt& N) const
  {
    RingElem arg(domain(*this), N);
    RingElem ans(codomain(*this));
    mySmartPtr->myApply(raw(ans), raw(arg));
    return ans;
  }

  RingElem PartialRingHom::operator()(const BigRat& q) const
  {
    RingElem arg(domain(*this), q);
    RingElem ans(codomain(*this));
    mySmartPtr->myApply(raw(ans), raw(arg));
    return ans;
  }


  PartialRingHom PartialRingHom::operator()(const RingHom& theta) const
  {
    if (codomain(theta) != domain(*this))
      CoCoA_ERROR(ERR::BadCompose, "composition PartialRingHom(RingHom)");
// BUG BUG BUG can be clever here!!!
    return sequential(*this, theta);
  }


  PartialRingHom PartialRingHom::operator()(const PartialRingHom& theta) const
  {
    if (codomain(theta) != domain(*this))
      CoCoA_ERROR(ERR::BadCompose, "composition PartialRingHom(PartialRingHom)");
    return sequential(*this, theta);
  }


  void RingHomBase::myOutputSelf(std::ostream& out) const
  {
    out << "PartialRingHom(" << myDomain << " --> " << myCodomain;
    myOutputSelfDetails(out);
    out << ")";
  }


  void RingHomBase::myOutputSelfDetails(std::ostream& /*out*/) const
  {
    // Default definition does nothing (as there are no extra details to print).
    // SHOULD THIS BE PURE VIRTUAL (with no default defn)???
  }


  std::ostream& operator<<(std::ostream& out, const PartialRingHom& phi)
  {
    out << "RingHom(" << domain(phi) << " --> " << codomain(phi);
    phi->myOutputSelfDetails(out);
    out << ")";
    return out;
  }


  // better to use indices or iterators in this fn???
// USE AN ALGORITHM!!!
  bool IsInKer(const ideal& I, const PartialRingHom& phi)
  {
    const vector<RingElem>& g = gens(I);
    const vector<RingElem>::const_iterator stop = g.end();
    for (vector<RingElem>::const_iterator it = g.begin(); it != stop; ++it)
    {
      if (!IsZero(phi(*it))) return false;
    }
    return true;
  }


  //-------------------------------------------------------
  // CLASS RingHom

  RingElem RingHom::operator()(ConstRefRingElem x) const
  {
    if (owner(x) != domain(*this))
      CoCoA_ERROR(ERR::BadRingHomArg, "Applying RingHom to RingElem");
    RingElem ans(codomain(*this));
    mySmartPtr->myApply(raw(ans), raw(x));
    return ans;
  }


  RingElem RingHom::operator()(const MachineInt& n) const
  {
    RingElem arg(domain(*this), n);
    RingElem ans(codomain(*this));
    mySmartPtr->myApply(raw(ans), raw(arg));
    return ans;
  }

  RingElem RingHom::operator()(const BigInt& N) const
  {
    RingElem arg(domain(*this), N);
    RingElem ans(codomain(*this));
    mySmartPtr->myApply(raw(ans), raw(arg));
    return ans;
  }

  RingElem RingHom::operator()(const BigRat& q) const
  {
    RingElem arg(domain(*this), q);
    RingElem ans(codomain(*this));
    mySmartPtr->myApply(raw(ans), raw(arg));
    return ans;
  }


  RingHom RingHom::operator()(const RingHom& theta) const
  {
    if (codomain(theta) != domain(*this))
      CoCoA_ERROR(ERR::BadCompose, "composition  RingHom(RingHom)");
    return domain(theta)->myCompose(*this, theta);
  }


  PartialRingHom RingHom::operator()(const PartialRingHom& theta) const
  {
    if (codomain(theta) != domain(*this))
      CoCoA_ERROR(ERR::BadCompose, "composition  RingHom(PartialRingHom)");
    return sequential(*this, theta);
  }


//   void RingHomBase::myOutputSelf(std::ostream& out) const
//   {
//     out << "RingHom(" << myDomain << " --> " << myCodomain;
//     myOutputSelfDetails(out);
//     out << ")";
//   }


//   void RingHomBase::myOutputSelfDetails(std::ostream& /*out*/) const
//   {
//     // Default definition does nothing (as there are no extra details to print).
//     // SHOULD THIS BE PURE VIRTUAL (with no default defn)???
//   }


  std::ostream& operator<<(std::ostream& out, const RingHom& phi)
  {
    out << "RingHom(" << domain(phi) << " --> " << codomain(phi);
    phi->myOutputSelfDetails(out);
    out << ")";
    return out;
  }


  // better to use indices or iterators in this fn???
// USE AN ALGORITHM!!!
  bool IsInKer(const ideal& I, const RingHom& phi)
  {
    const vector<RingElem>& g = gens(I);
    const vector<RingElem>::const_iterator stop = g.end();
    for (vector<RingElem>::const_iterator it = g.begin(); it != stop; ++it)
    {
      if (!IsZero(phi(*it))) return false;
    }
    return true;
  }



  //---------------------------------------------------------------------------

  class IdentityRingHomImpl: public RingHomBase
  {
  private:
    explicit IdentityRingHomImpl(const ring& R);
    friend RingHom IdentityHom(const ring& R); // The only function that calls the ctor.
  public:
    virtual void myApply(RingElemRawPtr image, RingElemConstRawPtr arg) const;
    virtual bool IamPartial() const { return false; }
    virtual void myOutputSelfDetails(std::ostream& out) const;
  };


  IdentityRingHomImpl::IdentityRingHomImpl(const ring& R):
      RingHomBase(R, R)
  {}


  void IdentityRingHomImpl::myApply(RingElemRawPtr image, RingElemConstRawPtr arg) const
  {
    myDomain->myAssign(image, arg);
  }


  void IdentityRingHomImpl::myOutputSelfDetails(std::ostream& out) const
  {
    out << " the identity";
  }


  RingHom IdentityHom(const ring& R)
  {
    return RingHom(new IdentityRingHomImpl(R));
  }


  //---------------------------------------------------------------------------

  class SequentialRingHomImpl: public RingHomBase
  {
  private:
    SequentialRingHomImpl(const PartialRingHom& phi, const PartialRingHom& theta);
    friend PartialRingHom sequential(const PartialRingHom& phi, const PartialRingHom& theta); // The only function that calls the ctor.
  public:
    virtual void myApply(RingElemRawPtr image, RingElemConstRawPtr arg) const;
    virtual bool IamPartial() const { return IsPartial(myPhi) || IsPartial(myTheta); }
    virtual void myOutputSelfDetails(std::ostream& out) const;
  private: // data members
    PartialRingHom myPhi;
    PartialRingHom myTheta;
  };


  SequentialRingHomImpl::SequentialRingHomImpl(const PartialRingHom& phi, const PartialRingHom& theta):
      RingHomBase(domain(theta), domain(phi)),
      myPhi(phi),
      myTheta(theta)
  {
    CoCoA_ASSERT(domain(myPhi) == codomain(myPhi));
  }


  void SequentialRingHomImpl::myApply(RingElemRawPtr image, RingElemConstRawPtr arg) const
  {
//    RingElem tmp = myPhi(myTheta(arg));
    RingElem tmp = myPhi(myTheta(RingElemAlias(domain(myTheta),arg)));
    myDomain->myAssign(image, raw(tmp));
  }


  void SequentialRingHomImpl::myOutputSelfDetails(std::ostream& out) const
  {
    out << " composite of " << myPhi << " applied to " << myTheta;
  }


  PartialRingHom sequential(const PartialRingHom& phi, const PartialRingHom& theta)
  {
    if (domain(phi) != codomain(theta))
      CoCoA_ERROR(ERR::BadCompose, "sequential(phi,theta)");
    return PartialRingHom(new SequentialRingHomImpl(phi, theta));
  }


  bool ImageLiesInSubfield(const RingHom& phi)
  {
    if (IsField(domain(phi)) || IsField(codomain(phi))) return true;
    return codomain(phi)->myImageLiesInSubfield(phi);
  }


  //---------------------------------------------------------------------------

  void RingHomEmbeddingBase::myOutputSelfDetails(std::ostream& out) const
  {
    out << " canonical embedding";
  }

  //---------------------------------------------------------------------------

  RingHomInducedBase::RingHomInducedBase(const ring& NewDomain, const RingHom& InducingHom):
      RingHomBase(NewDomain, codomain(InducingHom)),
      myInducingHom(InducingHom)
  {}


  void RingHomInducedBase::myOutputSelfDetails(std::ostream& out) const
  {
    out << " induced by " << myInducingHom;
  }


} // end of namespace CoCoA


// RCS header/log
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/RingHom.C,v 1.8 2013/02/21 14:14:42 abbott Exp $
// $Log: RingHom.C,v $
// Revision 1.8  2013/02/21 14:14:42  abbott
// First attempt at implementing PartialRingHom -- some problems remain!!
//
// Revision 1.7  2011/11/09 14:11:58  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.6  2011/08/24 10:28:49  bigatti
// -- renamed QQ --> BigRat
// -- sorted #include
//
// Revision 1.5  2011/08/14 15:52:16  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.4  2011/02/23 15:00:03  bigatti
// -- RingHom can now be applied also to BigRat
//
// Revision 1.3  2008/12/17 12:11:52  abbott
// Changed type from long to MachineInt in operations which use a machine integer
// in place of a RingElem.  The change is "superficial" but affects many files.
//
// Revision 1.2  2007/10/30 17:14:07  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.5  2007/03/07 14:07:58  bigatti
// -- minor: commented argument names for -Wextra
//
// Revision 1.4  2007/03/05 21:06:07  cocoa
// New names for homomorphism pseudo-ctors: removed the "New" prefix.
//
// Revision 1.3  2007/01/13 14:14:34  cocoa
// Overhaul of RingHom code: it nows uses SmartPtrIRC, and printing is more logical.
// Have not yet updated the documentation.
//
// Revision 1.2  2006/11/08 16:21:59  cocoa
// Structural cleaning of RingHom; many consequential changes.
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.3  2006/03/15 18:09:31  cocoa
// Changed names of member functions which print out their object
// into myOutputSelf -- hope this will appease the Intel C++ compiler.
//
// Revision 1.2  2006/03/12 21:28:33  cocoa
// Major check in after many changes
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.1.1.1  2005/05/03 15:47:31  cocoa
// Imported files
//
// Revision 1.4  2005/04/20 15:40:48  cocoa
// Major change: modified the standard way errors are to be signalled
// (now via a macro which records filename and line number).  Updated
// documentation in error.txt accordingly.
//
// Improved the documentation in matrix.txt (still more work to be done).
//
// Revision 1.3  2005/04/19 14:06:03  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.2  2005/02/11 14:15:20  cocoa
// New style ring elements and references to ring elements;
// I hope I have finally got it right!
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.5  2004/11/12 15:49:29  cocoa
// Tidying prior to 0.90 release.
// (a) documentation improved (or marked as poor)
// (b) sundry minor improvements to the code
//
// Revision 1.4  2004/07/20 15:04:06  cocoa
// The next step in the new "ring element" conversion process:
// handling the case of creating a "const RefRingElem" object
// (since C++ refuses to do this properly itself).
//
// Revision 1.3  2004/05/24 15:52:13  cocoa
// Major update:
//   new error mechanism
//   many fixes
//   RingHoms almost work now
//   RingFloat much improved
//
// Revision 1.2  2003/10/09 12:16:38  cocoa
// New coding convention for rings.
//
// Revision 1.2  2003/06/23 16:59:28  abbott
// Minor cleaning prior to public release.
//
// Revision 1.1  2003/05/14 17:16:40  abbott
// Initial revision
//
