#ifndef CoCoA_RingDistrMPolyClean_H
#define CoCoA_RingDistrMPolyClean_H

//   Copyright (c)  2005-2007,2010-2011  John Abbott

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


#include "CoCoA/DistrMPolyClean.H"

#include <vector>
//using std::vector;

namespace CoCoA
{

  class symbol;         // forward declaration -- defined in symbol.H
  class PPOrdering;     // forward declaration -- defined in PPOrdering.H
  class PPOrderingCtor; // forward declaration -- defined in PPOrdering.H
  class PPMonoid;       // forward declaration -- defined in PPMonoid.H
  class SparsePolyRing; // forward declaration -- defined in SparsePolyRing.H
  class ring;           // forward declaration -- defined in ring.H

  // Pseudo-ctors for (sparse) polynomial rings.

  SparsePolyRing NewPolyRing_DMP(const ring& CoeffRing, const PPMonoid& PPM);
  SparsePolyRing NewPolyRing_DMP(const ring& CoeffRing, const std::vector<symbol>& IndetNames, const PPOrdering& ord);
  SparsePolyRing NewPolyRing_DMP(const ring& CoeffRing, const std::vector<symbol>& IndetNames, const PPOrderingCtor& ord);
  SparsePolyRing NewPolyRing_DMP(const ring& CoeffRing, const std::vector<symbol>& IndetNames); ///< Default ordering is StdDegRevLex
  SparsePolyRing NewPolyRing_DMP(const ring& CoeffRing, long NumIndets); ///< Default ordering is StdDegRevLex
  SparsePolyRing NewPolyRing_DMP(const ring& CoeffRing, long NumIndets, const PPOrderingCtor& ord);

  class RingDistrMPolyCleanImpl: public SparsePolyRingBase
  {
  private:
    typedef DistrMPolyClean value_t; // DistrMPolyClean is the actual type of the values in a RingDistrMPolyCleanImpl
    static value_t& import(RingElemRawPtr rawf);
    static const value_t& import(RingElemConstRawPtr rawf);

  public:
    RingDistrMPolyCleanImpl(const ring& R, const PPMonoid& PPM);
    ~RingDistrMPolyCleanImpl();

  protected: // Data members of RingDistrMPolyCleanImpl
    const ring myCoeffRingValue;  ///< the coefficient ring
    const PPMonoid myPPMValue; ///< the monoid of the power-products
    mutable MemPool myDMPPool; ///< memory manager for polynomials
    mutable MemPool mySummandPool; ///< memory manager for summands; MemPool MUST COME BEFORE myZeroPtr, myOnePtr, and myIndetVector!
    std::unique_ptr<RingElem> myZeroPtr;  ///< Every ring stores its own zero.
    std::unique_ptr<RingElem> myOnePtr;   ///< Every ring stores its own one.
    std::vector<RingElem> myIndetVector;

  public:  // functions which every ring must implement
    virtual ConstRefRingElem myZero() const;
    virtual ConstRefRingElem myOne() const;
    using RingBase::myNew;    // disable warnings of overloading
    using RingBase::myAssign; // disable warnings of overloading
    using PolyRingBase::myIndets; // disable warnings of overloading
    virtual RingElemRawPtr myNew() const;
    virtual RingElemRawPtr myNew(const MachineInt& n) const;
    virtual RingElemRawPtr myNew(const BigInt& N) const;
    virtual RingElemRawPtr myNew(ConstRawPtr rawcopy) const;
    virtual void myDelete(RawPtr rawx) const;                             // destroys x (incl all resources)
    virtual void mySwap(RawPtr rawx, RawPtr rawy) const;                  // swap(x, y)
    virtual void myAssign(RawPtr rawlhs, ConstRawPtr rawx) const;         // lhs = x
    virtual void myAssign(RawPtr rawlhs, const MachineInt& n) const;  // lhs = n
    virtual void myAssign(RawPtr rawlhs, const BigInt& N) const;          // lhs = N
    virtual void myAssignZero(RawPtr rawlhs) const;                       // lhs = 0
    virtual void myRecvTwinFloat(RawPtr rawlhs, ConstRawPtr rawx) const;
    virtual void myNegate(RawPtr rawlhs, ConstRawPtr rawx) const;         // lhs = -x
    virtual void myAdd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const;        // lhs = x+y
    virtual void mySub(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const;        // lhs = x-y
    virtual std::string myImplDetails() const {return "RingDistrMPolyClean";}
    virtual bool myIsZero(ConstRawPtr rawx) const;                                // x == 0
    virtual bool myIsEqual(ConstRawPtr rawx, ConstRawPtr rawy) const;                // x == y

    // functions which every PolyRing must implement

    virtual long myNumIndets() const;
    virtual const ring& myCoeffRing() const;
    virtual const std::vector<RingElem>& myIndets() const;
    virtual void myIndetPower(RawPtr rawf, long var, long exp) const;

    ///@name Simple functions on polynomials
    //@{
    virtual long myNumTerms(ConstRawPtr rawf) const;
    virtual bool myIsMonomial(ConstRawPtr rawf) const;
    virtual RingElemAlias myLC(ConstRawPtr rawf) const;
    virtual void myMulByCoeff(RawPtr rawf, ConstRawPtr rawc) const; ///< WEAK EXCEPTION GUARANTEE
    virtual bool myDivByCoeff(RawPtr rawf, ConstRawPtr rawc) const; ///< WEAK EXCEPTION GUARANTEE
    //@}

    //----------------------------------------------------------------------
    // Functions which every SparsePolyRing must implement:
    //----------------------------------------------------------------------

    virtual const PPMonoid& myPPM() const;

    ///@name   Functions for creating/building polynomials
    //@{
    virtual RingElem myMonomial(ConstRawPtr rawc, PPMonoidElemConstRawPtr rawpp) const;
    virtual SparsePolyIter myBeginIter(ConstRawPtr rawf) const;
    virtual SparsePolyIter myEndIter(ConstRawPtr rawf) const;

    virtual void myPushFront(RawPtr rawf, ConstRawPtr rawc, const std::vector<long>& expv) const;
    virtual void myPushBack(RawPtr rawf, ConstRawPtr rawc, const std::vector<long>& expv) const;
    virtual void myPushFront(RawPtr rawf, ConstRawPtr rawc, PPMonoidElemConstRawPtr rawpp) const;
    virtual void myPushBack(RawPtr rawf, ConstRawPtr rawc, PPMonoidElemConstRawPtr rawpp) const;
    //@}

    ///@name Special functions on polynomials needed for implementing Buchberger's Algorithm
    //@{
    virtual ConstRefPPMonoidElem myLPP(ConstRawPtr rawf) const;
    virtual void myMulByPP(RawPtr rawf, PPMonoidElemConstRawPtr rawpp) const;
    virtual bool myIsZeroAddLCs(RawPtr rawf, RawPtr rawg) const; ///< f+=LM(g); g-=LM(g); assumes LPP(f)==LPP(g); returns LC(f)+LC(g)==0
    virtual void myMoveLM(RawPtr rawf, RawPtr rawg) const; ///< f+=LM(g); g-=LM(g); assumes LM(f)<LM(g)
    virtual void myDeleteLM(RawPtr rawf) const; // ????? right interface
    virtual void myDivLM(RawPtr rawlhs, ConstRawPtr rawf, ConstRawPtr rawg) const; ///< lhs=div(LM(f),LM(g)); assumes f!=0,g!=0
    virtual int  myCmpLPP(ConstRawPtr rawf, ConstRawPtr rawg) const; ///< cmp(LPP(f),LPP(g)); assumes f!=0,g!=0
    virtual void myAddClear(RawPtr rawf, RawPtr rawg) const; ///< f+=g; g=0;
    virtual void myAppendClear(RawPtr rawf, RawPtr rawg) const; ///< f+=g; g=0; appends g to f with no checks

    virtual void myAddMul(RawPtr rawf, ConstRawPtr rawh, ConstRawPtr rawg) const; ///<  f += LM(h)*g
    virtual void myAddMul(RawPtr rawf, ConstRawPtr rawh, ConstRawPtr rawg, SkipLMFlag) const; ///<  f += LM(h)*g
    virtual void myReductionStep(RawPtr rawf, ConstRawPtr rawg) const;
    // ??? aggiungere coefficiente
    virtual void myReductionStepGCD(RawPtr rawf, ConstRawPtr rawg, RingElem& FScale) const;
    // should it all be in ReductionStep ??? ANNA
    //@}

  };

} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/RingDistrMPolyClean.H,v 1.2 2011/03/10 16:39:35 abbott Exp $
// $Log: RingDistrMPolyClean.H,v $
// Revision 1.2  2011/03/10 16:39:35  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.1  2010/10/08 08:05:28  bigatti
// -- renamed (Ring)DistrMPoly --> (Ring)DistrMPolyClean
//
// Revision 1.3  2010/03/05 18:43:48  abbott
// Added pseudo-ctors allowing polynomial rings to be created specifying
// the ordering using a PPOrderingCtor object.
//
// Revision 1.2  2007/10/30 17:14:12  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.8  2007/03/08 11:07:12  cocoa
// Made pseudo ctors for polynomial rings more uniform.  This allowed me to
// remove an include of CoCoA/symbol.H  from the RingDistrM*.H files, but then
// I had to put the include in several .C files.
//
// Revision 1.7  2007/02/10 18:44:03  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.6  2006/12/06 17:15:26  cocoa
// -- removed #include "config.H"
//
// Revision 1.5  2006/11/24 17:25:00  cocoa
// -- reorganized includes of header files
//
// Revision 1.4  2006/11/22 14:46:26  cocoa
// -- added #include "CoCoA/symbol.H"
//
// Revision 1.3  2006/10/06 14:04:15  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.2  2006/06/08 16:45:28  cocoa
// -- RingDistrMPolyClean*.H  have been "moved" into RingDistrMPolyClean*.C
// -- some coding conventions fixed in DistrMPolyClean*
// -- functions wdeg and CmpWDeg have a common implementation in SparsePolyRing
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.13  2006/05/29 16:22:37  cocoa
// Third time lucky???
// Added myIsInteger member function to all rings (NYI for RingFloat).
//
// Revision 1.12  2006/05/12 17:06:09  cocoa
// -- moved myIsUnit, myGcd into SparsePolyRing (common implementation)
//
// Revision 1.11  2006/05/12 16:10:58  cocoa
// Added OpenMathFwd.H, and tidied OpenMath.H.
// Many consequential but trivial changes.
//
// Revision 1.10  2006/04/26 16:44:53  cocoa
// -- myMul has now a single implementation in SparsePolyRing
// -- myMul and mul in RingDistrMPolyClean* and DistrMPolyClean* have been disabled
//
// Revision 1.9  2006/03/21 09:43:14  cocoa
// Changed names of some member fns of ideals (dealing with setting and testing
// the flags for primeness and maximality).  Hope icc will complain less now.
//
// Revision 1.8  2006/03/17 18:10:27  cocoa
// -- changed: myMul --> myMulByPP
//
// Revision 1.7  2006/03/15 18:09:31  cocoa
// Changed names of member functions which print out their object
// into myOutputSelf -- hope this will appease the Intel C++ compiler.
//
// Revision 1.6  2006/03/14 15:01:49  cocoa
// Improved the implementation of ring member fns for computing powers.
// Should keep Intel C++ compiler quieter too.
//
// Revision 1.5  2006/03/12 21:28:34  cocoa
// Major check in after many changes
//
// Revision 1.4  2006/03/07 10:06:12  cocoa
// -- fixed: PPMonoidElem LPP(f) now returns ConstRefPPMonoidElem
//
// Revision 1.3  2006/02/20 22:41:20  cocoa
// All forms of the log function for power products now return SmallExponent_t
// (instead of int).  exponents now resizes the vector rather than requiring
// the user to pass in the correct size.
//
// Revision 1.2  2006/01/19 18:03:53  cocoa
// -- fixed coding conventions for myGBasis, myGBasisValue
// -- fixed myReduceMod (NF - tested)
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.6  2005/08/08 16:36:33  cocoa
// Just checking in before going on holiday.
// Don't really recall what changes have been made.
// Added IsIndet function for RingElem, PPMonoidElem,
// and a member function of OrdvArith.
// Improved the way failed assertions are handled.
//
// Revision 1.5  2005/07/15 16:34:33  cocoa
// Added iterators for sparse polynomials.
// The code compiles (and the old tests still run).
// It'd Friday evening -- I'm going home before
// getting any ideas about making the iterator code run.
//
// Revision 1.4  2005/07/08 15:09:29  cocoa
// Added new symbol class (to represent names of indets).
// Integrated the new class into concrete polynomial rings
// and PPMonoid -- many consequential changes.
// Change ctors for the "inline" sparse poly rings: they no
// longer expect a PPMonoid, but build their own instead
// (has to be a PPMonoidOv).
//
// Revision 1.3  2005/07/01 16:08:16  cocoa
// Friday check-in.  Major change to structure under PolyRing:
// now SparsePolyRing and DUPolyRing are separated (in preparation
// for implementing iterators).
//
// A number of other relatively minor changes had to be chased through
// (e.g. IndetPower).
//
// Revision 1.2  2005/06/22 14:47:56  cocoa
// PPMonoids and PPMonoidElems updated to mirror the structure
// used for rings and RingElems.  Many consequential changes.
//
// Revision 1.1.1.1  2005/05/03 15:47:30  cocoa
// Imported files
//
// Revision 1.5  2005/04/27 16:14:56  cocoa
// Cleaned up example programs -- added "free use" permit.
// Changed a couple of ErrorInfo object names, and added
// ERR::NotTrueGCDDomain.
//
// Revision 1.4  2005/04/21 15:12:19  cocoa
// Revised NewPolyRing as Dag Arneson suggested (perhaps just an interim
// measure).
// Brought example programs up to date (new name for CoCoA error
// information objects).
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
// Revision 1.12  2004/11/18 18:33:41  cocoa
// Now every ring know its own "one" element (as well as "zero").
// Several consequential changes.
//
// Revision 1.11  2004/11/11 13:11:09  cocoa
// -- minor changes for doxygen
//
// Revision 1.10  2004/11/04 18:47:43  cocoa
// (1) Ring member functions which previously expected mpz_t args
//     now expect ZZ args.  Numerous minor consequential changes.
// (2) Renamed function which gives access to the mpz_t value inside
//     a ZZ object: previously was raw(...), now is mpzref(...).
//     Plenty of calls had to be altered.
//
// Revision 1.9  2004/10/21 17:16:37  cocoa
// Fairly major change: new OrdvArith namspace with various members,
//   new global typedef  SmallExponent_t (defined in config.H).
//
// Revision 1.8  2004/07/27 16:03:39  cocoa
// Added IsCommutative test and IamCommutative member function
// to all rings.  Tidied geobuckets a little.
//
// Revision 1.7  2004/07/20 15:04:06  cocoa
// The next step in the new "ring element" conversion process:
// handling the case of creating a "const RefRingElem" object
// (since C++ refuses to do this properly itself).
//
// Revision 1.6  2004/05/24 15:52:14  cocoa
// Major update:
//   new error mechanism
//   many fixes
//   RingHoms almost work now
//   RingFloat much improved
//
// Revision 1.5  2004/04/08 15:33:34  cocoa
// Added function IsInteger, and the related RingBase::myIsInteger
// virtual function, plus all necessary implementations.
//
// Revision 1.4  2004/03/20 17:46:10  cocoa
// Check in prior to departure to RWCA
//
// Revision 1.3  2004/01/30 14:07:10  cocoa
// Tidied RingRawValue union: now it contains just two fields,
// and has no need of forward declarations of types used internally
// by the concrete rings -- it uses explicitly a void* instead.
//
// I have tidied the "import" functions used by most concrete rings.
//
// I have moved the choice of representation type for RingFp and RingFpLog
// into a typedef in config.H -- this is to recognise that different
// choices may work best on different platforms.
//
// Revision 1.2  2004/01/28 16:27:00  cocoa
// Added the necessary for CmpDeg to work.
//
// Revision 1.1  2003/11/21 14:33:54  cocoa
// -- First Import
//

#endif
