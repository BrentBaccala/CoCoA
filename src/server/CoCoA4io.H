
#ifndef CoCoA_CoCoA4io_H
#define CoCoA_CoCoA4io_H

//   Copyright (c)  2007  Anna Bigatti

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

//#include "CoCoA/library.H"

#include "CoCoA/TmpGTypes.H"  // PolyList and VectorList

#include <string>

namespace CoCoA
{
  class BigInt;           // forward declaration -- CoCoA/BigInt.H
  class RingElemAlias; // forward declaration -- CoCoA/ring.H
  typedef const RingElemAlias& ConstRefRingElem; // forward declaration -- CoCoA/ring.H
  class FreeModule;       // forward declaration -- CoCoA/FreeModule.H
  class PPMonoid;         // forward declaration -- CoCoA/PPMonoid.H
  class PPMonoidElem;     // forward declaration -- CoCoA/PPMonoid.H
  class PPOrdering;       // forward declaration -- CoCoA/PPOrdering.H
  class SparsePolyRing;   // forward declaration -- CoCoA/SparsePolyRing.H
  class ideal;            // forward declaration -- CoCoA/ideal.H
  class matrix;           // forward declaration -- CoCoA/matrix.H
  //----------------------------------------------------------------------
  enum SkipTagType {TagWasRead, GetTag};  // used in functions Read...(...)

  //----------------------------------------------------------------------
  std::string ReadOperationString(std::istream& in);
  void ReadPPs(std::istream& in, std::vector<PPMonoidElem>& PPs, const PPMonoid& PPM, SkipTagType ST);
  void ReadPolyList(std::istream& in, PolyList& PL, const SparsePolyRing& P, SkipTagType ST);
  void ReadVectorList(std::istream& in, VectorList& L, const FreeModule& FM, SkipTagType ST);

  // cannot be "void Read..(res, ...)" because there is no empty ctor
  // anyway, they are reference counted objects
  int ReadVerbosityLevel(std::istream& in, SkipTagType ST);
  matrix ReadIntegerMatrix(std::istream& in, SkipTagType ST);
  matrix ReadRationalMatrix(std::istream& in, SkipTagType ST);
  FreeModule ReadFreeModule(std::istream& in, SkipTagType ST);
  SparsePolyRing ReadPolyRing(std::istream& in, SkipTagType ST);
  PPOrdering ReadPolyRingOrdering(std::istream& in, long NumInds, SkipTagType ST);

  //------------------------------
  void PrintTimeToLog(double T);
  void PrintTimeToCoCoA4(double T);
  void PrintVersionToCoCoA4();
  void EndOfTransmissionToCoCoA4();

  void WritePolyRing(std::ostream& out, const SparsePolyRing& P);
  void WriteMonomials(std::ostream& out, ConstRefRingElem p);
  void WriteMatrix(std::ostream& out, const matrix &M);
  void WritePoly(std::ostream& out, ConstRefRingElem f);
  void WritePolyList(std::ostream& out, const PolyList &PL);
  void WriteIdeal(std::ostream& out, const ideal &I, const std::string& ClosingStr);

  void WriteMatrixInVar(std::ostream& out, const std::string& VarName, const matrix &M);
  void WritePolyListInVar(std::ostream& out, const std::string& VarName, const PolyList& PL);
  void WriteVectorListInVar(std::ostream& out, const std::string& VarName, const VectorList& VL);
  //----------------------------------------------------------------------
  void AssertTag(const std::string& input_tag, const std::string& expected_tag);
  void SkipTag(std::istream& in, const std::string& expected_tag);
  void ThrowInputError(const std::string& unknown_tag);

  //----------------------------------------------------------------------
  SparsePolyRing NewPolyRingServer(const BigInt& charact, long FloatPrecision, long NumParams, const PPOrdering& O);
}


// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/server/CoCoA4io.H,v 1.1 2013/05/27 12:57:39 abbott Exp $
// $Log: CoCoA4io.H,v $
// Revision 1.1  2013/05/27 12:57:39  abbott
// Moved all server-related code into src/server/
//
// Revision 1.13  2012/10/24 11:38:30  abbott
// Revised fwd decl of ConstRefRingElem.
//
// Revision 1.12  2011/08/14 15:52:18  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.11  2011/03/10 17:04:00  bigatti
// -- using long instead of size_t
//
// Revision 1.10  2009/10/29 19:12:59  abbott
// Added EndOfTransmissionToCoCoA4 function.
// Used this new function in CoCoAServer to implement "over-and-out" handshake
// at the end of a CoCoA4-CoCoAServer session.
// New version number to mark this change -- not backward compatible!
//
// Revision 1.9  2009/10/26 17:40:01  bigatti
// -- added missing std::
//
// Revision 1.8  2009/10/26 17:20:32  bigatti
// -- trying to avoid sending small strings to sockets...
//
// Revision 1.7  2009/05/20 14:31:18  abbott
// Changed first param of NewPolyRingServer into a const reference
// to ZZ, rather than pass by value.
//
// Revision 1.6  2009/01/26 15:55:23  bigatti
// -- minor cleanup
// -- added WriteIdeal, WriteMatrix, WritePoly, WritePolyList
// -- sorted included files
//
// Revision 1.5  2008/09/19 11:34:16  bigatti
// -- new mechanism for passing verbosity level (or StatLevel)
//    [only partially tested]
//
// Revision 1.4  2008/05/28 16:19:02  bigatti
// -- added ReadPPs for reading a list of power-products in PPMonoidEv
//    (fo exponent overflow in PPMonoidOv)
//
// Revision 1.3  2008/04/18 15:35:57  abbott
// (long overdue) Major revision to matrices
//
// Revision 1.2  2007/10/30 17:14:12  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1  2007/04/27 14:54:22  bigatti
// -- content of CoCoAServer.C split into dedicated files
// -- new registration mechanism (through include "RegisterServerOps.H")
//

#endif
